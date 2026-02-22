#!/usr/bin/env python3
"""
MitoReach 2.0 Standalone: State-of-the-Art Mitochondrial DNA Editing Prediction
Physics: Joint Dimerization Probability + Boltzmann Linker Elasticity + Volumetric Sterics
Features: SASD (70Å Cutoff) + Approach Angle (90° Peak) + Automated Spacer Sweep

This is a unified single-file version of the MitoReach 2.0 suite, ready for sharing.
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import json
import csv
import heapq
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass
from Bio.PDB import Structure, Atom, Polypeptide

# ==============================================================================
# CORE CONSTANTS & SETTINGS
# ==============================================================================
VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
    'P': 1.80, 'F': 1.47, 'CL': 1.75, 'BR': 1.85, 'I': 1.98,
    'SE': 1.90, 'ZN': 1.39, 'MG': 1.73, 'CA': 2.31, 'FE': 1.47
}

# ==============================================================================
# 1. SASD CALCULATOR ENGINE (Solvent Accessible Surface Distance)
# ==============================================================================

@dataclass
class GridCell:
    indices: Tuple[int, int, int]
    center: np.ndarray
    occupied: bool = False
    distance: float = float('inf')
    visited: bool = False
    
    def __lt__(self, other):
        return self.distance < other.distance

class Grid3D:
    def __init__(self, min_coords: np.ndarray, max_coords: np.ndarray, grid_spacing: float = 1.0):
        self.grid_spacing = grid_spacing
        self.min_coords = min_coords
        self.max_coords = max_coords
        self.dimensions = np.ceil((max_coords - min_coords) / grid_spacing).astype(int) + 1
        self.cells: Dict[Tuple[int, int, int], GridCell] = {}
        
    def get_lazy(self, i, j, k) -> GridCell:
        idx = (i, j, k)
        if idx not in self.cells:
            # We don't allocate np.array unless needed
            center = self.min_coords + np.array([i, j, k]) * self.grid_spacing
            self.cells[idx] = GridCell(indices=idx, center=center)
        return self.cells[idx]
        
    def get_cell(self, coords: np.ndarray) -> Optional[GridCell]:
        indices = np.floor((coords - self.min_coords) / self.grid_spacing).astype(int)
        # Check bounds before creating lazy cell to avoid infinite expanses
        if (0 <= indices[0] < self.dimensions[0] and
            0 <= indices[1] < self.dimensions[1] and
            0 <= indices[2] < self.dimensions[2]):
            return self.get_lazy(*indices)
        return None
    
    def get_neighbors(self, cell: GridCell) -> List[GridCell]:
        neighbors = []
        i, j, k = cell.indices
        # 26-Connectivity: Includes diagonals to prevent distance inflation
        for di in [-1,0,1]:
            for dj in [-1,0,1]:
                for dk in [-1,0,1]:
                    if di == 0 and dj == 0 and dk == 0: continue
                    ni, nj, nk = i + di, j + dj, k + dk
                    
                    if (0 <= ni < self.dimensions[0] and
                        0 <= nj < self.dimensions[1] and
                        0 <= nk < self.dimensions[2]):
                        n_cell = self.get_lazy(ni, nj, nk)
                        neighbors.append(n_cell)
        return neighbors

class SASDCalculator:
    def __init__(self, grid_spacing: float = 1.0, probe_radius: float = 1.4):
        self.grid_spacing = grid_spacing
        self.probe_radius = probe_radius

    def calculate_sasd(self, structure, atom1_pos, atom2_pos, max_distance=75.0) -> float:
        all_atoms = []
        # OPTIMIZATION: Only include atoms within reach of the path
        path_mid = (atom1_pos + atom2_pos) / 2
        path_radius = np.linalg.norm(atom2_pos - atom1_pos) / 2 + 15.0  # 15Å margin
        
        for model in structure:
            for chain in model:
                for res in chain:
                    for atom in res:
                        coord = atom.get_coord()
                        if np.linalg.norm(coord - path_mid) <= path_radius:
                            elem = atom.element.strip().upper() if atom.element else 'C'
                            all_atoms.append((coord, VDW_RADII.get(elem, 1.7)))
        
        if not all_atoms: return float('inf') # Should not happen with atoms of interest
        
        coords = np.array([a[0] for a in all_atoms])
        # Grid Padding: 10A for solvent navigation
        min_c, max_c = coords.min(axis=0) - 10, coords.max(axis=0) + 10
        grid = Grid3D(min_c, max_c, self.grid_spacing)
        
        for p, r in all_atoms:
            rad = r + self.probe_radius
            rad_sq = rad * rad
            px, py, pz = p
            min_idx = np.floor((p - rad - grid.min_coords) / grid.grid_spacing).astype(int)
            max_idx = np.ceil((p + rad - grid.min_coords) / grid.grid_spacing).astype(int)
            
            i_start = max(0, min_idx[0])
            i_end = min(grid.dimensions[0], max_idx[0]+1)
            j_start = max(0, min_idx[1])
            j_end = min(grid.dimensions[1], max_idx[1]+1)
            k_start = max(0, min_idx[2])
            k_end = min(grid.dimensions[2], max_idx[2]+1)
            
            for i in range(i_start, i_end):
                for j in range(j_start, j_end):
                    for k in range(k_start, k_end):
                        c = grid.get_lazy(i, j, k)
                        if not c.occupied:
                            cx, cy, cz = c.center
                            if (cx-px)**2 + (cy-py)**2 + (cz-pz)**2 <= rad_sq:
                                c.occupied = True
        
        source, target = grid.get_cell(atom1_pos), grid.get_cell(atom2_pos)
        if not source or not target: return float('inf')
        
        # FIX: If source/target are buried inside protein, find nearest unoccupied neighbor.
        # This handles terminal residue CA atoms that are partially buried.
        def nearest_unoccupied(grid, cell):
            if not cell.occupied:
                return cell
            # BFS over expanding shells to find nearest unoccupied cell
            from collections import deque
            seen = {cell.indices}
            q = deque([cell])
            while q:
                c = q.popleft()
                for nb in grid.get_neighbors(c):
                    if nb.indices in seen:
                        continue
                    seen.add(nb.indices)
                    if not nb.occupied:
                        return nb
                    q.append(nb)
            return None
        
        source = nearest_unoccupied(grid, source)
        target = nearest_unoccupied(grid, target)
        if not source or not target: return float('inf')
        
        pq = [(0.0, source)]
        source.distance = 0.0
        
        # Precompute constants for performance
        dist_1 = self.grid_spacing
        dist_2 = self.grid_spacing * 1.41421356
        dist_3 = self.grid_spacing * 1.73205081
        
        while pq:
            d, curr = heapq.heappop(pq)
            if curr.indices == target.indices: return d
            if d > curr.distance: continue
            
            for neighbor in grid.get_neighbors(curr):
                if neighbor.occupied: continue
                
                # Fast distance calculation
                dx = abs(curr.indices[0] - neighbor.indices[0])
                dy = abs(curr.indices[1] - neighbor.indices[1])
                dz = abs(curr.indices[2] - neighbor.indices[2])
                ds = dx + dy + dz
                
                if ds == 1: dist_inc = dist_1
                elif ds == 2: dist_inc = dist_2
                else: dist_inc = dist_3
                
                if d + dist_inc > max_distance: continue
                new_d = d + dist_inc
                
                if new_d < neighbor.distance:
                    neighbor.distance = new_d
                    heapq.heappush(pq, (new_d, neighbor))
        return float('inf')

# ==============================================================================
# 2. GEOMETRIC REACH ENGINE (MitoReach 2.0 Scoring Logic)
# ==============================================================================

class GeometricReachEngine:
    def __init__(self, pdb_template):
        if pdb_template.endswith('.cif'):
            from Bio.PDB.MMCIFParser import MMCIFParser
            parser = MMCIFParser(QUIET=True)
        else:
            from Bio.PDB.PDBParser import PDBParser
            parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure('template', pdb_template)
        
    def compute_favorability_score(self, sasd, angle, obstructions, variant='DddA11'):
        if sasd == float('inf'): return 0.0
        
        # 1. Angle Score (90 deg peak)
        if angle >= 70: angle_score = 1.0
        elif angle <= 35: angle_score = 0.1
        else: angle_score = 0.1 + 0.9 * (angle - 35) / (70 - 35)
        
        # 2. Obstruction Penalty
        obs_penalty = max(0.2, 1.0 - (obstructions * 0.25))
        
        # 3. SASD Decay (70A Hard Cutoff)
        if sasd <= 45: sasd_score = 1.0
        elif sasd >= 70: sasd_score = 0.0
        else: sasd_score = 1.0 - (sasd - 45) / (70 - 45)
        
        if variant == 'DddA11':
            return (sasd_score * 0.4 + angle_score * 0.4 + obs_penalty * 0.2) * 100
        return (sasd_score * 0.3 * angle_score * obs_penalty) * 100

# ==============================================================================
# 3. MITOREACH 2.0 SYSTEM (Orchestrator)
# ==============================================================================

class MitoReach2System:
    def __init__(self, template_pdb):
        self.engine = GeometricReachEngine(template_pdb)
        self.comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
        # Store original Chain A coordinates HERE (once at init) to prevent stale sweep data
        self._original_coords = {
            atom: atom.coord.copy() 
            for atom in self.engine.structure.get_atoms() 
            if atom.get_parent().get_parent().id == 'A'
        }
        
    def get_favorability_profile(self, spacer_seq, variant='DddA11'):
        results = []
        spacer_len = len(spacer_seq)
        spacer_ts = "".join([self.comp.get(b, 'N') for b in spacer_seq])
        
        # Helical constants (B-DNA)
        h_step, rot_deg = 3.38, 34.3
        
        # Helical Axis for 9jo8 (Experimental Reference)
        C_AXIS = np.array([168.962, 167.737, 181.423])
        V_AXIS = np.array([-0.019, -0.234, 0.972])
        V_AXIS = V_AXIS / np.linalg.norm(V_AXIS)
        
        # WLC Physics (16aa linker: sigma_sq = 2 * lp * L = 425.6 A^2)
        # RMS end-to-end distance ~20.6A. 
        # Total optimal reach (Linker + DddA geometry) ≈ 35.0A based on 9jo8.
        sigma_sq = 425.6 
        l_optimal = 35.0 

        calculator = SASDCalculator(grid_spacing=1.5) # Optimized for sweep

        for pos in range(1, spacer_len + 1):
            # 1. Structural Simulation: Rotate and Translate TALE
            # shift = 0 is original 9jo8 position (P=4)
            shift = pos - 4
            angle_deg = shift * rot_deg
            translation = shift * h_step * V_AXIS
            
            # Reset to original baseline and apply transform
            for atom, coord in self._original_coords.items():
                atom.coord = coord.copy()
            self._apply_helical_transform(self.engine.structure, C_AXIS, V_AXIS, angle_deg, translation)
            
            # 2. Real SASD Calculation
            # Anchor: TALE C-term (641 CA)
            # Targets: DddA N-terms (1290 and 1398)
            try:
                # Find CA atoms for anchor and targets
                pos_a = self._get_ca_pos(self.engine.structure, 'A', 641)
                pos_dn = self._get_ca_pos(self.engine.structure, 'D', 1290)
                pos_dc = self._get_ca_pos(self.engine.structure, 'D', 1398)
                
                sasd_n = calculator.calculate_sasd(self.engine.structure, pos_a, pos_dn)
                sasd_c = calculator.calculate_sasd(self.engine.structure, pos_a, pos_dc)
            except Exception as e:
                # print(f"Warning: SASD fail at P={pos}: {e}")
                sasd_n, sasd_c = float('inf'), float('inf')

            # 3. Joint Polymer Elasticity (Dimerization Probability)
            # Uses SASD as the constraint distance
            p_reach_n = np.exp(-max(0, sasd_n - l_optimal)**2 / sigma_sq) if sasd_n != float('inf') else 0
            p_reach_c = np.exp(-max(0, sasd_c - l_optimal)**2 / sigma_sq) if sasd_c != float('inf') else 0
            p_dimer = p_reach_n * p_reach_c
            
            # 4. Structural Accessibility (Angle & Sterics)
            # 9jo8-calibrated major groove access.
            angle = max(20.0, 90.0 - abs(pos - (spacer_len+1)/2.0) * 10.0)
            
            # Restore Steric Factor: Positions 1-2 and end-2 are blocked by TALE
            steric_factor = 1.0
            if pos <= 2 or pos >= spacer_len - 1:
                steric_factor = 0.2 # Significant penalty for steric zones
            
            obs = 1 if steric_factor < 1.0 else 0

            # 5. Score Calculations
            # N-term fragment (Anchor Left) -> Bottom Strand (TS) Access
            # C-term fragment (Anchor Right) -> Top Strand (NTS) Access
            score_t = self.engine.compute_favorability_score(sasd_n, angle, obs, variant)
            score_t *= (p_dimer * steric_factor)
            if spacer_ts[pos-1] != 'C': 
                score_t = 0.0
            elif variant == 'WT':
                # TS constraint (read 5' to 3', so pos-2 is preceding base)
                if pos == 1 or spacer_ts[pos-2] != 'T':
                    score_t = 0.0
            
            score_n = self.engine.compute_favorability_score(sasd_c, angle, obs, variant)
            score_n *= (p_dimer * steric_factor)
            if spacer_seq[pos-1] != 'C': 
                score_n = 0.0
            elif variant == 'WT':
                # NTS constraint
                if pos == 1 or spacer_seq[pos-2] != 'T':
                    score_n = 0.0
            
            results.append({
                'pos': pos, 'score_n': score_n, 'score_t': score_t,
                'sasd': min(sasd_n, sasd_c), 'angle': angle
            })
        return results

    def _apply_helical_transform(self, structure, center, axis, angle_deg, translation):
        """Rotate and translate Chain A of the structure."""
        def rotation_matrix(axis, angle_deg):
            angle_rad = np.radians(angle_deg)
            ux, uy, uz = axis
            c, s = np.cos(angle_rad), np.sin(angle_rad)
            return np.array([
                [c + ux**2 * (1-c), ux*uy*(1-c) - uz*s, ux*uz*(1-c) + uy*s],
                [uy*ux*(1-c) + uz*s, c + uy**2 * (1-c), uy*uz*(1-c) - ux*s],
                [uz*ux*(1-c) - uy*s, uz*uy*(1-c) + ux*s, c + uz**2 * (1-c)]
            ])
        
        R = rotation_matrix(axis, angle_deg)
        for model in structure:
            for chain in model:
                if chain.id == 'A':
                    for residue in chain:
                        for atom in residue:
                            rel_pos = atom.coord - center
                            rot_pos = np.dot(R, rel_pos)
                            atom.coord = rot_pos + center + translation

    def _get_ca_pos(self, structure, chain_id, res_id) -> np.ndarray:
        """Utility to get CA position for a specific residue."""
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                if res_id in chain:
                    res = chain[res_id]
                    if 'CA' in res:
                        return res['CA'].get_coord()
        raise ValueError(f"Could not find {chain_id}:{res_id}:CA")

    def generate_heatmap(self, heatmap_data, output_file):
        lengths = sorted(heatmap_data.keys())
        max_len = max(lengths)
        grid = np.zeros((len(lengths), max_len))
        
        for i, l in enumerate(lengths):
            for tp in range(1, l + 1):
                grid[i, tp-1] = heatmap_data[l].get(tp, 0)
        
        plt.figure(figsize=(12, 8))
        im = plt.imshow(grid, aspect='auto', cmap='RdYlGn', origin='lower')
        plt.colorbar(im, label='Favorability Score (0-100)')
        plt.yticks(range(len(lengths)), [f"{l} bp" for l in lengths])
        plt.xticks(range(max_len), range(1, max_len + 1))
        plt.xlabel('Target Position in Spacer')
        plt.ylabel('Spacer Length')
        plt.title('MitoReach 2.0: Physical Favorability Landscape')
        
        for i, l in enumerate(lengths):
            for tp in range(1, l + 1):
                val = int(heatmap_data[l].get(tp, 0))
                if val > 0:
                    plt.text(tp-1, i, str(val), ha='center', va='center', color='white' if val > 70 or val < 30 else 'black')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300)

    def plot_results(self, results, output_file, spacer_seq):
        pos = [r['pos'] for r in results]
        score_n = [r['score_n'] for r in results]
        score_t = [r['score_t'] for r in results]
        sasd = [r['sasd'] for r in results]
        angles = [r['angle'] for r in results]
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True, gridspec_kw={'height_ratios': [1, 1.2]})
        ax1.plot(pos, sasd, 'o-', color='black', linewidth=2, label='Surface Distance (SASD)')
        ax1.set_ylabel('SASD (Å)', fontsize=12, fontweight='bold')
        ax_angle = ax1.twinx()
        ax_angle.plot(pos, angles, 's--', color='darkorange', alpha=0.7, label='Approach Angle')
        ax_angle.set_ylabel('Approach Angle (°)', color='darkorange', fontsize=12, fontweight='bold')
        ax_angle.set_ylim(0, 100)
        ax1.set_title('MitoReach 2.0: Structural Metrics', fontsize=14, loc='left')
        ax1.legend(loc='upper right')

        width = 0.35
        ax2.bar(np.array(pos) - width/2, score_n, width=width, label='Top Strand (NTS)', color='#2c3e50')
        ax2.bar(np.array(pos) + width/2, score_t, width=width, label='Bottom Strand (TS)', color='#3498db')
        ax2.set_ylabel('Favorability Score (0-100)', fontsize=12, fontweight='bold')
        ax2.set_xlabel('Spacer Position (bp)', fontsize=12, fontweight='bold')
        ax2.set_ylim(0, 110)
        ax2.grid(axis='y', linestyle='--', alpha=0.6)
        
        spacer_ts = "".join([self.comp.get(b, 'N') for b in spacer_seq])
        xtick_labels = [f"{p}\n{n}/{t}" for p, n, t in zip(pos, spacer_seq, spacer_ts)]
        ax2.set_xticks(pos)
        ax2.set_xticklabels(xtick_labels, fontsize=10)
        ax2.set_title('Predicted Editing Likelihood', fontsize=14, loc='left')
        ax2.legend(loc='upper right')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300)

# ==============================================================================
# MAIN ENTRY POINT
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description="MitoReach 2.0: Professional Standalone Suite")
    parser.add_argument("--sequence", required=True, help="Full DNA sequence containing the target")
    parser.add_argument("--target-pos", type=int, required=True, help="Index of the target C in the sequence (0-indexed)")
    parser.add_argument("--variant", default="DddA11", choices=["WT", "DddA11"])
    parser.add_argument("--template", required=True, help="Path to structural template (PDB/CIF)")
    parser.add_argument("--output-dir", default="mitoreach_results")
    parser.add_argument("--min-len", type=int, default=10, help="Minimum spacer length")
    parser.add_argument("--max-len", type=int, default=18, help="Maximum spacer length")
    args = parser.parse_args()

    if not os.path.exists(args.template):
        print(f"Error: Template {args.template} not found.")
        return

    os.makedirs(args.output_dir, exist_ok=True)
    target_dir = os.path.join(args.output_dir, f"Target_{args.target_pos}")
    os.makedirs(target_dir, exist_ok=True)
    
    system = MitoReach2System(args.template)
    
    # 1. Spacer Length Sweep (10-18 bp)
    sweep_results = {}
    heatmap_data = {}
    
    print(f"Starting MitoReach 2.0 Physical Sweep for sequence index {args.target_pos}...")
    print(f"Data organized into: {target_dir}")
    
    for l in range(args.min_len, args.max_len + 1):
        heatmap_data[l] = {}
        sweep_results[l] = []
        for tp in range(1, l + 1):
            s_start = args.target_pos - (tp - 1)
            s_end = s_start + l
            if s_start >= 0 and s_end <= len(args.sequence):
                spacer_seq = args.sequence[s_start:s_end]
                print(f"  Analysing L={l}bp, Target Pos P={tp}...", end="\r")
                profile = system.get_favorability_profile(spacer_seq, args.variant)
                
                # Best for heatmap
                score = max(profile[tp-1]['score_n'], profile[tp-1]['score_t'])
                heatmap_data[l][tp] = score
                
                bystanders = sum(1 for i, p in enumerate(profile) 
                               if (i+1) != tp and (p['score_n'] > 40 or p['score_t'] > 40))
                
                sweep_results[l].append({
                    'spacer': spacer_seq, 'tp': tp, 'score': score, 
                    'bystanders': bystanders, 'profile': profile
                })

    # 2. Results Generation
    heatmap_png = os.path.join(target_dir, "structural_landscape_heatmap.png")
    system.generate_heatmap(heatmap_data, heatmap_png)
    
    # Generate profiles for transparency
    best_overall = None
    for l in sorted(sweep_results.keys()):
        if not sweep_results[l]: continue
        best_for_len = sorted(sweep_results[l], key=lambda x: (-x['score'], x['bystanders']))[0]
        if best_overall is None or (best_for_len['score'] > best_overall['score']):
            best_overall = best_for_len
            
        profile_png = os.path.join(target_dir, f"profile_L{l}.png")
        system.plot_results(best_for_len['profile'], profile_png, best_for_len['spacer'])
        print(f" - Profile for {l}bp saved: {profile_png}")

    print("\n--- MitoReach 2.0 Optimal Discovery ---")
    print(f"Best Spacer: {best_overall['spacer']} ({len(best_overall['spacer'])} bp)")
    print(f"Target Pos: {best_overall['tp']} | Peak Score: {best_overall['score']:.1f} | Bystanders: {best_overall['bystanders']}")
    print(f"Results can be found in: {target_dir}/")

if __name__ == "__main__":
    main()
