"""
MitoReach 2.0 — Web Interface
Ekker Lab, UT Austin

Minimal wrapper around MitoReach2_Standalone.py + genome_parser.py.
No changes to the engine. Just pick species, enter a position, and run.
"""

import streamlit as st
import subprocess
import sys
import os
import glob
from pathlib import Path
from PIL import Image

# ── Paths ──────────────────────────────────────────────────────────────────
APP_DIR  = os.path.dirname(os.path.abspath(__file__))
SRC_DIR  = os.path.join(APP_DIR, 'mitoreach_source')
DATA_DIR = os.path.join(APP_DIR, 'data')
ENGINE   = os.path.join(SRC_DIR, 'MitoReach2_Standalone.py')
PDB_FILE = os.path.join(DATA_DIR, '9jo8.chimera.pdb')
RESULTS  = os.path.join(APP_DIR, 'results')
os.makedirs(RESULTS, exist_ok=True)

sys.path.insert(0, SRC_DIR)

GENOMES = {
    "Zebrafish (D. rerio)":  os.path.join(DATA_DIR, 'zebrafish_mtDNA.gb'),
    "Human (H. sapiens)":    os.path.join(DATA_DIR, 'human_mtDNA.gb'),
    "Mouse (M. musculus)":   os.path.join(DATA_DIR, 'mouse_mtDNA.gb'),
}

# ── Page config ─────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="MitoReach 2.0 — Ekker Lab",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 MitoReach 2.0")
st.caption("Ekker Lab · UT Austin · Enter a target position and run")

# ── Genome loader ────────────────────────────────────────────────────────────
@st.cache_resource
def load_genome(path):
    from genome_parser import GenomeParser
    return GenomeParser(path)

# ══════════════════════════════════════════════════════════════════════════════
# SIDEBAR
# ══════════════════════════════════════════════════════════════════════════════
with st.sidebar:
    st.header("⚙️ Settings")

    species = st.selectbox("Species:", list(GENOMES.keys()))
    genome_path = GENOMES[species]

    variant = st.radio(
        "DddA variant:",
        ["DddA11 (any C)", "WT (TC-constrained)"],
        index=0
    )
    variant_flag = "WT" if variant.startswith("WT") else "DddA11"

    st.markdown("---")
    st.caption("Engine: MitoReach2_Standalone.py\nTemplate: 9jo8.chimera.pdb")

# ── Load genome ───────────────────────────────────────────────────────────────
try:
    genome = load_genome(genome_path)
except Exception as e:
    st.error(f"Could not load genome: {e}")
    st.stop()

# ══════════════════════════════════════════════════════════════════════════════
# TABS
# ══════════════════════════════════════════════════════════════════════════════
tab_run, tab_genes, tab_results = st.tabs([
    "▶ Run Analysis",
    "📋 Gene List",
    "📁 Past Results"
])

# ──────────────────────────────────────────────────────────────────────────────
# TAB 1 — RUN
# ──────────────────────────────────────────────────────────────────────────────
with tab_run:
    st.subheader("Target Selection")

    input_mode = st.radio(
        "Input by:",
        ["Gene name", "Genomic coordinate"],
        horizontal=True
    )

    target_position = None  # genomic position of target C (1-indexed)

    if input_mode == "Gene name":
        gene_list = sorted([g['name'] for g in genome.list_genes()])
        col1, col2 = st.columns([2, 1])
        with col1:
            gene_sel = st.selectbox("Select gene:", [""] + gene_list)
        with col2:
            offset = st.number_input(
                "Offset within gene (bp from start):",
                min_value=0, value=0,
                help="0 = use gene start position"
            )

        if gene_sel:
            ginfo = genome.get_gene_info(gene_sel)
            if ginfo:
                target_position = ginfo['start'] + offset
                st.info(
                    f"**{gene_sel}** → {ginfo['start']}–{ginfo['end']} "
                    f"({ginfo['end'] - ginfo['start'] + 1} bp)  \n"
                    f"Target position: **{target_position}**"
                )

    else:
        col1, col2 = st.columns([2, 1])
        with col1:
            target_position = st.number_input(
                "Target C position in mtDNA (1-indexed):",
                min_value=1,
                max_value=genome.genome_length,
                value=9007,
                help=f"e.g. 9007 — the mtDNA position of your target cytosine"
            )
        with col2:
            flank = st.number_input(
                "Flanking bp (each side):",
                min_value=10, max_value=200, value=50,
                help="Sequence context extracted around the target"
            )

    # ── Run ──────────────────────────────────────────────────────────────────
    if target_position:

        # Determine flanking window
        if input_mode == "Gene name" and gene_sel:
            flank = 50  # default for gene mode

        region_start = max(1, target_position - flank)
        region_end   = min(genome.genome_length, target_position + flank)

        # Extract sequence
        region_seq = genome.extract_sequence(region_start, region_end)

        # 0-indexed position of target within extracted sequence
        target_idx_0 = target_position - region_start

        st.caption(
            f"Extracted region: {region_start}–{region_end} "
            f"({len(region_seq)} bp) | "
            f"Target index in sequence: {target_idx_0} (0-indexed)"
        )

        if st.button("▶ Run MitoReach", type="primary"):

            out_label = f"pos{target_position}_{variant_flag}"
            out_dir   = os.path.join(RESULTS, out_label)
            os.makedirs(out_dir, exist_ok=True)

            cmd = [
                sys.executable, ENGINE,
                "--sequence",   region_seq,
                "--target-pos", str(target_idx_0),
                "--variant",    variant_flag,
                "--template",   PDB_FILE,
                "--output-dir", out_dir,
                "--min-len",    "10",
                "--max-len",    "18",
            ]

            with st.spinner(f"Running MitoReach on position {target_position}… (~2–5 min)"):
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=600
                )

            if result.returncode != 0:
                st.error("Engine error:")
                st.code(result.stderr[-2000:] if result.stderr else "(no stderr)")
            else:
                st.success("✅ Done!")
                st.code(result.stdout)

                # Show output PNGs generated by the engine
                target_subdir = os.path.join(out_dir, f"Target_{target_idx_0}")
                pngs = sorted(glob.glob(os.path.join(target_subdir, "*.png")))

                if pngs:
                    heatmaps = [p for p in pngs if 'heatmap' in p]
                    profiles = [p for p in pngs if 'profile' in p]

                    if heatmaps:
                        st.subheader("Structural Landscape Heatmap")
                        st.image(Image.open(heatmaps[0]), use_container_width=True)

                    if profiles:
                        st.subheader("Favorability Profiles by Spacer Length")
                        for i in range(0, len(profiles), 2):
                            c1, c2 = st.columns(2)
                            with c1:
                                st.image(
                                    Image.open(profiles[i]),
                                    caption=os.path.basename(profiles[i]),
                                    use_container_width=True
                                )
                            if i + 1 < len(profiles):
                                with c2:
                                    st.image(
                                        Image.open(profiles[i + 1]),
                                        caption=os.path.basename(profiles[i + 1]),
                                        use_container_width=True
                                    )
                else:
                    st.warning("Run completed but no plots were generated. Check the log above.")

# ──────────────────────────────────────────────────────────────────────────────
# TAB 2 — GENE LIST
# ──────────────────────────────────────────────────────────────────────────────
with tab_genes:
    st.subheader(f"{species} — mitochondrial genes ({genome.genome_length} bp)")

    import pandas as pd
    rows = []
    for g in genome.list_genes():
        rows.append({
            'Gene':       g['name'],
            'Start':      g['start'],
            'End':        g['end'],
            'Length (bp)': g['end'] - g['start'] + 1,
            'Strand':     '+' if g.get('strand', 1) == 1 else '-',
            'Product':    g.get('product', ''),
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, height=500)

# ──────────────────────────────────────────────────────────────────────────────
# TAB 3 — PAST RESULTS
# ──────────────────────────────────────────────────────────────────────────────
with tab_results:
    st.subheader("Previously run analyses")

    result_dirs = sorted(glob.glob(os.path.join(RESULTS, "*")), reverse=True)

    if not result_dirs:
        st.info("No results yet. Run an analysis in the ▶ Run Analysis tab.")
    else:
        for rdir in result_dirs:
            label = os.path.basename(rdir)
            pngs  = sorted(glob.glob(os.path.join(rdir, "**", "*.png"), recursive=True))

            with st.expander(f"📂 {label}  ({len(pngs)} plot(s))"):
                if pngs:
                    heatmaps = [p for p in pngs if 'heatmap' in p]
                    profiles = [p for p in pngs if 'profile' in p]

                    if heatmaps:
                        st.image(Image.open(heatmaps[0]), use_container_width=True)

                    for i in range(0, len(profiles), 2):
                        c1, c2 = st.columns(2)
                        with c1:
                            st.image(Image.open(profiles[i]),
                                     caption=os.path.basename(profiles[i]),
                                     use_container_width=True)
                        if i + 1 < len(profiles):
                            with c2:
                                st.image(Image.open(profiles[i + 1]),
                                         caption=os.path.basename(profiles[i + 1]),
                                         use_container_width=True)
                else:
                    st.caption("No plots found.")
