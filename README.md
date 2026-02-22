# MitoReach 2.0 — Web Interface
**Ekker Lab, UT Austin**

Streamlit web interface for MitoReach 2.0 mitochondrial base editing predictor.

## Usage

1. Select species (Zebrafish / Human / Mouse)
2. Select a gene OR enter genomic coordinates
3. Click **Run Physics Engine** on any target
4. Results appear inline (heatmap + favorability profiles)

## Deploy on Streamlit Community Cloud

1. Push this folder to a GitHub repo
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect your GitHub account → select this repo → `app.py` as main file
4. Deploy — get a permanent public URL

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Files

```
app.py                      ← Streamlit UI (wrapper only, engine untouched)
requirements.txt            ← Python dependencies
data/
  9jo8.chimera.pdb          ← TALE-DddA crystal structure + DNA (physics template)
  zebrafish_mtDNA.gb        ← D. rerio mitochondrial genome
  human_mtDNA.gb            ← H. sapiens mitochondrial genome
  mouse_mtDNA.gb            ← M. musculus mitochondrial genome
mitoreach_source/
  MitoReach2_Standalone.py  ← MitoReach 2.0 physics engine (unchanged)
  genome_parser.py          ← GenBank genome parser (unchanged)
  tc_motif_scanner.py       ← TC motif utilities (unchanged)
```

## Citation

If you use this tool, please cite the original MitoReach 2.0 paper (Antigravity/tool authors)
and acknowledge the Ekker Lab TC-constraint fix (Feb 2026).
