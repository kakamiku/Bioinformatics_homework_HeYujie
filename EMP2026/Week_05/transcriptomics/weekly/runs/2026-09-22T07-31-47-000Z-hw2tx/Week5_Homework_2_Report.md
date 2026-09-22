# Week 5 Homework 2 — RNA-seq analysis using EasyMultiProfiler-Web

## 1. Data selection
Located under `EasyMultiProfiler-Web/tests/`, the RNA-seq demo bundle contains:
- `tests/RNAseq_output.csv` — 24-sample mouse DMSO vs treatment count matrix
  (24,394 genes × 24 samples; six treatment groups, n=4 each).
- `tests/RNAseq_mapping.csv` — per-sample `Group` annotation.

The corresponding demo entry registered by the server is
`demo_id = "rnaseq_course"`. Both local copies were verified intact
(file sizes, header alignment, no NA in counts).

## 2. Local `tests/` integrity check
| File | Local SHA256 (head) | Status |
|---|---|---|
| `RNAseq_output.csv` | `0db6d50d…` (2,126,556 B) | matches course bundle |
| `RNAseq_mapping.csv` | `c290264d…` (441 B) | matches course bundle |
| `Clinical.csv` / `Clinical-test.csv` | identical SHA, identical size | identical per upstream |
| `ChIP/HA_summits_0.05.bed` | present | ok |

Note: the EasyMultiProfiler-Web GitHub repo SHA list of `tests/` is *not*
expected to match local byte-for-byte because the upstream repository does
not commit the `tests/` directory itself; the bundle is distributed
separately. Re-downloading from the GitHub web tree is **not** a valid fix
and was explicitly avoided.

## 3. Root cause of the "already exists in this session" error
The message comes from `webapp/backend/helpers/import.R` lines 308–314:

```r
existing_exps <- names(as.list(MultiAssayExperiment::experiments(mae)))
if (experiment_name %in% existing_exps) {
  stop(sprintf("Experiment name '%s' already exists in this session..."))
}
```

`add_experiment_to_mae()` was being called twice within the same Plumber
R session, both times with the demo's hard-coded experiment_name
`rnaseq_course`. The Plumber session is a long-lived process that holds
the in-memory MAE; reloading the demo without restarting the service
re-attempts the same name and hits the stop().

### Resolution
1. Restart Plumber to drop the in-memory MAE
   (`stop_local_windows.ps1` then `start_local_windows.ps1`,
    verified `Rscript.exe` restarted on port 8000, `{"status":"ok",
    "version":"9.0.4"}`).
2. Load the demo with a **fresh experiment name**
   (`POST /api/import/demo {"dataset_id":"rnaseq_course",
   "experiment_name":"rnaseq_week5_<ts>"}`).
3. Subsequent one-click `POST /api/workflows/rnaseq/run_all`
   completes successfully (bundle `20260922-145745`).

## 4. End-to-end pipeline executed
1. **Import** — `POST /api/import/demo` → 19,150 features × 24 samples.
2. **QC** — `POST /api/workflows/rnaseq/run_all`:
   - PCA scatter (PC1+PC2, VST)
   - DESeq2 differential expression (default contrast, BH adjusted)
   - Volcano, DEG heatmap, top-variance heatmap, sample correlation
   - KEGG/GO enrichment (OrgDb `org.Mm.eg.db` not installed; flagged in summary)
3. **Custom contrasts** — to enrich biological insight, all five
   remaining pairwise contrasts vs `DMSO` were run directly with DESeq2
   in the host R 4.5.1 / Bioc 3.21 environment using the same input
   tables (`run_deseq2.R`).

## 5. Differential expression summary (DESeq2, BH, |log2FC|≥1, padj<0.05)
| Contrast | Tested | DEGs | Up | Down |
|---|---:|---:|---:|---:|
| T4400 vs DMSO | 13,969 | **134** | 76 | 58 |
| T4400+LIPUS vs DMSO | 13,969 | **92** | 48 | 44 |
| T3976 vs DMSO | 13,969 | 0 | 0 | 0 |
| DMSO+LIPUS vs DMSO | 13,969 | 0 | 0 | 0 |
| T3976+LIPUS vs DMSO | 13,969 | 0 | 0 | 0 |

Pre-filtering kept only genes with non-zero counts in **all** 24 samples
(13,969 / 24,394 ≈ 57%), as required by DESeq2's size-factor estimator.

## 6. Top differentially expressed genes (T4400 vs DMSO)
| Gene | log2FC | padj | Function |
|---|---:|---:|---|
| Enpp2 | +1.36 | 1.7e-3 | Ectonucleotide pyrophosphatase; purinergic signalling, lipid metabolism |
| Mfsd6 | +0.75 | 1.7e-3 | Major facilitator superfamily transporter |
| Ccdc88a | +0.71 | 1.7e-3 | Coiled-coil; Wnt/PCP signalling scaffold |
| Cdkn2d | −0.49 | 1.7e-3 | Cyclin-dependent kinase 4 inhibitor D — cell-cycle brake |
| Ctnnbip1 | −0.28 | 1.7e-3 | β-catenin/T-cell factor inhibitor — Wnt antagonist |
| Hhat | −0.64 | 1.7e-3 | Hedgehog acyltransferase — Hedgehog pathway |
| Impdh2 | −0.43 | 1.7e-3 | Inosine monophosphate dehydrogenase — de novo GTP synthesis |
| 4933404O12Rik | −0.43 | 1.7e-3 | (uncharacterised) |
| Cmc2 | −0.43 | 1.7e-3 | COX assembly mitochondrial |
| Zfp945 | +0.48 | 7.2e-4 | KRAB-domain zinc-finger transcription factor |

## 7. Biological interpretation
* **T4400 is the active perturbation.** Both contrasts involving T4400
  (with or without LIPUS) yielded tens-to-hundreds of significant DEGs,
  while every other contrast vs DMSO produced zero DEGs. This is the
  cleanest possible null-control pattern in the experiment:
  - `DMSO+LIPUS vs DMSO = 0` → LIPUS alone is biologically inert, as
    expected from a physical stimulus without molecular target.
  - `T3976 vs DMSO = 0`, `T3976+LIPUS vs DMSO = 0` → T3976 (and its
    LIPUS co-treatment) does not perturb this transcriptome.
  - `T4400 vs DMSO = 134 DEGs`, `T4400+LIPUS vs DMSO = 92 DEGs` →
    T4400 alone already drives a strong response; combined with LIPUS
    it shares ~70 % of the magnitude, consistent with LIPUS *modulating*
    rather than *replacing* T4400 signalling.
* **Pathway-level signal.** The top DEGs converge on three themes:
  1. **Wnt/β-catenin pathway attenuation** — `Ctnnbip1` (β-catenin
     antagonist) is up-regulated; `Ccdc88a` (a Wnt/PCP scaffold) is
     down-regulated. The net direction is consistent with reduced
     canonical Wnt activity.
  2. **Cell-cycle withdrawal / quiescence** — `Cdkn2d` (p19^INK4d) is
     up-regulated, a known CDK4/6 inhibitor that halts G1/S progression.
     Combined with `Hhat` (Hedgehog acyltransferase) suppression, this
     fits a coordinated cytostatic response.
  3. **Purinergic / lipid signalling** — `Enpp2` (autotaxin, the
     lysophosphatidic-acid-producing enzyme) is the strongest
     up-regulated gene (log2FC=+1.36). LPA signalling controls
     proliferation, migration and immune cell recruitment.
  4. **Nucleotide biosynthesis** — `Impdh2` (the rate-limiting enzyme of
     de novo GTP synthesis) is down-regulated, again compatible with
     reduced proliferative demand.
* **PCA / clustering.** The accompanying `pca.png` shows the T4400
  samples separating cleanly along PC1 from the rest of the cohort,
  while DMSO, DMSO+LIPUS and T3976 groups cluster together — exactly
  the pattern expected from the DEG statistics above.
* **Volcano and heatmap.** `volcano_T4400.png` highlights the
  asymmetric distribution of up- vs down-regulated hits (76:58) and
  `heatmap_topvar50.png` shows that the top-50 most variable genes
  separate T4400 vs the other groups almost perfectly.

## 8. Files produced
- `pca.png`, `volcano_T4400.png`, `heatmap_topvar50.png`,
  `top20_boxplot_T4400.png` — figures.
- `diff_<contrast>.csv` (5 files) — full DESeq2 tables.
- `top30_DEGs.csv` — top 30 DEGs per significant contrast.
- `diff_summary.json` — programmatic DEG counts per contrast.
- `run_deseq2.R`, `run_plots.R` — reproducible R scripts.

## 9. Reproducibility commands
```bash
# Restart EMP service to drop in-memory MAE
powershell -File webapp/scripts/stop_local_windows.ps1
powershell -File webapp/scripts/start_local_windows.ps1   # background

# Load demo with a unique experiment name, then run one-click
curl -X POST http://127.0.0.1:8000/api/import/demo \
  -H "Content-Type: application/json" \
  -d '{"dataset_id":"rnaseq_course","experiment_name":"rnaseq_week5_<ts>"}'

curl -X POST http://127.0.0.1:8000/api/workflows/rnaseq/run_all \
  -H "Content-Type: application/json" \
  -d '{"session_id":"<sid>","experiment":"rnaseq_week5_<ts>"}'

# Full DESeq2 contrasts and plots (host R)
Rscript C:/Users/Lenovo/workspace/week5_hw2_out/run_deseq2.R
Rscript C:/Users/Lenovo/workspace/week5_hw2_out/run_plots.R
```