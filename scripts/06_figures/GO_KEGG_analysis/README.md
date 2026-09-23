# GO and KEGG enrichment analysis — OsHV-1 infection in *Crassostrea gigas* single-nucleus RNA-seq

Functional enrichment pipeline accompanying the analysis of single-nucleus transcriptomes from Pacific oysters (*Crassostrea gigas*) sampled through a time course of OsHV-1 infection.

The pipeline covers two independent questions:

- **What changes during infection?** Control versus infected nuclei within each transcriptomic cluster, at 6, 24, 72 and 96 hours post-infection.
- **What defines Cluster 1?** An unassigned cluster characterised using marker genes computed on uninfected animals only.

---
```
00a → 00b → 01 ─┬→ 02 ─┬→ 03b → 07        (Figure 7)
                │      └→ 03a              (optional, supplementary)
                ├→ 03d → 03c               (diagnostics)
                └→ 05a → 05b → 06          (Cluster 1)
```
---

## Requirements

**HPC (steps 00a–00b)** — `gffread`, `seqkit`, eggNOG-mapper v2 with the eggNOG database, DIAMOND.

**R, Seurat environment** — R ≥ 4.2, Seurat v5, `Matrix`, `dplyr`, `readr`, `tibble`, `ggplot2`.

**R, enrichment environment** — `clusterProfiler`, `GO.db`, `AnnotationDbi`, `fgsea`, `dplyr`, `readr`, `tidyr`, `ggplot2`, `stringr`.

Seurat and clusterProfiler are kept in separate conda environments; the two never need to be loaded together. Scripts exchange data as files on disk.

---

## Pipeline

| Script | Environment | Purpose |
|---|---|---|
| `00a_prepare_proteome.sh` | HPC | Generate a protein FASTA from the genome and GFF3, select the longest isoform per gene, remove mitochondrial sequences containing invalid characters |
| `00b_run_eggnog_mapper.sh` | HPC | Run eggNOG-mapper on the full proteome (Metazoa scope, one-to-one orthologues), producing `full_proteome.emapper.annotations` |
| `01_build_reference_tables.R` | enrichment | Build GO and KEGG gene set tables from the eggNOG output |
| `02_rerun_de_full_ranked.R` | Seurat | Differential expression between control and each infection stage, within each cluster, with no fold-change threshold |
| `03a_ora_across_stages.R` | enrichment | Over-representation analysis of thresholded differentially expressed gene lists |
| `03b_gsea_across_stages.R` | enrichment | Gene set enrichment analysis on full ranked gene lists |
| `03c_gsea_across_stages_pseudoRepl.R` | enrichment | The same GSEA applied to the between-animal diagnostic comparisons |
| `03d_pseudoreplication_checks_1.R` | Seurat | Control-versus-control, animal-versus-animal and animal-level pseudobulk comparisons |
| `04_dotplot_matrix.R` | enrichment | Dot plot matrices of enriched terms, per-cluster supplementary panels, supplementary tables |
| `05a_cluster1_identity_markers.R` | Seurat | Cluster 1 marker genes from control nuclei only, with matched backgrounds and a per-animal reproducibility check |
| `05b_cluster1_identity_enrichment.R` | enrichment | GO and KEGG enrichment of cluster marker genes, and identification of terms unique to Cluster 1 |
| `06_cluster1_figure_table.R` | Seurat | Cluster 1 supplementary figure and table |

### Execution order

```
00a → 00b → 01 → 02 → 03a, 03b → 04
                  ↓
                 03d → 03c
01 → 05a → 05b → 06
```

Steps `02–04` and `05a–06` are independent of one another; both depend on `01`.

---

## Methodological notes

### Gene set construction

There is no GO or KEGG annotation resource for *C. gigas*, so both were built from the eggNOG-mapper output and supplied to `clusterProfiler` as custom `TERM2GENE` tables. Three points affect the results:

- **GO ancestor propagation.** `enricher()` performs no propagation with a custom `TERM2GENE`, unlike `enrichGO()` with an OrgDb. Terms are therefore propagated explicitly using `GO.db`. In this dataset eggNOG had already supplied a largely propagated annotation, so the step increased the table size by only ~10%.
- **Ontology separation.** GO terms are split into biological process, molecular function and cellular component, so multiple-testing correction is applied within an ontology rather than across all three.
- **KEGG deduplication.** eggNOG reports both `ko#####` and `map#####` identifiers for the same pathway. These are folded into a single KO identifier, and the global overview maps (for example `ko01100 Metabolic pathways`) are excluded. 373 pathways are retained.

### Enrichment background

Over-representation analysis is tested against a **cluster-specific background** comprising the genes eligible for differential expression testing within that cluster, rather than the full annotated transcriptome. Using the whole transcriptome instead returns 1,188 significant terms, of which 628 (53%) are not supported under the cluster-specific background — the difference reflects cell-type-specific expression profiles rather than biology. `03a` reports both so the comparison is auditable.

Gene set enrichment analysis requires no background: the ranked list is itself the universe.

### Choice of method

Differentially expressed genes are annotated at a lower rate (30–34%) than the expressed background within the same clusters (46–50%), leaving thresholded gene lists too small for reliable over-representation testing — a median of 25 annotated genes per comparison against 5,234 testable GO biological process terms. GSEA on full ranked lists is therefore used as the primary pathway-level method, with ORA reported alongside it.

### Differential expression is run twice

`02` is run into two output directories:

| `min.pct` | Directory | Use |
|---|---|---|
| 0.25 | `de_full_ranked` | Gene-level results, matching the differentially expressed gene counts reported in the manuscript |
| 0.1 | `de_full_ranked_minpct01` | Pathway-level analysis; the relaxed threshold gives ranked lists 3–4× longer, which GSEA requires |

`03a`, `03b` and `04` use the `min.pct = 0.1` output.

### Pseudoreplication

One animal was sampled per infection stage, so nuclei within a comparison are not biological replicates. Nucleus-level p-values are used for ranking and candidate selection, not for inference. `03d` and `03c` implement three checks:

- the two uninfected control animals compared against one another;
- the two animals sampled at 6 hpi compared against one another, and each against the pooled controls;
- animal-level pseudobulk profiles, two controls against two 6 hpi animals — the only comparison in the dataset with biological replication on both sides.

### Cluster 1 identity

Marker genes for Cluster 1 are recomputed on control nuclei only. Pooling control and infected nuclei introduces infection-responsive genes into the identity signature, because 71% of nuclei in the dataset are from infected animals. Reproducibility is assessed by repeating marker detection within each control animal separately.

Enrichment terms are cross-referenced against the same analysis run for all 18 clusters, to identify terms enriched in Cluster 1 and in no other cluster.

### Interpreting orthology-transferred annotations

GO terms assigned by orthology can carry functional labels that do not apply in the target species. In this dataset, inhibitor-of-apoptosis paralogues are returned under antimicrobial peptide terms, because their insect orthologues act in the IMD pathway upstream of antimicrobial peptide production. Similarly, KEGG human-disease maps such as `ko05145 Toxoplasmosis` and `ko01524 Platinum drug resistance` are populated by conserved signalling and apoptosis genes. **Leading-edge genes should be inspected before any term is reported.**

---

## Configuration

Paths are set at the top of each script and need editing for a new environment. The values that must match the upstream differential expression settings are:

| Parameter | Value | Where |
|---|---|---|
| `min.pct` (pathway analysis) | 0.1 | `02` |
| `min.pct` (gene-level) | 0.25 | `02` |
| `min.pct` (cluster markers) | 0.20 | `05a` |
| Gene set size | 5–500 (GSEA), 10–500 (ORA) | `03a`, `03b` |
| Multiple testing | Benjamini–Hochberg, adjusted p < 0.05 | all |
| GSEA seed | 42 | `03b` |

`03b` is run twice: once on the main comparisons, and once via `03c` on the diagnostic comparisons from `03d`, with `de_dir` and `outdir` changed accordingly.

---

## Outputs

**Main text**

- Dot plot matrices of GO and KEGG terms depleted in infected nuclei (GSEA) and enriched among upregulated genes (ORA), across clusters and infection stages

**Supplementary**

- Per-cluster enrichment panels, all significant terms
- Full GSEA results with leading-edge genes
- Cluster-specific versus global background comparison
- Transcript complexity by cluster and condition
- Cluster 1 IAP paralogue expression across clusters
- Cluster 1 enriched terms with driving genes

---

## Citation

If you use this pipeline, please cite the accompanying manuscript.
