# MERFISH Foxf2 analysis

Analysis of MERFISH data from Foxf2 endothelial knockout (ECKO) and control (WT) mice: NeuN-low glutamatergic neurons, endothelial cells (ECs) near them, and cell-level differential expression.

## Setup
- Data (LRZ): `/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/analysis/foxf2/Cellpose_1_nuclei_model/`
- Input `adatas/adata_integrated.h5ad.gz` (segmentation, QC, integration, cell types) is from [cellseg-benchmark](https://github.com/simonmfr/cellseg-benchmark), which the notebooks also import.
- Environment: container `misc/enroot_images/benchmark_new2.sqsh` (env `analysis`) plus
  `mamba install -n analysis -c conda-forge -c bioconda mygene xlsxwriter openpyxl papermill shapely "numpy<2.5"`
- `05` needs an env with `spatialdata spatialdata-plot scanpy matplotlib-scalebar`.

## Notebooks
| | |
|---|---|
| `00_spatial_overview` | spatial plots, cell types |
| `01_foxf2_expression` | Foxf2 knockout check |
| `02_neuron_annotation` | neuron subclustering and annotation (WT, ECKO) |
| `03_neuron_annotation_pcko_glko` | same for PCKO/GLKO, labels for the joint DEA |
| `04_neun_low` | NeuN-low definition, % and density per sample, EC distance to NeuN-low neurons, dotplots |
| `05_boundary_plots` | cell boundaries with the EC distance cutoff |
| `06_neun_low_cutoffs` | NeuN-low labels at 5-25% cutoffs (input for 07) |
| `07_dea` | cell-level DEA, one comparison per run |
| `08_dea_sensitivity` | DEA across NeuN-low and distance cutoffs |

Notebooks read the dated `.h5ad.gz` files written by earlier notebooks; update the dates when rerunning.

## DEA
MAST with a random effect per sample, ≤5,000 cells per condition.
- ECKO vs WT: `~ genotype + slide + (1 | sample)`, fitted on all genotypes, testing ECKO vs WT
- within ECKO (NeuN-low vs NeuN+, proximal vs distal ECs): `~ condition + (1 | sample)`

Tables: FC and FDR from MAST. `08` adds log2FC from group means (`log2FC_scanpy`).

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
pm="papermill -k python3 --log-output 07_dea.ipynb"
mast="-p method MAST-RE -p max_cells 5000"
joint="-p batch slide -p joint True"
mkdir -p runs

$pm runs/ko.ipynb -p comp KO_vs_WT $mast $joint -p pooled False
for d in 18 31 58 103; do
  $pm runs/F$d.ipynb -p comp ECs_proximal_vs_distal_ECKO -p dist $d $mast
done
for p in 5 10 15 20 25; do
  $pm runs/E2_p$p.ipynb -p comp lowNeuN_ECKO_vs_Glut_WT -p low_pct $p $mast $joint -y "subsets: [Neurons_Glut]"
  $pm runs/E3_p$p.ipynb -p comp lowNeuN_vs_normalNeuN_ECKO -p low_pct $p $mast
  $pm runs/F18_p$p.ipynb -p comp ECs_proximal_vs_distal_ECKO -p dist 18 -p low_pct $p $mast
  $pm runs/G_p$p.ipynb -p comp ECs_proximal_ECKO_vs_WT -p dist 18 -p low_pct $p $mast $joint
done
```
`-p method wilcoxon_scanpy` gives wilcoxon tables. `-p overwrite False` reuses cached MAST fits but skips existing xlsx files.

## Outputs
In the data folder:
- `dea/final/`: tables used in the paper
- `dea/MAST-RE_5k/`, `dea/wilcoxon_scanpy/`: all DEA tables and volcano plots; `sensitivity/`: heatmaps and gene tables (08)
- `neuronal-subclustering/`: dotplots, UMAPs, NeuN-low % and density, spatial plots
