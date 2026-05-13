# Single Cell RNAseq Workshop

Workshop materials for *scRNAseq Analysis in R with Seurat*, run by members of QCIF, SIH, Monash, and the Australian BioCommons. Topics covered include QC, normalisation, dimensionality reduction, batch correction (Harmony), clustering, cell type annotation (SingleR), and differential expression.

Visit the [Workshop page](https://swbioinf.github.io/scRNAseq_Workshop/).

## Repository structure

| Path | Description |
|---|---|
| `part1/` and `part2/` | Workshop sections (*.Rmd files) that learners work through sequentially with follow-along and self-guided exercises  |
| `data/` | Input datasets used during the workshop |
| `scripts/` | Helper and preprocessing scripts |
| `setup/` | Training environment (VM) setup files |
| `_bookdown.yml`, `_output.yml`, `_common.R`, etc. | Auxiliary files for rendering the bookdown site |

## Environment setup 

1. Install CRAN, bioconductor, and github packages:

```bash
bash setup/setup.sh install_packages
```

Review any failed package installs and resolve manually in an
interactive R session. Usually due to missing dependencies in the output
log.

2. Move relevant files and folder structure for the VM:

```bash
bash setup/setup.sh prep_data
```

3. Remove unneeded files and folders:

```bash
bash setup/setup.sh clean
```

4. Start up the rstudio server:

```
module load rstudio
sudo rstudio-server start
```

5. Open rstudio in the browser with the url provided from `rstudio-server start`.

