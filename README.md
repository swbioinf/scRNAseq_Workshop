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

## Rendering the bookdown

```r
install.packages(c("bookdown", "downlit", "pander")) # if required
Rscript -e 'rmarkdown::render("scripts/kang2018_preprocessing.Rmd")'
Rscript -e 'bookdown::render_book'
```

Output is written to `docs/`.

## Environment setup 

1. Install CRAN, bioconductor, and github packages:

```bash
bash setup/setup.sh install
```

Review any failed package installs and resolve manually in an
interactive R session. Usually due to missing dependencies in the output
log.

2. Move relevant files and folder structure for the VM.

For small files packaged with the repo:

```bash
bash setup/setup.sh prep_data
```

Transfer the rest from a local machine, stored on 2026-07 scRNAseq/data Google Drive:

```bash
rsync -vPrla *.qs2 tdevXX@XXX.XXX.XXX.XX:/home/tdevXX/data
```

3. Start up the rstudio server:

```
module load rstudio
sudo rstudio-server start
```

* Set a password if required

4. Open rstudio in the browser with the url provided from `rstudio-server start`.

5. Test everything is ready by running the workshop end-to-end

6. Move files and libraries to `/etc/skel/`

7. Clear the rstudio server session history before snapshotting and cloning:

```
sudo rm -rfv /etc/skel/.local/share/rstudio/session
```
