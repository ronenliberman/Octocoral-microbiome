# Octocoral microbiomes in Florida

Analysis code for:

> Liberman R, Gómez-Corrales M, Lopez JV. Host species identity drives
> temporally stable octocoral microbiomes in Florida. *FEMS Microbiology
> Ecology*, fiag124.

The bacterial communities (16S rRNA, V4) of four octocorals were sampled at four
reef sites off southeast Florida (Broward and Miami, shallow and deep) in 2023,
2024 and 2025: the branching *Eunicea flexuosa* and *Muricea muricata* and the
encrusting *Briareum asbestinum* and *Erythropodium caribaeorum*, with seawater
controls. Communities differ mainly by host species and are stable across years
and sites. *Endozoicomonas* and *Mycoplasma* dominate the branching species,
and the encrusting species have richer, more mixed communities.

Raw reads are in NCBI BioProject
[PRJNA1373357](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1373357).

## Getting started

You need R 4.4 or later. Clone **outside any OneDrive or Dropbox folder**:

```bash
git clone https://github.com/ronenliberman/Octocoral-microbiome.git
cd Octocoral-microbiome
git config core.hooksPath .githooks       # once per clone
```

Open `Octocoral-microbiome.Rproj`, then run the scripts in the order given in
[HELPER.md](HELPER.md). Every path is relative to the repository, so nothing
needs editing. Outputs are written to `output/`, which is not tracked.

## What is here

| Folder | Contents |
|---|---|
| `raw files/` | QIIME2 outputs the analysis starts from (PEE2 feature table, taxonomy, representative sequences) and the sample metadata |
| `blast/` | Saved NCBI BLAST results used to verify the *Endozoicomonas* and *Mycoplasma* ASVs |
| `phylogeny/bacteria/` | Alignments and IQ-TREE results of the two bacterial trees (Figs 5–6) |
| `phylogeny/` | Host mtMutS phylogenies (Fig S2) |
| `maps/` | Site map (Fig 1) |
| `script/` | Community analysis: diversity, PERMANOVA, core microbiome, ANCOM-BC2 (Figs 2–4, Table 2) |
| `script/Editor_round3/` | ASV verification, stable labels, trees, Figs 3–6 and the supplementary workbook |
| `documentation/` | Project context and analysis notes |

## Where to look

| File | What it is |
|---|---|
| `HELPER.md` | Layout, run order, and what each script reads and writes. |
| `RULES.md` | Every data-handling and analytical decision, with the reason. |
| `AGENTS.md` | How we work here: branches, pull requests, raw data, attribution. Read before committing. |
| `DATA_PROVENANCE.md` | Where each dataset came from. |
| `documentation/PROJECT_CONTEXT.md` | The science, the review history and the open items. |

## Two things to know before you commit

**Raw data and saved BLAST results are read-only.** Nothing under `raw files/`
or `blast/` changes. If a file looks wrong, open an issue.

**Nothing in the history credits an AI tool.** `tools/check-attribution.sh`
enforces it locally (via the hook above) and on every pull request.

## People

Ronen Liberman (Nova Southeastern University) leads the project and wrote the
analysis code. Matías Gómez-Corrales made the site map and the host phylogenies.
Jose V. Lopez supervised the work.

Contact: Ronen Liberman, rliberma@nova.edu

## License

The code is released under the [MIT License](LICENSE).
