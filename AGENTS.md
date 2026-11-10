# Working rules for this repository

This is how we work in this repository. Ronen Liberman leads the project and
owns the analysis code. Matías Gómez-Corrales contributed the site map
(`maps/`) and the host phylogenies (`phylogeny/*.Rmd`), and Jose V. Lopez
supervises. Anyone can change anything below by opening a pull request against
this file, and that is the intended way to disagree with any of it.

The repository is the code record for a published paper: Liberman R,
Gómez-Corrales M, Lopez JV, "Host species identity drives temporally stable
octocoral microbiomes in Florida", *FEMS Microbiology Ecology* (fiag124). The
paper's Data availability statement points readers here. A change to a script
that produced a published number changes what the paper's code record says, so
treat the analysis scripts with that in mind.

Three of the rules are absolute rather than preferential, and they are marked
that way. Breaking them destroys something that cannot be recovered: a
permanent branch name on a pull request, the only copy of an original data
file, or the match between the published paper and the code that produced it.
Everything else is a convention, and it can be argued with.

If you are an AI coding assistant working for any of us, read this file in full
before your first edit and follow it ahead of your own defaults. Where your
habits and this file disagree, follow this file, and raise the disagreement with
the person you are working for rather than resolving it yourself. Then read
`HELPER.md` (how the project is put together), `RULES.md` (how each dataset is
handled and why) and `documentation/PROJECT_CONTEXT.md` (the science, and the
decisions that are settled).

---

## 1. Never credit an AI tool, anywhere (absolute)

The history of this repository must read as the work of its human authors. We
use assistants, and nobody is coy about that. The point is narrower: the commit
log, the contributor graph, and the pull request record are how a reader works
out who is accountable for the science, and a tool is not accountable for
anything.

There are no exceptions to this and no place it does not apply.

Do not write the names Claude, Anthropic, Codex, Copilot, Cursor, ChatGPT,
Gemini, or any equivalent, and do not write phrases such as "generated with",
"AI-assisted", "written with the help of", or a `Co-Authored-By:` trailer naming
any tool or bot, in any of these:

- commit messages, including trailers at the bottom
- branch names
- pull request titles, bodies, and comments
- issue titles, bodies, and comments
- code comments, documentation, and every other committed file

Three configuration paths may have tool names in them, and naming those paths is
fine: `CLAUDE.md`, the `.claude/` directory, and
`.github/copilot-instructions.md`. The automated check allows those exact
filenames and nothing else.

**Branch names need care before your first push.** GitHub writes the head branch
name into a pull request's permanent record, and nobody can edit it afterward.
Assistants name new branches with their own prefixes by default. Check and
rename before you push:

```bash
git symbolic-ref --short HEAD
git branch -m <bad-name> ronen/what-youre-actually-doing
```

Turn on the local hook once per clone, so that a bad name or message is caught
before the push:

```bash
git config core.hooksPath .githooks
```

The hook and the pull request check both run `tools/check-attribution.sh`, so
the rule has one definition rather than two that can drift.

Human contributors are credited, and that is the other half of the same rule.
When squash merging a pull request someone else wrote, add a
`Co-authored-by: Full Name <email>` trailer, one per human author.

## 2. Raw data and saved NCBI results are read-only (absolute)

- `raw files/` holds the QIIME2 outputs the scripts start from (the PEE2
  feature table, taxonomy and representative sequences) and the sample
  metadata. Only these files belong there: no processed tables, no figures.
- `blast/` holds the NCBI BLAST hit tables, search results and taxonomy cache
  that the validation scripts read. NCBI's databases change, so these files are
  the only way to reproduce the published ASV validation.

**Never edit, reformat, clean, or re-save a file in either directory.** Not to
fix an encoding, not to correct an obvious typo, not to tidy a column name. When
a file looks wrong, open a Data question issue, handle the problem in the script
that reads it, and write down what you did in `RULES.md`.

The one exception is deliberate: re-running a BLAST script (`01b`, `01b_Myco`,
`01c`, `01d`) rewrites its own files in `blast/`. Do that only on purpose,
through an issue, because it replaces the searches the paper is based on.

`DATA_PROVENANCE.md` records where each dataset came from.

## 3. The published results stay reproducible (absolute)

`script/` and `script/Editor_round3/` produced the published figures, tables and
numbers. Do not change what any of them computes (a threshold, a filter, a
seed, a reference list, a model) on `main` without an Analysis issue that says
which published result moves and by how much.

If the paper ever needs a different analysis, add a new script or a new
versioned output rather than silently changing the old one, so the published
record can still be rerun. Fixing a path, a missing `library()` call or a crash
that does not change any output is fine; say so in the commit message.

## 4. Settled decisions stay settled until someone reopens them on the record

`RULES.md` lists the analytical decisions that were deliberated and closed
during three rounds of review: the PEE2 table, the three excluded samples, the
normalization per test, the BLAST validation thresholds (Yarza et al. 2014), the
fixed core ASV labels (Endo_1–6, Myco_1–10), the type-strain references, the
tree settings, and others. They are not irreversible, but they were argued out
once and the reasoning is written down.

Do not quietly change one because a different value looks better. Open an
Analysis issue, say what the current decision is and why it should move, and
change `RULES.md` in the same pull request as the code.

## 5. Branches and pull requests

After the initial upload, nothing lands on `main` directly. Every change
arrives through a pull request.

Name branches `<yourname>/<short-topic>`, lowercase with hyphens. Ronen uses
`ronen/`, Matías uses `matias/`. Never a tool name.

```bash
git checkout main && git pull && git checkout -b ronen/what-youre-doing
git push -u origin ronen/what-youre-doing
```

- **One pull request does one thing.** A data fix and an analysis change are two
  pull requests.
- **Link the issue** it closes, with `Closes #12` in the body.
- **Fill in the template.** It is short and every box matters.
- **Squash merge**, then delete the branch.
- A change to `maps/` or the host phylogeny Rmds waits for Matías's review.

Write commit messages as a summary line under about seventy characters, then a
blank line, then prose explaining why the change was needed. Describe the
change, not the process that produced it.

## 6. Issues come first

Open an issue before starting anything beyond a typo fix.

| Template | Use it for |
|----------|-----------|
| Data question | Something looks wrong in a raw file, a saved BLAST result, or the metadata; a sample whose identity is in doubt |
| Analysis | A statistical choice, a change to any script that writes a published result, or a figure or number in the paper |

Say what you expected, what you got, and which script or file you were in.

## 7. The pipeline runs in order, and outputs are not committed

`HELPER.md` gives the run order and what each script reads and writes. Every
script writes into `output/`, which git ignores: figures, tables and logs are
rebuilt by running the scripts, never committed and never edited by hand.

Two kinds of file are committed although a script can write them, because the
paper depends on the exact versions: the saved NCBI results in `blast/`, and the
tree inputs and IQ-TREE results in `phylogeny/bacteria/`. Re-running `03`/`04`
rewrites `*_for_tree.fasta` and `*_reference_names.csv` there; `git diff` should
show no change. If it does, stop and find out why before committing.

Start R sessions empty. Some round-2 scripts (`Supplementary_Materials.R`, the
community composition script) use objects created by
`Octocoral_Revised_Analysis_COMPLETE.R` and must be run in the same session
after it; `HELPER.md` says which.

## 8. Keep the written record true

Update these in the same commit as the change they describe:

- `RULES.md` when you change how a dataset is handled or reopen a decision.
- `HELPER.md` when you change the structure of the project or the run order.
- `documentation/PROJECT_CONTEXT.md` when a result it quotes changes.
- The README in `raw files/`, `blast/` or `phylogeny/bacteria/` when a file
  there is added or replaced.

## 9. What to do when you are unsure

Ask. Open an issue, or leave the question in the pull request. Do not invent a
workaround for something you do not understand, do not weaken a check to make
something pass, and do not widen a pull request past what its issue asked for.
