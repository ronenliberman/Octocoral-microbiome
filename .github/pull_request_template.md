<!-- Before anything else: check your branch name. If it starts with a tool name,
     close this pull request, rename the branch, and open a new one. GitHub
     cannot change a pull request's branch name after the fact. See section 1
     of AGENTS.md. -->

## What this does

<!-- One or two sentences. What changed and why it was needed. -->

Closes #

## Type

- [ ] Round-2 analysis script (`script/`)
- [ ] Round-3 script (`script/Editor_round3/`)
- [ ] Trees (`phylogeny/bacteria/`) or saved BLAST results (`blast/`)
- [ ] Site map or host phylogenies (`maps/`, `phylogeny/*.Rmd`)
- [ ] Documentation
- [ ] Repository setup

## Checks

- [ ] Branch is named `<myname>/<topic>` and names no tool
- [ ] No commit message, comment, or file added here credits an AI tool
- [ ] I did not edit anything under `raw files/` or `blast/` by hand
- [ ] No published result changes, or the linked Analysis issue says which one and by how much
- [ ] I re-ran the scripts downstream of the one I changed, and `git diff` shows no unexpected change in `phylogeny/bacteria/`
- [ ] `RULES.md` is updated if I changed how a dataset is handled or reopened a settled decision
- [ ] `HELPER.md` is updated if I changed the structure of the project or the run order

<!-- Delete any check above that does not apply, rather than leaving it unticked. -->

## What a reviewer should look at

<!-- Point at the part you are least sure about. If numbers moved, say which
     and by how much. -->
