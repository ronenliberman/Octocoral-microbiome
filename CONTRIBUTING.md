# Contributing

The working rules live in [AGENTS.md](AGENTS.md), and that is the file to read
rather than this one. GitHub links this page from the issue and pull request
pages, so it is a signpost rather than a second copy.

## The short version

Open an issue first. Branch as `<yourname>/<topic>`, never with a tool name in
it. One pull request does one thing, links its issue, and gets squash merged.

```bash
git checkout main && git pull && git checkout -b yourname/what-youre-doing
```

## Three things that cannot be undone

1. **Nothing in the history credits an AI tool.** Branch names are the sharp
   edge, because GitHub records a pull request's branch name permanently. A
   pull request check enforces this.
2. **Nothing under `raw files/` or `blast/` gets edited.** Those are the
   original QIIME2 outputs and the saved NCBI searches. Problems get handled in
   the script that reads them.
3. **The published results stay reproducible.** These scripts produced the
   figures and tables of a published paper. A change to what a script computes
   needs an Analysis issue first.

## Before you open the pull request

Re-run the scripts downstream of the one you changed (`HELPER.md` has the
order) and check that `git diff` shows no unexpected change in
`phylogeny/bacteria/`. Update `RULES.md` if you changed how a dataset is
handled, and `HELPER.md` if you changed how the project is put together, in the
same commit as the change itself.
