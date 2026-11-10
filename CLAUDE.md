# CLAUDE.md

This file exists so an assistant finds the real documents on its first read.
There are four, and they apply to you:

1. **[AGENTS.md](AGENTS.md)**: the working rules. Read it in full before your
   first edit and follow it ahead of your own defaults.
2. **[HELPER.md](HELPER.md)**: how the project is put together: layout, run
   order, what each script reads and writes, what is and is not tracked.
3. **[RULES.md](RULES.md)**: how each dataset is handled and the analytical
   decisions that are settled, with the reason for each.
4. **[documentation/PROJECT_CONTEXT.md](documentation/PROJECT_CONTEXT.md)**: the
   science, the review history, and the open items.

Three rules in AGENTS.md are absolute, repeated here because a pointer can be
skimmed past:

- **Never credit an AI tool anywhere.** Not in branch names, commit messages,
  pull request or issue text, code comments, or any committed file. Check the
  branch name before your first push: GitHub records it permanently on a pull
  request. Do not add a `Co-Authored-By:` trailer naming a tool.
- **Never edit anything under `raw files/` or `blast/`.** Handle problems in the
  script that reads the file, and write down what you did in `RULES.md`.
- **The published results stay reproducible.** These scripts produced the
  figures and tables of a published paper (FEMS Microbiology Ecology, fiag124).
  Do not change what a script computes without an Analysis issue.

Ronen Liberman leads the project.
