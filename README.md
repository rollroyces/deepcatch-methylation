# deepcatch-methylation — ARCHIVED on 2026-10-16

> ⚠️ **This repository is archived.** The cfRNA persistence coupling and
> methylation-array reader have been folded into
> [**`rollroyces/deepcatch`**](https://github.com/rollroyces/deepcatch)
> as part of the v2.3.0 release (CANCER-SEEK integration).

## Scope conflict resolution (final)

This repo's `README.md` originally described a **cfRNA read-side companion
to `mrnavax-rare`** while the GitHub repo **description** and the
`METHYLATION_PROJECT.md` content described **CpG methylation-based
cfDNA cancer detection**. Both scopes are real, but the consolidation
opted for the **CpG / FinaleMe** track documented in
[`METHYLATION_PROJECT.md`](https://github.com/rollroyces/deepcatch/blob/main/METHYLATION_PROJECT.md).

The cfRNA coupling becomes `deepcatch/src/methylation/cfrna_coupling.py`
(a structural companion, no runtime dependency on `mrnavax-rare`); the
methylation-array reader becomes `deepcatch/src/methylation/methylation_reader.py`.

| What lived here                                                 | Where it lives now                                |
|-----------------------------------------------------------------|---------------------------------------------------|
| `deepcatch_methylation/cfrna_coupling.py`                         | `deepcatch/src/methylation/cfrna_coupling.py`     |
| `deepcatch_methylation/methylation_reader.py`                     | `deepcatch/src/methylation/methylation_reader.py` |
| `METHYLATION_PROJECT.md`                                         | `deepcatch/METHYLATION_PROJECT.md` (verbatim; SHA `3282f90b`) |
| `PHASE0_PLAN.md`                                                  | `deepcatch/docs/ultrasensitive/PHASE0_PLAN.md`    |
| `BIORXIV_PAPER_FRAGMENTOMICS.{md,pdf}`                            | `deepcatch/paper/methylation/BIORXIV_PAPER_FRAGMENTOMICS.{md,pdf}` |
| `CITATION.cff` (v0.10.0)                                          | `deepcatch/CITATION.cff` (v2.3.0; v0.10.0 block merged) |
| `.github/workflows/methylation-tests.yml`                        | superseded by `deepcatch`'s existing test infra (no `[methylation]` mirror required) |

## Why archived?

`METHYLATION_PROJECT.md` §9 already concluded this:

> "A separate `rollroyces/deepcatch-methylation` repo is **not** warranted
> at this stage — the project is an extension of DeepCatch, not a separate
> framework. Splitting repos adds CI/maintenance burden without research
> benefit. If the project grows past Phase 4 (the GNN training), we can
> revisit repo splitting."

We did. Phase 0/1/2 work shipped under `deepcatch/src/methylation/` and
this repo has no remaining deliverables.

## What still works here?

- git history intact; SHA references continue to resolve.
- LICENSE (MIT) added on 2026-10-16; until now the repo had no LICENSE
  file -- the MIT license below applies retroactively to all v0.10.0
  and v0.11.0 commits.
- Workflow file `methylation-tests.yml` is left in place but disabled
  by archiving; admins should disable branch-protection rules.

## Where to file a new issue?

- **Methylation** → https://github.com/rollroyces/deepcatch/issues with
  the `methylation` label.
- **cfRNA coupling** (cross-cutting with mRNA therapy design)
  → https://github.com/rollroyces/deepcatch/issues with the `mrnavax` label.

---

*Archived 2026-10-16 as part of the CANCER-SEEK consolidation into `rollroyces/deepcatch` v2.3.0.*

## License

This project is licensed under the MIT License — see the [`LICENSE`](LICENSE) file
in this repository (added 2026-10-16 as part of the consolidation).
