# Canonical OMIX Module Source

## Canonical module

- **Module:** [OMIX GSEA Visualization Legacy](https://github.com/NIDAP-Community/OMIX/tree/main/modules/OMIX-GSEA-Visualization-Legacy)
- **Canonical path:** `modules/OMIX-GSEA-Visualization-Legacy/`
- **Canonical module version:** `4.0.0`
- **Canonical interface version:** `1`
- **Canonical release tag:** Pending — baseline tag not yet established.
- **Canonical source reference:** [`eee4433cd74d0ce9bd1daba1571ca6312cba7ea2`](https://github.com/NIDAP-Community/OMIX/commit/eee4433cd74d0ce9bd1daba1571ca6312cba7ea2)
- **Interface schema:** [schemas/interface.yml](https://github.com/NIDAP-Community/OMIX/blob/main/modules/OMIX-GSEA-Visualization-Legacy/schemas/interface.yml)
- **Module contract:** [OMIX module contract](https://github.com/NIDAP-Community/OMIX/blob/main/docs/module-contract.md)

## Exported scientific files

| Canonical file | Adapter copy | SHA-256 | Purpose |
| --- | --- | --- | --- |
| `R/gsea_enrichment_plot.R` | `code/functions/gsea_enrichment_plot.R` | `5fdae88b3968443a492b654d1d0810292fea4cb8ed70242acb4de17f1ffe4e5a` | Preserved legacy GSEA ES, RNK, and LE plotting implementation. |

The listed export was verified byte-for-byte against the canonical source
reference above. The canonical `schemas/interface.yml` used for App Panel and
adapter-contract reconciliation has SHA-256
`f322dcfbb48616aba755cd51c94772e83a24f61998fc75b69562732fc36342f4`.

## Platform translations

<!-- omix-adapter-contract: {"hidden_inputs":[{"canonical":"msigdb_database","binding":"dataset:eb51e72e-bcb9-400d-bf94-fc8cad85c7f1 mount=msigdb"},{"canonical":"gsea_filter_results","binding":"dataset:e5c4e095-f2e5-48ac-8148-d7fdb8673e05 mount=gsea_filter_results"}]} -->

- The required canonical `msigdb_database` input is a platform-managed data
  asset mounted below `/data/msigdb`; it is intentionally hidden from the
  named-parameter panel. The adapter recursively requires exactly one
  `msigdb*.csv` or `msigdb*.rds` candidate and rejects ambiguity.
- The required canonical `gsea_filter_results` input is a platform-managed
  GSEA Filters result mounted below `/data/gsea_filter_results`; it is
  intentionally hidden from the named-parameter panel. The adapter prefers
  the stable `filtered_gsea_results.csv` output and accepts the historical RDS
  equivalent only when no CSV is present.
- The canonical `deg_table` and `sample_metadata` inputs are exposed as an
  optional explicit pair. When both selectors are blank, the adapter resolves
  them together from the single combined DEG Analysis asset mounted below
  `/data/deg-training`. Supplying only one file is an error; supplying both
  makes the pair authoritative.
- Canonical internal `output_dir` is translated to Code Ocean's `/results` and
  is not shown in the App Panel.
- All other public and advanced canonical controls use their canonical CLI
  names, types, allowed values, order, and scientific defaults. There are no
  scientific aliases or hidden user-settable parameters.
- The retired `code/GSEA_Visualization_Local_v1.R` was an unreferenced,
  divergent historical copy. Git history preserves it; no launcher, entry
  point, documentation, or test sourced it. It is removed so the only active
  scientific implementation is the managed canonical export under
  `code/functions/`.

## Adapter release record

| Field | Recorded value |
| --- | --- |
| Adapter version | Pending — baseline tag not yet established. |
| Adapter release tag | Pending. |
| Platform release | Pending validation record. |
| Runtime tag | `codeocean/omix-r-pathway:r4.4.3-bioconductor3.20-v1` (from `.codeocean/environment.json`). |
| Runtime digest | Pending — no immutable digest has been validated or recorded. |
| Named-parameter platform run | Pending — GitHub/local checks are not Code Ocean validation. |
| Syncweaver lock | Pending — no generated `.syncweaver-lock.json` is recorded. |

See the [OMIX versioning and release policy](https://github.com/NIDAP-Community/OMIX/blob/main/docs/versioning-and-releases.md). The source commit, adapter tag, platform release, and runtime identity are separate records.

## Adapter-only support code

`code/main.R` performs App Panel translation, recursive workflow discovery,
paired DEG/metadata resolution, and result-path handling. It is deployment
support code, not part of the canonical portable CLI.

## Code Ocean validation plan

Before an adapter tag or release claim, run the capsule once with its three
default data assets and once with a workflow-connected GSEA Filters result plus
combined DEG result. Confirm the stable PDF, RunningES CSV, bubble-plot PNGs,
manifest are present, confirm the input-consistency summary in the run log, and visually inspect at
least one ES+RNK+LE page and one bubble plot. Repeat with the explicit DEG and
metadata selectors, and verify that a partial explicit pair and every
ambiguous asset configuration fail with candidate paths. Record the capsule
release identifier and the resolved runtime digest here. These external checks
remain pending; local and GitHub evidence does not satisfy them.

## Ownership and synchronization

The canonical module owns scientific functions, portable CLI behavior, schemas,
tests, and scientific documentation. This adapter owns deployment UI, input
discovery, result paths, runtime setup, and the platform entry point.

Make reusable changes in the canonical module, update its tests and interface,
record the next canonical version and immutable source reference here, and then re-export the listed
file without unreviewed behavioral changes. Validate the adapter with
representative deployment inputs before release.
