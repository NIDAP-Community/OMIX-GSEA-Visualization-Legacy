# Canonical OMIX Module Source

## Canonical module

- **Module:** [OMIX GSEA Visualization Legacy](https://github.com/NIDAP-Community/OMIX/tree/main/modules/OMIX-GSEA-Visualization-Legacy)
- **Canonical path:** `modules/OMIX-GSEA-Visualization-Legacy/`
- **Released source reference:** [`d7ff38dd5849de2e63698540f477e60a06933bb9`](https://github.com/NIDAP-Community/OMIX/commit/d7ff38dd5849de2e63698540f477e60a06933bb9)
- **Interface schema:** [schemas/interface.yml](https://github.com/NIDAP-Community/OMIX/blob/main/modules/OMIX-GSEA-Visualization-Legacy/schemas/interface.yml)
- **Module contract:** [OMIX module contract](https://github.com/NIDAP-Community/OMIX/blob/main/docs/module-contract.md)

## Exported scientific files

| Canonical file | Adapter copy | Purpose |
| --- | --- | --- |
| `R/gsea_enrichment_plot.R` | `code/functions/gsea_enrichment_plot.R` | Preserved legacy GSEA ES, RNK, and LE plotting implementation. |

The listed export was verified byte-for-byte against the released source
reference above.

## Adapter-only support code

`code/main.R` performs App Panel translation, recursive workflow discovery,
and result-path handling. It is deployment support code, not part of the
canonical portable CLI.

## Ownership and synchronization

The canonical module owns scientific functions, portable CLI behavior, schemas,
tests, and scientific documentation. This adapter owns deployment UI, input
discovery, result paths, runtime setup, and the platform entry point.

Make reusable changes in the canonical module, update its tests and interface,
record the next released source reference here, and then re-export the listed
file without unreviewed behavioral changes. Validate the adapter with
representative deployment inputs before release.
