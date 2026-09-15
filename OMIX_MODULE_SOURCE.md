# Canonical OMIX Module Source

## Canonical module

- **Module:** [OMIX GSEA Visualization Legacy](https://github.com/NIDAP-Community/OMIX/tree/main/modules/OMIX-GSEA-Visualization-Legacy)
- **Canonical path:** `modules/OMIX-GSEA-Visualization-Legacy/`
- **Canonical module version:** `4.0.0`
- **Canonical interface version:** `1`
- **Canonical release tag:** Pending — baseline tag not yet established.
- **Canonical source reference:** [`b65eeb3c231eb9dfbdbc09326ab03c23bbe84c8f`](https://github.com/NIDAP-Community/OMIX/commit/b65eeb3c231eb9dfbdbc09326ab03c23bbe84c8f)
- **Interface schema:** [schemas/interface.yml](https://github.com/NIDAP-Community/OMIX/blob/main/modules/OMIX-GSEA-Visualization-Legacy/schemas/interface.yml)
- **Module contract:** [OMIX module contract](https://github.com/NIDAP-Community/OMIX/blob/main/docs/module-contract.md)

## Exported scientific files

| Canonical file | Adapter copy | Purpose |
| --- | --- | --- |
| `R/gsea_enrichment_plot.R` | `code/functions/gsea_enrichment_plot.R` | Preserved legacy GSEA ES, RNK, and LE plotting implementation. |

The listed export was verified byte-for-byte against the canonical source
reference above.

## Adapter release record

| Field | Recorded value |
| --- | --- |
| Adapter version | Pending — baseline tag not yet established. |
| Adapter release tag | Pending. |
| Platform release | Pending validation record. |
| Runtime identity | Not yet recorded as an immutable image digest or lockfile reference. |

See the [OMIX versioning and release policy](https://github.com/NIDAP-Community/OMIX/blob/main/docs/versioning-and-releases.md). The source commit, adapter tag, platform release, and runtime identity are separate records.

## Adapter-only support code

`code/main.R` performs App Panel translation, recursive workflow discovery,
and result-path handling. It is deployment support code, not part of the
canonical portable CLI.

## Ownership and synchronization

The canonical module owns scientific functions, portable CLI behavior, schemas,
tests, and scientific documentation. This adapter owns deployment UI, input
discovery, result paths, runtime setup, and the platform entry point.

Make reusable changes in the canonical module, update its tests and interface,
record the next canonical version and immutable source reference here, and then re-export the listed
file without unreviewed behavioral changes. Validate the adapter with
representative deployment inputs before release.
