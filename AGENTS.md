# Deployment Adapter Agent Instructions

This repository is a deployment adapter for the canonical OMIX module recorded
in [OMIX_MODULE_SOURCE.md](OMIX_MODULE_SOURCE.md). Read that file, this
repository's README, and the canonical [OMIX module contract](https://github.com/NIDAP-Community/OMIX/blob/main/docs/module-contract.md)
before editing.

## Ownership

- The canonical module owns scientific functions, scientific defaults,
  portable CLI behavior, schemas, and tests.
- This repository owns deployment UI, attached-input discovery, result paths,
  runtime configuration, and the platform entry point.
- Preserve the established GSEA plotting behavior. Do not permanently change
  exported scientific functions here; backport scientific or reusable-interface
  changes to the canonical module first.

## Working rules

1. Inspect Git status, the app-panel definition, `OMIX_MODULE_SOURCE.md`, and
   the canonical schema before editing.
2. Keep the three logical inputs unambiguous: MSigDB, filtered GSEA result,
   and the paired DEG table plus sample metadata bundle.
3. Keep explicit DEG/metadata uploads higher priority than workflow discovery.
   If discovery finds more than one suitable candidate, report every path and
   require explicit selection.
4. Preserve stable result names and the input-consistency validation that
   prevents plotting incompatible DEG and GSEA analyses.
5. Use the named pinned runtime. Do not rebuild unrelated shared environments
   or use an unpinned `latest` image.
6. Do not commit input data, generated results, credentials, package caches,
   or package-inventory archives.
7. Validate adapter changes with representative DEG-to-GSEA workflow inputs
   and report deployment validation separately from local tests.

## Release discipline

- Keep the README and `OMIX_MODULE_SOURCE.md` current with the canonical
  module link and source reference.
- Before release, confirm Git is clean, inputs are unambiguous, the app panel
  matches the intended interface, the environment is pinned, and workflow
  handoff works.
