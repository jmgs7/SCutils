## SCutils v0.10.3

### 🛠️ Bug fixes and improvements.

- Fixed a bug in `CalculateQC()` layer autodetection where `scaled.data` could be incorrectly treated as a log-normalized `data` layer.
- Layer selection in `CalculateQC()` now only accepts true log-normalized layers named like `data`, `data.*`, or `data_*`.
- Added stricter validation for user-provided `layers` in `CalculateQC()` to prevent invalid layer inputs from causing downstream metadata collisions.
- Prevented duplicate `cell.id` row-name failures in subset/scaled Seurat workflows caused by unintended inclusion of non-log layers.
