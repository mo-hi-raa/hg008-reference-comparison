# metadata

Small text files that describe the data. They are committed to Git.

Typical files:

- `samples.tsv`: one row per sample (ID, source accession, condition, file names).
- `sources.tsv`: where each raw file is downloaded from, with the provider's checksum.
- `checksums.sha256`: checksums of the raw files you keep, so anyone can check
  they have the same data (`shasum -a 256 -c metadata/checksums.sha256`).
- `README.md` (this file): a data dictionary that explains every column.

## Data dictionary

### sources.tsv

One row per raw file, tab-separated. Read by `scripts/s01_download.sh`.

| Column | Meaning |
|--------|---------|
| `file_name` | Name of the file as saved in `data/raw/` |
| `description` | What the file is, and which arm of the project uses it |
| `source` | Who provides the file |
| `url` | Exact web address the download script uses |
| `provider_checksum` | Where the provider publishes a checksum for the file, or a note if there is none |

### checksums.sha256

One line per file in `data/raw/`: its sha256 code, then its path. Made with
`shasum -a 256` after downloading; check with `shasum -a 256 -c metadata/checksums.sha256`.

## Checks done (2026-11-04)

- The 9 benchmark files matched NIST's md5 codes (`checksums.md5` in the benchmark folder).
- The GRCh38 GIABv3 reference matched GIAB's md5 code (`checksum.md5` in the GIAB references folder).
- T2T-CHM13v2.0 and HG008-N v6.3 have no provider checksum; their sha256 codes are recorded in `checksums.sha256`.
