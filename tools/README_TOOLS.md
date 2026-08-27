# Tools

## Integrity check

Copy these two files into the completed `result` folder and run
`VERIFY_RESULT.bat`:

- `VERIFY_RESULT.bat`
- `verify_result_integrity.ps1`

The completed production run returned `OVERALL: PASS`.

## Zenodo archive builder

After the integrity pass, copy

- `build_zenodo_archives.py`

into the same completed `result` folder and run

```text
python build_zenodo_archives.py
```

It creates `result/zenodo_upload/` with five logically partitioned prime-field
archives and one prime-power archive, plus the run metadata/manifests.

The archives use ZIP "store" mode because the certificate CSV files are already
gzip-compressed. No raw result file is altered.
