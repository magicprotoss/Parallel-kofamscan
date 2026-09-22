# VM validation results — verified 2026-09-22

The stock workflow and full dataset both completed successfully and passed
their output audits. Full-data audit finished at 05:20:08 Asia/Shanghai on
2026-09-22; checked via SSH at 09:36. Controller stage is `complete`, with exit
code 0. An attempt to pause the completed monitor returned "Automation does not
exist in the app"; no replacement monitor was created.
All test execution is on `navi@192.168.77.63`; none was performed locally.

| Check | Observed status |
| --- | --- |
| Regression suite | 15 tests passed in 3.484 seconds; exit code 0 |
| Real SeqKit integration | Passed FNA translation, gzipped FAA normalization, shared SHA-256 IDs, and deduplication/shuffle/splitting; exit code 0 |
| Stock workflow | Passed; 2 samples, 3,246 sequence occurrences, 1,623 unique hashes, 41,235 HMM hit rows; both Excel sheets audited |
| Full dataset | Passed; 1,402 FAA files, 2,620,226 sequence occurrences, 1,986,303 unique hashes, 54,147,991 HMM hit rows |

## Performance and validation scope

| Measurement | Stock run | Full dataset |
| --- | --- | --- |
| Workflow elapsed time | 13m 30.75s | 11h 06m 11s |
| Collection elapsed time | 0.88s | 633.80s (10m 33.80s) |
| Collection sampled peak RSS | 144.41 MiB | 258.39 MiB |

Collection measurements come from Snakemake's benchmark, which divides bytes
by 1024 squared for memory values. They are sampled peaks, not a hard memory
ceiling. No comparable run of the old collector was performed, so no measured
fold-reduction or speedup over the old implementation is claimed.

The full audit checked every metadata row against the original FAA inputs,
verified that all hit hashes exist in metadata, and confirmed that the raw
annotation chunks and Parquet contain exactly the same 54,147,991 hit rows.
All sample-sheet mappings passed. Excel top-hit identities were fully checked
for five samples, including comparison with raw Kofam results; the remaining
1,397 sheets were not exhaustively checked cell by cell. Stock Excel identities
were checked for both sheets; its raw chunks were removed by the normal cleanup
hook, so raw/Parquet row-count comparison was not available for that run.

Persistent full-run outputs are under
`/home/navi/pfamscan_tests/700genome_test/run-20260921/KEGG_annotations/final_results/`:

| Output | Bytes |
| --- | ---: |
| `hmm-hits.parquet` | 1,287,278,194 |
| `sequence-metadata-map.parquet` | 121,123,450 |
| `top-hits.xlsx` | 31,838,416 |

## Job records and files

- Repository: `/home/navi/pfamscan_tests/Parallel-kofamscan`
- Regression evidence: `/home/navi/pfamscan_tests/repo_simple_test/regression-20260921/`
- SeqKit integration log: `/home/navi/pfamscan_tests/repo_simple_test/preprocess-20260921.log`
- Completed stock tmux session: `kofam-stock-20260921`
- Stock run: `/home/navi/pfamscan_tests/repo_simple_test/stock-20260921/`
- Completed continuation tmux session: `kofam-validation-20260921`
- Full run: `/home/navi/pfamscan_tests/700genome_test/run-20260921/`
- Full input directory: `/home/navi/pfamscan_tests/700genome_test/input_bakta_faas/`

The continuation is `test/run_vm_jobs.sh`. It stops on a failed stock job or
audit, then caches the supplied database and runs the full workflow with
`test/vm_workflow.py`. That harness uses the normal CLI's option validation and
Snakemake workflow, but disables cleanup hooks to retain raw annotation chunks
for comparison. `test/verify_vm_results.py` compares every metadata occurrence
with the original FAA, checks that every hit hash exists in metadata, compares
raw/Parquet hit counts when raw chunks are retained, and independently checks
Excel top hits for up to five samples. It validates all sample-sheet mappings.

Read `stage`, `controller.log`, and `controller.exit-code` in the full-run
directory for current status. `stage=complete` and exit code 0 indicate both
the workflow and its audit succeeded. An absent exit-code file means the job
may still be running; confirm via tmux/process state. Do not launch a second
copy of this dated job over existing output directories.

Collection-only timing and sampled peak RSS are in each workflow's
`KEGG_annotations/logs/collect_results.benchmark.tsv`. Workflow-wide GNU time
measurements are in `workflow-time.txt`. Audit outcomes are in each run's
`verification.json`; the full audit took 5m 55.34s and exited with code 0.

## Storage

The previously blank 1 TB `/dev/sdb` was initialized as ext4 and mounted at
`/mnt/sdb`, with scratch directories owned by `navi` beneath
`/mnt/sdb/pfamscan-tests`. Full-run intermediates, temporary files, and a cache
of `/mnt/kofamscan-db/` use this non-persistent drive. Source input files, logs,
configuration, audit reports, and final outputs remain on the persistent root
disk under `/home/navi/pfamscan_tests`. The scratch mount is not added to fstab;
a VM reset may require reinitialization and re-running affected stages.
