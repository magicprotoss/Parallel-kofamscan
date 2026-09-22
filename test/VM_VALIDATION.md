# VM validation

No local test execution has been performed. The VM was restored on 2026-09-21;
see [VM_TEST_STATUS.md](VM_TEST_STATUS.md) for observed results and active jobs.

Target: `navi@192.168.77.63`, with test runs under `/home/navi/pfamscan_tests`
and the Kofam database at `/mnt/kofamscan-db/`. Copy the updated repository to the
VM after IT restores access. The examples below assume that checkout is at
`/home/navi/pfamscan_tests/Parallel-kofamscan`; adjust `repo` if necessary.

## 1. Environment and regression tests

Run these commands only on the VM, after the environment is available:

```bash
repo=/home/navi/pfamscan_tests/Parallel-kofamscan
conda env update -f "$repo/parallel_kofam_scan/envs/main.yaml"
conda activate parallel_kofam_scan
cd "$repo"
python -m unittest discover -s test -p 'test_*.py' -v
```

The suite includes the bundled `SM1510S13_test.faa` for metadata/collection
validation using synthetic annotation rows. This is not a substitute for the
real Kofam search below. Other cases cover sequence copies across samples,
full headers, case/wrapping normalization, preserved stop residues, long and
colliding names, quoted Kofam definitions, empty thresholds/results, cutoff
boundaries, deterministic ties, and rejection of old UUID intermediates.

## 2. Full workflow with the stock FASTA

Use a fresh, recoverable run directory. Do not use the older shell examples that
delete existing result folders or refer to untracked test inputs.

```bash
stock_run=$(mktemp -d /home/navi/pfamscan_tests/stock.XXXXXX)
/usr/bin/time -v -o "$stock_run/workflow-time.txt" \
  python "$repo/parallel_kofam_scan/pkofamscan" \
    -i "$repo/test/data/SM1510S13_test.faa" \
    -o "$stock_run" -d "$stock_run" \
    -db /mnt/kofamscan-db/ -p 8 -t 4
```

Check that both Parquets and `top-hits.xlsx` are produced, every input sequence
appears in the metadata map, and the Excel top hits agree with the requested
cutoffs and ranking. Compare a small selection with direct Kofam results.
Repeat with a copy of the stock FASTA whose basename exceeds 31 characters and
with duplicated proteins across samples. Confirm the mapping sheet and shared
hashes. Verify gzipped FAA and translated FNA input paths on the VM as well.

## 3. Real 700-genome dataset

Proceed only after the stock run and regression tests pass. Locate the user's
prepared 1,400+ Bakta CDS FAA files on the VM; the exact input directory has not
yet been verified. Record the input manifest, file/sequence counts, database
version, Git revision/diff, environment versions, and worker settings. Use a
fresh run directory under `/home/navi/pfamscan_tests` and the same database.

For collection-specific memory measurements and repeatable comparisons, run the
Snakemake workflow through annotation with `--until exec_annotation` and
`--no-hooks` using a saved config. `--no-hooks` is essential because the normal
success hook deletes the working directory. Once all chunk outputs exist, run
the collector directly under `/usr/bin/time -v`, passing the saved working
directory, a fresh result directory, and the chosen cutoffs:

```bash
/usr/bin/time -v -o "$run_dir/collection-time.txt" \
  python "$repo/parallel_kofam_scan/workflow/scripts/collect_results.py" \
    -wd "$working_dir" -o "$run_dir/collection-results" -e 0.001 -s 100
```

Resolve `run_dir` and `working_dir` from the actual VM run before using this
example. The workflow-wide time measurement from the stock example is not an
isolated collection-memory measurement.

Validate using streamed reads or a disk-backed audit, so verification itself
does not create the old global join:

- Metadata row count equals the total sequence count across the input files;
  full headers and sample identities are preserved.
- Hit row count equals the number of actual HMM hit rows across annotation
  chunks, not the number of hits multiplied by the number of samples.
- Every reported hit hash resolves to input metadata. Shared hashes map to all
  corresponding input occurrences, including copies within one sample.
- Workbook sample count, names, and `sample-id-map` entries agree with the
  input manifest. Samples without hits still have a worksheet.
- Independently ranked hits for selected samples match Excel; there are no
  blank-header rows caused by unrelated samples' hits.
- Record collection peak RSS, elapsed time, output sizes, temporary-disk usage,
  and file-descriptor usage across the 1,400+ worksheets. Check smaller sample
  subsets to establish the scaling behavior. No numeric memory improvement
  should be claimed until these measurements exist.

Old UUID-based annotations cannot be used for these tests with the new
collector. Regenerate sequence hashes and annotation chunks first.
