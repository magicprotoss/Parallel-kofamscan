# Parallel-kofamscan

A snakemake workflow to split large unigenes.faa into chucks in order to prevent [kofam_scan](https://github.com/takaram/kofam_scan) from stuck at formatting results using ruby. This workflow also filters the results down to top-hits and saves it to a excel sheet for the user.

## Installation

1.  clone the repo

``` bash
git clone https://github.com/magicprotoss/Parallel-kofamscan
```

2.  install create 2 conda envs for the main shell as well as dependencies

``` bash
# main env, named parallel_kofam_scan
conda env create -f Parallel-kofamscan/parallel_kofam_scan/envs/main.yaml
# dependency env, named parallel-kofamscan.dependency.kofamscan
conda env create -f Parallel-kofamscan/parallel_kofam_scan/envs/kofamscan.yaml
```

3.  move the repo folder under the main env, soft link the scripts to the bin sub-dir, and add exec permission

``` bash
mv Parallel-kofamscan <your-path-to-miniconda3>/envs/parallel-kofamscan/
cd <your-path-to-miniconda3>/envs/parallel-kofamscan/Parallel-kofamscan
ln -s parallel_kofam_scan/pkofamscan ../bin/ && chmod u+x ../bin/pkofamscan
ln -s parallel_kofam_scan/workflow ../bin && chmod u+x ../bin/workflow/scripts/*.py
```

To use it, simply activate the conda env and run pkofamscan

``` bash
conda activate parallel_kofam_scan

pkofamscan --help
usage: Parallel-kofamscan [-h] -i INPUT_PATH [INPUT_PATH ...] -o OUTPUT_PATH -db PATH_TO_KEGG_DATABASE
                          [-d WORKING_DIR] [-p WORKERS] [-t THREADS_PER_WORKER] [-c CHUNK_SIZE] [-f FORCE]
                          [-e MIN_E_VALUE] [-s MIN_SCORE]

A snakemake workflow to split unigenes into chucks in order to reduce kofamscan's runtime

options:
  -h, --help            show this help message and exit
  -i INPUT_PATH [INPUT_PATH ...], --input_path INPUT_PATH [INPUT_PATH ...]
                        Unigenes to annotate, filenames will be converted to sample-ids, gzipped files are
                        supported as well
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        Directory to save the annotation results, the output would be stored in an sub-dir
                        called KEGG_annotations
  -db PATH_TO_KEGG_DATABASE, --path_to_KEGG_database PATH_TO_KEGG_DATABASE
                        Path to the kofamscan database, which contains a 'ko_list' file and a 'profiles'
                        directory
  -d WORKING_DIR, --working_dir WORKING_DIR
                        Directory to store intermediate results, default is
                        /tmp/parallel_kofamscan/<your_uid>/<hash_of_output_dir>_<hash_of_current_datetime>,
                        if set, the temp dir would be <working_dir>/parallel_kofamscan/<your_uid>/<hash_of_ou
                        tput_dir>_<hash_of_current_datetime>
  -p WORKERS, --workers WORKERS
                        Number of cores to use for when running the workflow, default is 0, which indicates
                        the workflow will use the num of cores on the machine - 1, or the num of chuncks *
                        16, whichever is smaller
  -t THREADS_PER_WORKER, --threads_per_worker THREADS_PER_WORKER
                        Number of theards to use when running the exec_annotation cmd, default is 8
  -c CHUNK_SIZE, --chunk_size CHUNK_SIZE
                        Number of unigenes per chunk, default is 100000
  -f FORCE, --force FORCE
                        whether to overwrite existing result, default is no
  -e MIN_E_VALUE, --min_e_value MIN_E_VALUE
                        Maximum E-value threshold to retain a KO hit during top-hit filtering, default is
                        0.001
  -s MIN_SCORE, --min_score MIN_SCORE
                        Minimum hmmsearch score threshold to retain a KO hit during top-hit filtering,
                        default is 100
```

## Results

Results are written under `<output_path>/KEGG_annotations/final_results/`:

| File | Contents |
| --- | --- |
| `hmm-hits.parquet` | One row per reported HMM hit: `sequence_hash`, `KO`, `threshold`, `hmmsearch_score`, `E_value`, `KO_definition`, `if_score_higher_than_threshold`. All hits are retained, including those below the top-hit cutoffs. A hash can have multiple hits. Missing KO thresholds are null. |
| `sequence-metadata-map.parquet` | One row per input sequence occurrence: `sample_id`, `sequence_header`, `sequence_hash`. Full headers (without `>`) and sequences without hits are preserved, including duplicate sequences within or across samples. |
| `top-hits.xlsx` | One worksheet per sample, with `unigene_id` (the full sequence header) and `KO`. Samples with no qualifying hits have a header-only sheet. |

`sequence_hash` is the lowercase hexadecimal SHA-256 of the protein sequence
actually emitted for searching: sequence whitespace is removed and residues are
uppercased. Terminal and internal `*` characters are preserved. FNA inputs are
translated with genetic code 11 before hashing. Files ending in `.faa` or
`.faa.gz` are treated as proteins; other FASTA inputs are classified by content.
Identical canonical proteins receive the same hash and are searched once.

Top-hit filtering uses `E_value <= --min_e_value` and
`hmmsearch_score >= --min_score`. Among qualifying hits for a hash, the lowest
E-value wins, followed by the highest score, then the lexically smallest KO for
an exact tie. Every sequence occurrence with that hash receives the same KO.
Duplicate full headers remain separate occurrences, even when their sequences
differ. The Kofam `*` marker is retained in the hit table but is not an additional
top-hit filter.

Before exporting Excel, all sample names are checked against the 31-character
worksheet limit, invalid characters, reserved names, and case-insensitive
collisions. If any name fails, **all** sample sheets use `FAA1`, `FNA2`, etc.,
with the prefix indicating the input type and a global sequence number assigned
in sorted sample-ID order. A `sample-id-map` sheet records `sample_id` and
`sheet_name` for every sample. Otherwise, the original sample IDs are retained.
There are no literal brackets or slashes in the generated sheet names.

The collector writes Parquet in batches of 8,192 rows. It keeps only the best
qualifying hit per hash in a temporary on-disk SQLite index with an 8 MiB page
cache, and uses batched lookups to stream sample sheets through openpyxl's
write-only mode. It does not build a sample-by-hit table or load all workbook
cells into memory. Temporary SQLite and staged output files are placed under
the result directory; openpyxl worksheet XML uses the system temporary
directory (respect `TMPDIR` when selecting storage on the VM). The `-t` collector
option is retained for compatibility; collection uses one worker.

### Migration

The two Parquet files replace `annotations.parquet`. Readers should use
`sequence_hash` to associate hits with sequence metadata, restricting to the
samples/hashes needed rather than reconstructing a joined table for all samples.
Old UUID-based working directories cannot be resumed with this version: rerun
validation and annotation in a fresh output directory. The collector rejects
legacy header-map schemas. The historical `test/merge_results.ipynb` uses the old
intermediate format and is not a reader for the new results.

The main environment no longer needs DuckDB; `lxml` is included for streaming
Excel export. Update the environment from `parallel_kofam_scan/envs/main.yaml`
on the VM before validation.

## Testing

Regression tests and the VM validation procedure are in
[`test/test_collect_results.py`](test/test_collect_results.py) and
[`test/VM_VALIDATION.md`](test/VM_VALIDATION.md). Tests have **not been run locally**.
The VM regression suite, real SeqKit preprocessing check, stock workflow, and
1,402-file dataset passed their checks. Full collection processed 54,147,991
hits in 633.80 seconds with 258.39 MiB sampled peak RSS. See
[`test/VM_TEST_STATUS.md`](test/VM_TEST_STATUS.md) for results, audit coverage,
measurement limitations, and output locations (verified 2026-09-22).
