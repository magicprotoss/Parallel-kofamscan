# this script does the following things:
# 1. load and concat per-sample header_map.parquet under the validate_input dir
# validate_input/<sample_id>/header_map.parquet
# header_map_df
# 2. load and concat splitted annotation tsv under the exec_annotation dir
# i.e. exec_annotation/part_*/annotations.tsv
# annotation_df
import os
import sys
import glob
import shutil
import argparse
import pandas as pd
import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        prog='Parallel-kofamscan',
        description="A sub-script for spliting unigenes into chunks"
    )

    parser.add_argument('-wd', '--working_dir', required=True,
                        help="Working dir for the whole workflow, which should be the outdir of the exec_annotation rule")
    parser.add_argument('-o', '--output_dir', required=True,
                        help="Directory to write the final results")
    parser.add_argument('-t', '--threads', required=False, default=1,
                        help="Number of threads to use for polars")
    parser.add_argument('-e', '--min_e_value', default=0.001, type=float,
                        help="Minimum E-value threshold to retain a KO hit during top-hit filtering, default is 0.001")
    parser.add_argument('-s', '--min_score', default=100, type=int,
                        help="Minimum hmmsearch score threshold to retain a KO hit during top-hit filtering, default is 100")

    return parser.parse_args()


def collect_results(working_dir, output_dir, threads, min_e_value, min_score):
    header_map_fps = glob.glob(os.path.join(
        working_dir, "validate_input", "*", "header_map.parquet"))
    # sample_id -> header_map_df
    header_maps = {
        os.path.basename(os.path.dirname(fp)): pd.read_parquet(fp) for fp in header_map_fps
    }
    annotation_fps = glob.glob(os.path.join(
        working_dir, "exec_annotation", "part_*", "annotations.tsv"))
    annotation_dfs = [pd.read_csv(
        fp, sep="\t", header=None, index_col=None, comment="#") for fp in annotation_fps]
    annotation_df = pd.concat(annotation_dfs, ignore_index=True)
    del annotation_dfs

    annotation_df.columns = ["if_score_higher_than_threshold", "gene_name",
                             "KO", "threshold", "hmmsearch_score", "E_value", "KO_definition"]
    # annotation_df_splits sample_id -> annotation_df
    annotation_df_splits = {}
    for sample_id, header_map_df in header_maps.items():
        annotation_df_splits[sample_id] = pd.merge(annotation_df, header_map_df,
                                                   left_on="gene_name", right_on="uuid", how="left")
        annotation_df_splits[sample_id]['sample_id'] = sample_id
    del header_maps, annotation_df

    for sample_id, annotation_df in annotation_df_splits.items():
        annotation_df_splits[sample_id] = annotation_df[['sample_id', "id", "KO", "threshold", "hmmsearch_score",
                                                         "E_value", "KO_definition", "if_score_higher_than_threshold"]]
        annotation_df_splits[sample_id].columns = ['sample_id', "unigene_id", "KO", "threshold", "hmmsearch_score",
                                                   "E_value", "KO_definition", "if_score_higher_than_threshold"]

    # dump the full annotation table to disk
    print(f"Dumping the full annotation table to parquet, this would take a while ...")
    annotation_df_all = pd.concat(
        annotation_df_splits.values(), ignore_index=True)
    annotation_df_all.to_parquet(os.path.join(output_dir, "annotations.parquet"), index=False,
                                 engine='pyarrow',
                                 compression='zstd',
                                 # Zstd (Zstandard) offers highly tunable compression levels ranging from -7 (fastest) to 22 (highest ratio), with default level 3 providing a strong balance for general use. Higher levels (10-19) improve ratios significantly, while ultra-levels (20-22) maximize compression at the cost of memory.
                                 compression_level=19
                                 )
    del annotation_df_all

    # filt res and apply top-hit here
    # yhz: apply cut-off
    # top_hits_dfs: sample_id -> top_hits_df (after cut-off and top-hit filtering)
    print(f"Dumping the per sample top-hits and dumping results to Excel, this would take a while ...")
    top_hits_dfs = {}
    for sample_id, annotation_df in annotation_df_splits.items():
        annotation_df = annotation_df.query(
            f'(E_value <= {min_e_value}) and (hmmsearch_score >= {min_score})')
        # yhz: for each orf, tophit base on e-value and score
        # id is the groupby key, so it is not in the columns to be selected. We use reset_index() to get it back.
        # yhz: KO only
        annotation_df = annotation_df.sort_values(["E_value", 'hmmsearch_score'], ascending=[True, False]).groupby(
            'unigene_id')[["KO"]].first().reset_index()
        top_hits_dfs[sample_id] = annotation_df
    with pd.ExcelWriter(os.path.join(output_dir, "top-hits.xlsx"), engine='openpyxl') as writer:
        for sample_id, top_hits_df in top_hits_dfs.items():
            top_hits_df.to_excel(writer, sheet_name=sample_id, index=False)


def main():
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # not need to overkill here
    if int(args.threads) > 20:
        threads = 20
    else:
        threads = args.threads

    collect_results(args.working_dir, args.output_dir,
                    threads, args.min_e_value, args.min_score)


if __name__ == "__main__":
    sys.exit(main())
