# this script does the following things:
# 1. load and concat per-sample header_map.parquet under the validate_input dir
# validate_input/<sample_id>/header_map.parquet
# header_map_df
# 2. load and concat splitted annotation tsv under the exec_annotation dir
# i.e. exec_annotation/part_*/annotations.tsv
# annotation_df
#
# DuckDB rewrite: uses DuckDB's columnar SQL engine to process data with minimal RAM
# instead of loading all DataFrames into pandas memory simultaneously.
import os
import sys
import argparse
import duckdb
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        prog='Parallel-kofamscan',
        description="A sub-script for spliting unigenes into chunks"
    )

    parser.add_argument('-wd', '--working_dir', required=True,
                        help="Working dir for the whole workflow, which should be the outdir of the exec_annotation rule")
    parser.add_argument('-o', '--output_dir', required=True,
                        help="Directory to write the final results")
    parser.add_argument('-t', '--threads', required=False, default=1, type=int,
                        help="Number of threads to use for DuckDB")
    parser.add_argument('-e', '--min_e_value', default=0.001, type=float,
                        help="Minimum E-value threshold to retain a KO hit during top-hit filtering, default is 0.001")
    parser.add_argument('-s', '--min_score', default=100, type=int,
                        help="Minimum hmmsearch score threshold to retain a KO hit during top-hit filtering, default is 100")

    return parser.parse_args()


def collect_results(working_dir, output_dir, threads, min_e_value, min_score):
    header_maps_glob = os.path.join(
        working_dir, "validate_input", "*", "header_map.parquet")
    annotations_glob = os.path.join(
        working_dir, "exec_annotation", "part_*", "annotations.tsv")

    con = duckdb.connect()
    if threads > 1:
        con.execute(f"SET threads TO {threads}")

    # Register header_map parquet files with filename to extract sample_id.
    # Each file path is: .../validate_input/<sample_id>/header_map.parquet
    con.execute(f"""
        CREATE OR REPLACE TABLE header_maps AS
        SELECT
            regexp_extract(filename, 'validate_input/([^/\\\\]+)/header_map\\.parquet', 1) AS sample_id,
            id,
            uuid
        FROM read_parquet('{header_maps_glob}', filename=true)
    """)

    # Register all annotation TSV files.
    # Columns: '*' marker, gene_name, KO, threshold, score, E_value, definition
    # kofamscan output starts with '#' comment lines; comment='#' skips them.
    con.execute(f"""
        CREATE OR REPLACE TABLE annotations AS
        SELECT column0 AS if_score_higher_than_threshold,
               column1 AS gene_name,
               column2 AS KO,
               CAST(column3 AS DOUBLE) AS threshold,
               CAST(column4 AS DOUBLE) AS hmmsearch_score,
               CAST(column5 AS DOUBLE) AS E_value,
               column6 AS KO_definition
        FROM read_csv('{annotations_glob}', header=false, delim='\\t',
                      all_varchar=true, comment='#', quote='"', escape='"')
    """)

    # Build the full annotation table matching original pandas behaviour:
    # for each sample, LEFT JOIN all annotations with that sample's header_map.
    # This produces num_samples × num_annotations rows; unigene_id is filled only
    # when gene_name matches the sample's uuid.
    con.execute("""
        CREATE OR REPLACE TABLE full_annotations AS
        SELECT
            s.sample_id,
            hm.id AS unigene_id,
            a.KO,
            a.threshold,
            a.hmmsearch_score,
            a.E_value,
            a.KO_definition,
            a.if_score_higher_than_threshold
        FROM annotations a
        CROSS JOIN (SELECT DISTINCT sample_id FROM header_maps) s
        LEFT JOIN header_maps hm
            ON a.gene_name = hm.uuid AND s.sample_id = hm.sample_id
    """)

    # Dump full annotation table to parquet
    print("Dumping the full annotation table to parquet, this would take a while ...")
    con.execute(f"""
        COPY (SELECT * FROM full_annotations)
        TO '{os.path.join(output_dir, "annotations.parquet")}'
        (FORMAT PARQUET, COMPRESSION ZSTD, COMPRESSION_LEVEL 19)
    """)

    # Per-sample top-hit filtering
    print("Dumping the per sample top-hits and dumping results to Excel, this would take a while ...")
    top_hits_sql = f"""
        CREATE OR REPLACE TABLE top_hits AS
        SELECT sample_id, unigene_id, KO
        FROM (
            SELECT
                sample_id,
                unigene_id,
                KO,
                ROW_NUMBER() OVER (
                    PARTITION BY sample_id, unigene_id
                    ORDER BY E_value ASC, hmmsearch_score DESC
                ) AS rn
            FROM full_annotations
            WHERE E_value <= {min_e_value} AND hmmsearch_score >= {min_score}
        ) sub
        WHERE rn = 1
    """
    con.execute(top_hits_sql)

    samples = con.execute(
        "SELECT DISTINCT sample_id FROM top_hits ORDER BY sample_id").fetchall()

    with pd.ExcelWriter(os.path.join(output_dir, "top-hits.xlsx"), engine='openpyxl') as writer:
        for (sample_id,) in samples:
            df = con.execute(
                "SELECT unigene_id, KO FROM top_hits WHERE sample_id = ?",
                [sample_id]
            ).df()
            df.to_excel(writer, sheet_name=sample_id, index=False)

    con.close()


def main():
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # not need to overkill here
    if int(args.threads) > 32:
        threads = 32
    else:
        threads = args.threads

    collect_results(args.working_dir, args.output_dir,
                    threads, args.min_e_value, args.min_score)


if __name__ == "__main__":
    sys.exit(main())
