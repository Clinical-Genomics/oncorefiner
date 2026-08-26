#!/usr/bin/env python3

import click
import pandas as pd

### tumor and normal baf takes the analysis_type parameter to know if normal analysis is needed.

__version__ = "1.0"

DEFAULT_LEVELS = {
    "o": 135,
    "a": 30,
    "b": 8,
    "c": 3,
    "d": 1,
}

# Updating chromosome names to be compatible with GENS. For example, chr1 becomes o_1, chrX becomes o_X, etc.
def chrom_to_output_name(chromosome: str, level: str) -> str:
    chrom = str(chromosome).replace("chr", "")
    return f"{level}_{chrom}"

# Generate a BAF file for GENS from the input DataFrame, converting to the specified levels.
def generate_baf_file_for_gens(
    df: pd.DataFrame,
    baf_column: str,
    output_path: str,
    levels: dict,
) -> None:
    """
    Generate a BAF file for GENS from the input DataFrame, converting to the specified levels.

    The input DataFrame is expected to have the following columns:
    - chromosome: chromosome name
    - position: genomic position
    - baf_column: BAF value column specified by the baf_column parameter (ex normalBAF or tumorBAF)
    """
    rows = []

    for level, step in levels.items():
        for chrom, chrom_df in df.groupby("chromosome", sort=False):
            sampled = chrom_df.iloc[::step]

            for _, row in sampled.iterrows():
                pos = int(row["position"])

                rows.append(
                    [
                        chrom_to_output_name(chrom, level),
                        pos,
                        pos + 1,
                        f"{float(row[baf_column]):.4f}",
                    ]
                )

    pd.DataFrame(rows).to_csv(
        output_path,
        sep="\t",
        header=False,
        index=False,
    )


@click.command()
@click.option(
    "--input-file",
    required=True,
    type=click.Path(exists=True),
    help="Input BAF TSV file",
)
@click.option(
    "--analysis-type",
    show_default=True,
    type=click.Choice(["tumor_only", "tumor_normal"]),
    help="Analysis type (e.g., tumor_only or tumor_normal)",
)
@click.option(
    "--output-file-prefix",
    default="*.bed",
    show_default=True,
    help="Output file prefix for BAF",
)

@click.version_option(__version__, prog_name="prepare_amber_for_gens")

def main(
    input_file: click.Path,
    output_file_prefix: click.STRING,
    analysis_type: click.STRING,
) -> None:
    """Convert BAF table into zoom-level files."""

    df = pd.read_csv(
        input_file,
        sep=r"\s+",
        engine="python",
    )

    required_columns = {"chromosome", "position"}
    missing = required_columns - set(df.columns)

    if missing:
        raise click.ClickException(
            f"Missing required columns: {', '.join(sorted(missing))}"
        )

    generate_baf_file_for_gens(
        df=df,
        baf_column="tumorBAF",
        output_path=f"{output_file_prefix}.tumor.bed",
        levels=DEFAULT_LEVELS,
    )
    click.echo(f"Wrote tumor output: {output_file_prefix}.tumor.bed")

    if "normal" in analysis_type:
        generate_baf_file_for_gens(
            df=df,
            baf_column="normalBAF",
            output_path=f"{output_file_prefix}.normal.bed",
            levels=DEFAULT_LEVELS,
        )
        click.echo(f"Wrote normal output: {output_file_prefix}.normal.bed")

    click.echo(f"prepare_amber_for_gens {__version__}")


if __name__ == "__main__":
    main()
