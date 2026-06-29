#!/usr/bin/env python3

import click
import pandas as pd

### tumor and normal baf take tumor/normal parameter to know which file to check.

DEFAULT_LEVELS = {
    "o": 135,
    "a": 30,
    "b": 8,
    "c": 3,
    "d": 1,
}


def chrom_to_output_name(chromosome: str, level: str) -> str:
    chrom = str(chromosome).replace("chr", "")
    return f"{level}_{chrom}"


def write_baf_file(
    df: pd.DataFrame,
    baf_column: str,
    output_path: str,
    levels: dict,
) -> None:
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
    "--tumor-output",
    default="tumor_baf_zoom.tsv",
    show_default=True,
    help="Output file for tumorBAF",
)
@click.option(
    "--normal-output",
    default="normal_baf_zoom.tsv",
    show_default=True,
    help="Output file for normalBAF",
)
def main(
    input_file: str,
    tumor_output: str,
    normal_output: str,
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

    if "tumorBAF" in df.columns:
        write_baf_file(
            df=df,
            baf_column="tumorBAF",
            output_path=tumor_output,
            levels=DEFAULT_LEVELS,
        )
        click.echo(f"Wrote tumor output: {tumor_output}")
    else:
        click.echo("No tumorBAF column found.")

    if "normalBAF" in df.columns:
        write_baf_file(
            df=df,
            baf_column="normalBAF",
            output_path=normal_output,
            levels=DEFAULT_LEVELS,
        )
        click.echo(f"Wrote normal output: {normal_output}")
    else:
        click.echo("No normalBAF column found.")


if __name__ == "__main__":
    main()
