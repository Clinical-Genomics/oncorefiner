#!/usr/bin/env python3

import math
import click
import pandas as pd

__version__ = "1.0"

DEFAULT_SPACING = {
    "o": 50_000,
    "a": 20_000,
    "b": 10_000,
    "c": 5_000,
    "d": 1_000,
}


def chrom_to_output_name(chromosome: str, level: str) -> str:
    chrom = str(chromosome).replace("chr", "")
    return f"{level}_{chrom}"


def segment_to_points(start: int, end: int, spacing: int) -> list[int]:
    """
    Split a segment into new segments (chunks) based on the level from the DEFAULT_SPACING dictionary,
    and return the midpoint of each chunk.

    Always returns at least one point.

    """
    length = end - start + 1
    n_points = max(1, math.ceil(length / spacing))

    points = []

    for i in range(n_points):
        chunk_start = start + int(i * length / n_points)
        chunk_end = start + int((i + 1) * length / n_points) - 1
        midpoint = (chunk_start + chunk_end) // 2
        points.append(midpoint)

    return points


def write_segment_zoom_file(
    df: pd.DataFrame,
    output_path: str,
    levels: dict[str, int],
) -> None:
    """
    Write a zoom-level file for GENS from the input DataFrame, converting to the specified levels.
    The input DataFrame is expected to have the following columns:
    - chrom: chromosome name
    - start.pos: start position of the segment
    - end.pos: end position of the segment
    - mean: mean value of the segment
    """

    rows = []

    for level, spacing in levels.items():
        for _, row in df.iterrows():
            chrom = row["chrom"]
            start = int(row["start.pos"])
            end = int(row["end.pos"])
            value = float(row["mean"])

            for pos in segment_to_points(start, end, spacing):
                rows.append(
                    [
                        chrom_to_output_name(chrom, level),
                        pos,
                        pos + 1,
                        value,
                    ]
                )

    out_df = pd.DataFrame(rows)
    out_df.to_csv(output_path, sep="\t", header=False, index=False)


@click.command()
@click.option(
    "--input-file",
    required=True,
    type=click.Path(exists=True),
    help="Input segment PCF file from Cobalt. Can be plain text or .gz.",
)
@click.option(
    "--output-file",
    required=True,
    help="Output BED-like zoom file for GENS.",
)

@click.version_option(__version__, prog_name="prepare_cobalt_for_gens")

def main(input_file: str, output_file: str) -> None:
    """Convert segment mean file into GENS zoom-level data."""

    df = pd.read_csv(
        input_file,
        sep=r"\s+",
        engine="python",
        compression="infer",
    )

    required_columns = {
        "chrom",
        "start.pos",
        "end.pos",
        "mean",
    }

    missing = required_columns - set(df.columns)

    if missing:
        raise click.ClickException(
            f"Missing required columns: {', '.join(sorted(missing))}"
        )

    df = df.sort_values(["chrom", "start.pos", "end.pos"])

    write_segment_zoom_file(
        df=df,
        output_path=output_file,
        levels=DEFAULT_SPACING,
    )

    click.echo(f"Wrote segment zoom file: {output_file}")

    click.echo(f"prepare_cobalt_for_gens {__version__}")


if __name__ == "__main__":
    main()
