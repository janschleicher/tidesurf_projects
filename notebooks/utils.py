from itertools import product
from typing import Tuple, Union

import numpy as np
import pandas as pd


def overlapping(
    chr_1: Union[int, str],
    start_1: int,
    end_1: int,
    strand_1: str,
    chr_2: Union[int, str],
    start_2: int,
    end_2: int,
    strand_2: str,
    opposing: bool = True,
) -> Tuple[bool, int, float]:
    """
    Check if two genes overlap. Optionally return the length of the overlap.

    Parameters
    ----------
    chr_1 : Union[int, str]
        Chromosome of the first gene.
    start_1 : int
        Start position of the first gene.
    end_1 : int
        End position of the first gene.
    strand_1 : str
        Strand of the first gene.
    chr_2 : Union[int, str]
        Chromosome of the second gene.
    start_2 : int
        Start position of the second gene.
    end_2 : int
        End position of the second gene.
    strand_2 : str
        Strand of the second gene.
    opposing : bool, optional
        Consider genes on opposing strands. By default True.

    Returns
    -------
    Tuple[bool, int, float
        Whether the two genes are overlapping, the number of overlapping
        bases, and relative overlap (compared to the shorter one of the
        two genes).
    """
    if (
        chr_1 != chr_2
        or (opposing and strand_1 == strand_2)
        or (not opposing and strand_1 != strand_2)
    ):
        return False, 0, 0.0
    assert start_1 < end_1 and start_2 < end_2, "Start must be less than end"
    if start_1 > end_2 or start_2 > end_1:
        return False, 0, 0.0

    # Compute overlap interval (handles partial overlaps and containment)
    overlap_start = max(start_1, start_2)
    overlap_end = min(end_1, end_2)
    overlap_bases = max(0, overlap_end - overlap_start + 1)
    overlap_relative = overlap_bases / (max(end_1, end_2) - min(start_1, start_2) + 1)

    return True, overlap_bases, overlap_relative


def compute_gene_overlaps(gene_positions: pd.DataFrame) -> pd.DataFrame:
    """
    Compute overlap between pairs of genes on opposite strands

    Parameters
    ----------
    gene_positions : pd.DataFrame
        Dataframe containing the positions of all genes. Must contain
        the following columns: Chromosome, Start, End, Strand

    Returns
    -------
    pd.DataFrame
        All pairs of overlapping genes with their absolute and relative
        overlaps. Columns are suffixed by `_x` and `_y` for the plus and
        minus strand genes, respectively.
    """
    # Get all gene combinations of genes on opposite strands
    gene_pairs = []
    for chr_name, chr_df in gene_positions.groupby("Chromosome"):
        genes_plus = chr_df[chr_df.Strand == "+"].index
        genes_minus = chr_df[chr_df.Strand == "-"].index
        gene_pairs.extend(list(product(genes_plus, genes_minus)))
    gene_pairs = np.asarray(gene_pairs)

    # Construct comparison dataframe
    #   Add each gene pair twice to enable computing correlation of Cell
    #   Ranger counts with velocyto unspliced in either direction
    overlapping_genes_df = pd.merge(
        gene_positions.loc[gene_pairs[:, 0], :].reset_index(),
        gene_positions.loc[gene_pairs[:, 1], :].reset_index(),
        left_index=True,
        right_index=True,
        suffixes=["_x", "_y"],
    )

    # Compute overlaps
    overlapping_genes_df["overlap_bases"] = (
        np.min(overlapping_genes_df[["End_x", "End_y"]], axis=1)
        - np.max(overlapping_genes_df[["Start_x", "Start_y"]], axis=1)
        + 1
    )
    overlapping_genes_df = overlapping_genes_df[overlapping_genes_df.overlap_bases > 0]
    overlapping_genes_df["overlap_relative_total"] = (
        overlapping_genes_df.overlap_bases
        / (
            np.max(overlapping_genes_df[["End_x", "End_y"]], axis=1)
            - np.min(overlapping_genes_df[["Start_x", "Start_y"]], axis=1)
            + 1
        )
    )
    overlapping_genes_df[["overlap_relative_gene_x", "overlap_relative_gene_y"]] = (
        overlapping_genes_df.overlap_bases.values[:, None]
        / (
            overlapping_genes_df[["End_x", "End_y"]].values
            - overlapping_genes_df[["Start_x", "Start_y"]].values
            + 1
        )
    )
    overlapping_genes_df["overlap_relative_average"] = (
        overlapping_genes_df.overlap_bases
        / np.mean(
            overlapping_genes_df[["End_x", "End_y"]].values
            - overlapping_genes_df[["Start_x", "Start_y"]].values
            + 1,
            axis=1,
        )
    )
    return overlapping_genes_df


def cosine(x: np.ndarray, y: np.ndarray, axis=1):
    """
    Compute the row- or column-wise cosine similarity.

    Parameters
    ----------
    x : np.ndarray
        First array.
    y : np.ndarray
        Second array.
    axis : int, optional
        Axis for computing the similarity. Use 0 for computing similarity
        per column and 1 for computing similarity per row. By default 1.

    Returns
    -------
    _type_
        _description_
    """
    return np.sum(x * y, axis=axis) / (
        np.sqrt((x**2).sum(axis=axis)) * np.sqrt((y**2).sum(axis=axis))
    )


def euclidean(x: np.ndarray, y: np.ndarray, axis=1):
    """
    Compute the row- or column-wise Euclidean distance.

    Parameters
    ----------
    x : np.ndarray
        First array.
    y : np.ndarray
        Second array.
    axis : int, optional
        Axis for computing the distance. Use 0 for computing distance
        per column and 1 for computing distance per row. By default 1.
    """
    return np.sqrt(((x - y) ** 2).sum(axis=axis))
