import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pysam
from tqdm.auto import tqdm


def reverse_complement(sequence: str):
    complement_dict = defaultdict(lambda: "N", {"A": "T", "C": "G", "G": "C", "T": "A"})
    rc = ""
    for base in sequence[::-1]:
        rc += complement_dict[base]
    return rc


def main():
    parser = argparse.ArgumentParser("umi_matching")
    parser.add_argument(
        "-cr",
        "--cellranger-bam",
        type=Path,
        help="Cell Ranger BAM file for UMI matching.",
    )
    parser.add_argument(
        "-pb", "--pacbio-bam", help="PacBio IsoSeq mapped BAM file for UMI matching."
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        type=Path,
        help="Output directy for saving barcodes and UMIs.",
    )
    args = parser.parse_args()

    bam_paths = {"cellranger": args.cellranger_bam, "pacbio": args.pacbio_bam}

    # Extract cell barcodes and UMIs from the BAM files
    cell_barcodes = {}
    umis = {}
    tags = {}
    umi_tags = {"cellranger": "UB", "pacbio": "XM"}
    for key, path in tqdm(bam_paths.items(), desc="Processing BAM files"):
        # Open BAM file
        bam_file = pysam.AlignmentFile(path, mode="r", check_sq=False)

        # Prepare sets for saving barcodes and UMIs
        cell_barcodes[key] = set()
        umis[key] = set()
        tags[key] = set()

        # Process reads
        for bam_read in tqdm(
            bam_file.fetch(until_eof=True), desc="Processing reads", leave=False
        ):
            if bam_read.has_tag("CB") and bam_read.has_tag(umi_tags[key]):
                cbc = str(bam_read.get_tag("CB"))
                umi = str(bam_read.get_tag(umi_tags[key]))
                if key == "pacbio":
                    cbc = reverse_complement(cbc)
                    umi = reverse_complement(umi)
                elif key == "cellranger":
                    cbc = cbc.strip("-1")
                tag = f"{cbc}_{umi}"
                cell_barcodes[key].add(cbc)
                umis[key].add(umi)
                tags[key].add(tag)

    # Save the barcode and UMI lists
    out_dir = args.out_dir

    for key, barcodes in cell_barcodes.items():
        pd.Series(list(barcodes)).to_csv(
            out_dir / f"{key}_barcodes.csv", index=False, header=False
        )

    for key, umi_set in umis.items():
        pd.Series(list(umi_set)).to_csv(
            out_dir / f"{key}_umis.csv", index=False, header=False
        )

    for key, tag_set in tags.items():
        pd.Series(list(tag_set)).to_csv(
            out_dir / f"{key}_tags.csv", index=False, header=False
        )

    common_tags = tags["cellranger"].intersection(tags["pacbio"])
    pd.Series(list(common_tags)).to_csv(
        out_dir / "common_tags.csv", index=False, header=False
    )


if __name__ == "__main__":
    main()
