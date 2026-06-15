import argparse
from pathlib import Path

import pandas as pd
import pysam
from tqdm.auto import tqdm
from umi_matching import reverse_complement


def main():
    # Set up command line arguments
    parser = argparse.ArgumentParser("filter_bam")
    parser.add_argument(
        "-f",
        "--filter-list",
        type=Path,
        help="Path to text file with one tag (CBC_UMI) per line for read filtering, "
        "without header and index.",
    )
    parser.add_argument("-i", "--input", type=Path, help="Input BAM file path.")
    parser.add_argument(
        "-o", "--output", type=Path, help="Filtered output BAM file path."
    )
    parser.add_argument(
        "-pb",
        "--pacbio",
        action="store_true",
        help="Whether the BAM file is from PacBio. If set, the reverse complement"
        "of CBC and UMI will be used, and the UMI tag will be changed.",
    )
    args = parser.parse_args()

    # Set up input and output
    filter_list = pd.read_csv(args.filter_list, header=None, index_col=None)
    filter_list = set(filter_list[0])
    input_bam = pysam.AlignmentFile(args.input, mode="rb")
    output_bam = pysam.AlignmentFile(args.output, mode="wb", template=input_bam)

    # Filter reads
    umi_tag = "XM" if args.pacbio else "UB"
    for read in tqdm(input_bam.fetch(until_eof=True), desc="Processing reads"):
        if read.has_tag("CB") and read.has_tag(umi_tag):
            cbc = str(read.get_tag("CB"))
            umi = str(read.get_tag(umi_tag))
        else:
            continue

        if args.pacbio:
            cbc = reverse_complement(cbc)
            umi = reverse_complement(umi)
        else:
            cbc = cbc.strip("-1")
        tag = f"{cbc}_{umi}"
        if tag in filter_list:
            output_bam.write(read)


if __name__ == "__main__":
    main()
