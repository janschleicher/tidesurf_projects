import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--orientation",
        type=str,
        default="sense",
        choices=["sense", "antisense"],
        help="Orientation of reads with respect to transcripts. For 10x"
        " Genomics, use 'sense' for three prime and 'antisense' for "
        "five prime.",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="tidesurf_out", help="Output directory."
    )
    parser.add_argument(
        "--min_intron_overlap",
        type=int,
        default=5,
        help="Minimum number of bases that a read must overlap with an "
        "intron to be considered intronic.",
    )
    parser.add_argument(
        "--multi_mapped_reads",
        action="store_true",
        help="Take reads mapping to multiple genes into account "
        "(default: reads mapping to more than one gene are discarded).",
    )
    parser.add_argument(
        "multi_sample_dir",
        metavar="MULTI_SAMPLE_DIR",
        help="Directory containing multiple Cell Ranger output directories.",
    )
    parser.add_argument(
        "gtf_file", metavar="GTF_FILE", help="GTF file with transcript information."
    )
    args = parser.parse_args()
    
    sample_dirs = [
        item for item in os.listdir(args.multi_sample_dir)
        if os.path.isdir(os.path.join(args.multi_sample_dir, item))
    ]
    
    for sample_dir in sample_dirs:
        os.system(
            f"tidesurf -o {os.path.join(args.output, sample_dir)} "
            f"--orientation {args.orientation} "
            f"--min_intron_overlap {args.min_intron_overlap} "
            f"{'--multi_mapped_reads ' if args.multi_mapped_reads else ''}"
            f"{os.path.join(args.multi_sample_dir, sample_dir)} {args.gtf_file}"
        )


if __name__ == "__main__":
    main()