import os
from typing import List, Literal

import argparse


def run_alevin_fry(
    fastq_prefixes: List[str],
    idx_path: str,
    tgm_path: str,
    out_dir: str,
    output_prefix: str = "",
    chemistry: Literal["gemcode", "chromium", "chromiumV3"] = "chromium",
    sketch: bool = True,
    whitelist: str = None,
    five_prime: bool = False,
    orient: str = None,
    valid_bc: str = None,
    threads: int = os.cpu_count(),
):
    """
    Run alevin-fry pipeline for 10x scRNA-seq data.
    :param fastq_prefixes: List of prefixes for FastQ files.
    :param idx_path: Path to salmon index.
    :param tgm_path: Path to transcript-to-gene map.
    :param out_dir: Output directory.
    :param output_prefix: Prefix for map and quant output directories. Defaults to empty string.
    :param chemistry: 10x reaction chemistry version. Defaults to "chromium".
    :param sketch: Whether to use pseudoalignment instead of selective alignment. Defaults to True.
    :param whitelist: Barcode whitelist, e.g. from CellRanger output. Defaults to None.
    :param five_prime: Switching to 5' protocol. Defaults to False.
    :param orient: Mapping orientation; defaults to fw for 3' and rc for 5'. Defaults to None.
    :param valid_bc: File with a list of valid barcodes instead of knee-distance method. Defaults to None.
    :param threads: Number of threads for multiprocessing. Defaults to number of CPUs.
    :raises FileNotFoundError: If FastQ files are not found.
    :return: The function does not return anything.
    """
    # Create list of R1 and R2 FastQs
    if os.path.isfile(f"{fastq_prefixes[0]}_R1_001.fastq.gz"):
        r1_str = " ".join([f"{x}_R1_001.fastq.gz" for x in fastq_prefixes])
        r2_str = " ".join([f"{x}_R2_001.fastq.gz" for x in fastq_prefixes])
    elif os.path.isfile(f"{fastq_prefixes[0]}_R1_001.fastq"):
        r1_str = " ".join([f"{x}_R1_001.fastq" for x in fastq_prefixes])
        r2_str = " ".join([f"{x}_R2_001.fastq" for x in fastq_prefixes])
    else:
        raise FileNotFoundError("FastQ files not found")

    # Specify library, orientation etc
    library = (
        "ISF" if five_prime else "ISR"
    )  # inward stranded forward/reverse (read 1 from f/r strand)
    if orient is None:
        orient = "rc" if five_prime else "fw"
    dist_method = "--knee-distance" if valid_bc is None else f"--valid-bc {valid_bc}"
    sketch_str = "--sketch" if sketch else f"--noQuant --tgMap {tgm_path}"

    # Run salmon alevin
    wl = "" if whitelist is None else f"--whitelist {whitelist} "
    os.system(
        f"salmon alevin -i {idx_path} -p {threads} -l {library} --{chemistry} {sketch_str} {wl}"
        f"-1 {r1_str} -2 {r2_str} -o {os.path.join(out_dir, f'{output_prefix}map')}"
    )

    os.system(
        f"alevin-fry generate-permit-list -d {orient} {dist_method}"
        f" -i {os.path.join(out_dir, f'{output_prefix}map')}"
        f" -o {os.path.join(out_dir, f'{output_prefix}quant')}"
    )
    os.system(
        f"alevin-fry collate -t {threads} -i {os.path.join(out_dir, f'{output_prefix}quant')}"
        f" -r {os.path.join(out_dir, f'{output_prefix}map')}"
    )
    os.system(
        f"alevin-fry quant -t {threads} -i {os.path.join(out_dir, f'{output_prefix}quant')}"
        f" -o {os.path.join(out_dir, f'{output_prefix}quant_res')} --tg-map {tgm_path}"
        f" --resolution cr-like --use-mtx"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fastq-prefix",
        required=True,
        nargs="+",
        type=str,
        help="Path to FastQ files preceding '_R[12]_001.fastq.gz'",
    )
    parser.add_argument(
        "--idx-path", required=True, type=str, help="Path to salmon index"
    )
    parser.add_argument(
        "--tgm-path", required=True, type=str, help="Path to transcript-to-gene map"
    )
    parser.add_argument(
        "--chemistry",
        default="chromium",
        type=str,
        choices=["gemcode", "chromium", "chromiumV3"],
        help="10x reaction chemistry version",
    )
    parser.add_argument(
        "--output-dir", default="alevin_fry", type=str, help="Output directory"
    )
    parser.add_argument(
        "--selective-alignment",
        action="store_true",
        help="Use selective alignment instead of pseudoalignment",
    )
    parser.add_argument(
        "--whitelist",
        default=None,
        type=str,
        help="Barcode whitelist, e.g. from CellRanger output",
    )
    parser.add_argument(
        "--five-prime", action="store_true", help="Switching to 5' protocol"
    )
    parser.add_argument(
        "--orient",
        default=None,
        choices=["fw", "rc", "both"],
        help="Mapping orientation; defaults to fw for 3' and rc for 5'",
    )
    parser.add_argument(
        "--valid-bc",
        default=None,
        type=str,
        help="File with a list of valid barcodes instead of knee-distance method",
    )
    parser.add_argument(
        "--output-prefix",
        default="",
        type=str,
        help="Prefix for map and quant output directories",
    )
    parser.add_argument(
        "--threads",
        default=os.cpu_count(),
        type=int,
        help="Number of threads for multiprocessing",
    )

    args = parser.parse_args()

    run_alevin_fry(
        fastq_prefixes=args.fastq_prefix,
        idx_path=args.idx_path,
        tgm_path=args.tgm_path,
        out_dir=args.output_dir,
        output_prefix=args.output_prefix,
        chemistry=args.chemistry,
        sketch=not args.selective_alignment,
        whitelist=args.whitelist,
        five_prime=args.five_prime,
        orient=args.orient,
        valid_bc=args.valid_bc,
        threads=args.threads,
    )
