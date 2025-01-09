import os
import argparse

from pyroe import load_fry


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-q",
        "--quant-res-path",
        required=True,
        type=str,
        help="Path to alevin-fry quantification results",
    )
    parser.add_argument(
        "-e2g",
        "--ens-to-gene",
        required=True,
        type=str,
        help="File with Ensemble ID to gene name mappings",
    )
    parser.add_argument(
        "-f", "--filename", type=str, default="adata", help="Name of output file"
    )
    parser.add_argument(
        "-of",
        "--output-format",
        type=str,
        default="velocity",
        choices=["scRNA", "all", "raw", "raw_all", "S+A", "velocity", "velocity_all"],
        help="Output format of AnnData count matrices",
    )
    args = parser.parse_args()

    # Load into AnnData object and use gene names as var names
    if args.output_format == "velocity_all":
        output_format = "velocity"
    elif args.output_format == "raw_all":
        output_format = "raw"
    else:
        output_format = args.output_format
    adata = load_fry(args.quant_res_path, output_format=output_format)
    if args.output_format == "velocity_all":
        adata.X = adata.layers["spliced"] + adata.layers["unspliced"]
    elif args.output_format == "raw_all":
        adata.X = (
            adata.layers["spliced"]
            + adata.layers["unspliced"]
            + adata.layers["ambiguous"]
        )
    e2g = dict([line.rstrip().split() for line in open(args.ens_to_gene).readlines()])
    adata.var["gene_id"] = adata.var_names.copy()
    adata.var_names = [e2g[e] for e in adata.var_names]
    adata.var["gene_name"] = adata.var_names.copy()
    adata.var_names_make_unique()

    adata.write_h5ad(
        os.path.join(
            args.quant_res_path,
            f"{args.filename}{'' if args.filename.endswith('.h5ad') else '.h5ad'}",
        )
    )


if __name__ == "__main__":
    main()
