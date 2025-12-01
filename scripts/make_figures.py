import os
import string
from typing import Iterable, Optional

import matplotlib as mpl
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import seaborn as sns
from matplotlib import pyplot as plt

# Set up plotting parameters
sc.set_figure_params(vector_friendly=True, dpi_save=400, scanpy=False, fontsize=10)
scv.set_figure_params(vector_friendly=True, dpi_save=400, fontsize=10)
plt.style.use("matplotlibrc.txt")
AXLAB_KWS = {"fontsize": 12, "fontweight": "bold", "va": "top", "ha": "left"}
BOXPLOT_PROPS = {
    prop: {"edgecolor" if prop == "boxprops" else "color": "k"}
    for prop in ["boxprops", "medianprops", "whiskerprops", "capprops"]
}
FIG_WIDTH = 19.05 / 2.54  # in inches
FIG_WIDTH_SINGLE = 13.2 / 2.54  # in inches
FIG_HEIGHT_MAX = 22.23 / 2.54
FIG_DIR = "../figures"
PIL_KWARGS = {"compression": "tiff_lzw"}  # needed if saving as TIFF
# Paul Tol's muted palette (https://sronpersonalpages.nl/~pault/)
PALETTE = [
    "#CC6677",
    "#332288",
    "#DDCC77",
    "#117733",
    "#88CCEE",
    "#882255",
    "#44AA99",
    "#999933",
    "#AA4499",
    "#DDDDDD",
]
# Paul Tol's vibrant palette (https://sronpersonalpages.nl/~pault/)
PALETTE_COMP = [
    "#EE7733",
    "#0077BB",
    "#33BBEE",
    "#EE3377",
    "#CC3311",
    "#009988",
    "#BBBBBB",
]

DATASETS_ORDER = ["pancreas", "retina", "stewart", "fu"]
DATASETS_LABELS = {
    "pancreas": "Pancreas",
    "retina": "Retina",
    "stewart": "B cells",
    "fu": "T cells",
    "fu_pe": "T cells (PE)",
    "fu_se": "T cells (SE)",
}
METHODS_ORDER = [
    "velocyto",
    "alevin-fry",
    "starsolo",
    "tidesurf",
    "cellranger",
    "alevin-fry_starsolo",
    "tidesurf_alevin-fry",
    "tidesurf_starsolo",
    "tidesurf_velocyto",
    "velocyto_alevin-fry",
    "velocyto_starsolo",
]
METHODS_DICT = {
    "velocyto": "velocyto",
    "alevin-fry": "alevin-fry",
    "starsolo": "STARsolo",
    "tidesurf": "tidesurf",
    "cellranger": "Cell Ranger",
}
SEQ_TYPES = ["3'", "5'"]
SPLICE_STATES = ["spliced", "unspliced", "ambiguous"]


def read_dataframes(name: str, datasets: Iterable = DATASETS_ORDER):
    df = (
        pd.concat(
            {
                dataset: pd.read_parquet(f"../out/{dataset}/{name}.parquet.gz")
                for dataset in datasets
            },
            names=["dataset", "idx"],
        )
        .reset_index()
        .drop(columns="idx")
    )
    return df


def split_boxplot(
    sub_fig: mpl.figure.SubFigure,
    data: pd.DataFrame,
    y: str,
    y_label: Optional[str] = None,
    hue: str = "method",
    palette: str = "Dark2",
):
    axs = sub_fig.subplots(1, 2, sharey=True)
    methods = data[hue].unique()
    for i, (seq_type, d_order) in enumerate(
        zip(SEQ_TYPES, [DATASETS_ORDER[:2], DATASETS_ORDER[2:]])
    ):
        sns.boxplot(
            data[data["dataset"].isin(d_order)],
            x="dataset",
            y=y,
            hue=hue,
            order=d_order,
            hue_order=[x for x in METHODS_ORDER if x in methods],
            palette=palette,
            showfliers=False,
            ax=axs[i],
            **BOXPLOT_PROPS,
        )
        axs[i].text(
            0.5, 1, seq_type, ha="center", va="bottom", transform=axs[i].transAxes
        )
        axs[i].legend().remove()
        axs[i].set_xlabel("")
        axs[i].set_xticks(
            ticks=axs[i].get_xticks(),
            labels=[
                DATASETS_LABELS[lab.get_text()] for lab in axs[i].get_xticklabels()
            ],
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )
    sns.despine(ax=axs[-1], left=True)
    axs[-1].tick_params(axis="y", which="both", left=False)
    if y_label is not None:
        axs[0].set_ylabel(y_label)


def figure_1():
    methods_order = ["velocyto", "cellranger"]
    datasets_order = ["stewart", "fu"]
    left_label_x = 0.007

    # Load all DataFrames
    total_counts = read_dataframes("total_counts", datasets=datasets_order)
    counts_diff = read_dataframes("counts_diff", datasets=datasets_order)
    counts_corr = read_dataframes("counts_corr_cellranger", datasets=datasets_order)
    counts_corr = counts_corr[counts_corr.alignment_type != "both"]
    spliced_unspliced_ratios = read_dataframes(
        "spliced_unspliced_ratios", datasets=datasets_order
    )
    for df in [total_counts, counts_diff, counts_corr, spliced_unspliced_ratios]:
        df.loc[df.dataset == "fu", "dataset"] = (
            df.loc[df.dataset == "fu", "dataset"].astype(str)
            + "_"
            + df.loc[df.dataset == "fu", "alignment_type"].str.lower()
        )
    overlapping_genes = read_dataframes("gene_overlaps", datasets=datasets_order)

    adata = sc.read_h5ad("../data/stewart/adata/adata_raw_velocyto.h5ad")

    # Subset to genes with mean expression >= 1
    overlapping_genes = overlapping_genes[
        np.asarray(adata[:, overlapping_genes.Gene_x].X.mean(axis=0) >= 1).ravel()
        & (overlapping_genes.alignment_type != "both")
    ]

    overlapping_genes_scatter = (
        overlapping_genes[overlapping_genes.dataset == "stewart"]
        .sort_values(by=["overlap_relative_gene_x", "overlap_bases"], ascending=False)
        .reset_index(drop=True)
    )

    datasets_order = ["stewart", "fu_se", "fu_pe"]

    fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT_MAX), layout="constrained")
    sub_figs = fig.subfigures(3, 1, height_ratios=[7, 7, 12])

    # Top plots
    sub_figs_top = sub_figs[0].subfigures(1, 4)
    axs_top = [sub_fig.subplots(1, 1) for sub_fig in sub_figs_top]
    sns.boxplot(
        total_counts[total_counts.method.isin(methods_order)],
        x="dataset",
        y="counts",
        hue="method",
        hue_order=methods_order,
        order=datasets_order,
        palette=[PALETTE[0], PALETTE[4]],
        showfliers=False,
        log_scale=True,
        ax=axs_top[0],
        **BOXPLOT_PROPS,
    )
    axs_top[0].set_ylabel("counts per cell")

    sns.boxplot(
        counts_diff[counts_diff.method.isin(methods_order)],
        x="dataset",
        y="difference",
        palette=[PALETTE[0]],
        order=datasets_order,
        showfliers=False,
        ax=axs_top[1],
        **BOXPLOT_PROPS,
    )
    axs_top[1].axhline(0, c="grey", linestyle=":")

    sns.boxplot(
        counts_corr[
            counts_corr.method.isin(methods_order)
            & counts_corr.dataset.isin(datasets_order)
        ],
        x="dataset",
        y="pearsonr",
        palette=[PALETTE[0]],
        order=datasets_order,
        showfliers=False,
        ax=axs_top[2],
        **BOXPLOT_PROPS,
    )
    axs_top[2].set_ylabel("Pearson corr.")

    spliced_unspliced_ratios.ratio = spliced_unspliced_ratios.ratio * 100
    sns.boxplot(
        spliced_unspliced_ratios[
            spliced_unspliced_ratios.method.isin(methods_order)
            & spliced_unspliced_ratios.dataset.isin(datasets_order)
            & (spliced_unspliced_ratios.genes == "all")
            & (spliced_unspliced_ratios.splice_state == "unspliced")
        ],
        x="dataset",
        y="ratio",
        palette=[PALETTE[0]],
        order=datasets_order,
        showfliers=False,
        ax=axs_top[3],
        **BOXPLOT_PROPS,
    )
    axs_top[3].set_ylabel("% unspliced")

    for ax in axs_top:
        ax.set_xlabel("")
        ax.set_xticks(
            ticks=ax.get_xticks(),
            labels=[DATASETS_LABELS[lab.get_text()] for lab in ax.get_xticklabels()],
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )
        if ax.get_legend():
            ax.get_legend().remove()

    axs_top[0].text(
        left_label_x,
        0.98,
        "A",
        transform=mpl.transforms.blended_transform_factory(
            sub_figs[0].transSubfigure, sub_figs_top[0].transSubfigure
        ),
        **AXLAB_KWS,
    )
    for i, label in enumerate(string.ascii_uppercase[1 : len(sub_figs_top)]):
        axs_top[i + 1].text(
            0.02,
            0.98,
            label,
            transform=sub_figs_top[i + 1].transSubfigure,
            **AXLAB_KWS,
        )

    handles, labels = axs_top[0].get_legend_handles_labels()
    sub_figs[0].legend(
        handles,
        [METHODS_DICT[lab] for lab in labels],
        loc="outside lower center",
        ncols=2,
        frameon=False,
    )

    # Center plots
    axs_center = sub_figs[1].subplots(1, 2)

    # Positions of overlapping transcripts RPL35 (+) and IQCG (-)
    transcript_rpl35a = [197_950_190, 197_956_610]
    exons_rpl35a = [
        [197_950_190, 197_950_221],
        [197_950_936, 197_950_936],
        [197_951_159, 197_951_311],
        [197_954_003, 197_954_147],
        [197_955_750, 197_956_610],
    ]
    cds_rpl35a = [
        [197_950_968, 197_950_978],
        [197_951_159, 197_951_311],
        [197_954_003, 197_954_147],
        [197_955_750, 197_955_770],
    ]

    transcript_iqcg = [197_943_778, 197_959_991]  # real start: 197_889_077
    exons_iqcg = [
        [197_959_845, 197_959_991],
        [197_959_529, 197_959_719],
        [197_945_620, 197_945_686],
        [197_943_778, 197_944_051],
        # [197_938_548, 197_938_780],
        # [197_932_174, 197_932_302],
        # [197_926_024, 197_926_105],
        # [197_913_881, 197_914_042],
        # [197_912_675, 197_912_749],
        # [197_892_628, 197_892_759],
        # [197_891_446, 197_891_547],
        # [197_889_077, 197_889_714],
    ]
    cds_iqcg = [
        [197_945_620, 197_945_627],
        [197_943_778, 197_944_051],
        # [197_938_548, 197_938_780],
        # [197_932_174, 197_932_302],
        # [197_926_024, 197_926_105],
        # [197_913_881, 197_914_042],
        # [197_912_675, 197_912_749],
        # [197_892_628, 197_892_759],
        # [197_891_446, 197_891_547],
        # [197_889_583, 197_889_714],
    ]

    y_pos = 0.2
    x_quiver = np.linspace(transcript_rpl35a[0], transcript_rpl35a[1], 10)
    y_quiver = np.ones_like(x_quiver) * y_pos
    axs_center[0].quiver(
        x_quiver[:-1],
        y_quiver[:-1],
        x_quiver[1:] - x_quiver[:-1],
        y_quiver[1:] - y_quiver[:-1],
        headaxislength=3,
        headlength=3,
        color="grey",
    )
    axs_center[0].plot(
        np.asarray(exons_rpl35a).T,
        [y_pos, y_pos],
        lw=10,
        c=mpl.colormaps["tab20"].colors[1],
    )
    axs_center[0].plot(
        np.asarray(cds_rpl35a).T,
        [y_pos, y_pos],
        lw=10,
        c=mpl.colormaps["tab20"].colors[0],
    )
    axs_center[0].text(
        (transcript_rpl35a[0] + transcript_rpl35a[1]) / 2,
        y_pos + 0.2,
        "RPL35A",
        va="bottom",
        ha="center",
        size=8,
    )
    x_quiver = np.linspace(transcript_iqcg[1], transcript_iqcg[0], 20)
    y_quiver = np.ones_like(x_quiver) * (-y_pos)
    axs_center[0].quiver(
        x_quiver[:-1],
        y_quiver[:-1],
        x_quiver[1:] - x_quiver[:-1],
        y_quiver[1:] - y_quiver[:-1],
        headaxislength=3,
        headlength=3,
        color="grey",
    )
    axs_center[0].plot(
        [transcript_iqcg[0] - 1_000, transcript_iqcg[0]],
        [-y_pos, -y_pos],
        c="grey",
        ls=":",
    )
    axs_center[0].plot(
        np.asarray(exons_iqcg).T,
        [-y_pos, -y_pos],
        lw=10,
        c=mpl.colormaps["tab20"].colors[3],
    )
    axs_center[0].plot(
        np.asarray(cds_iqcg).T,
        [-y_pos, -y_pos],
        lw=10,
        c=mpl.colormaps["tab20"].colors[2],
    )
    axs_center[0].text(
        (transcript_iqcg[0] + transcript_iqcg[1]) / 2,
        -y_pos - 0.2,
        "IQCG",
        va="top",
        ha="center",
        size=8,
    )
    # Draw forward and reverse strand
    xlim = np.asarray(axs_center[0].get_xlim())
    arrow_len = xlim[1] - xlim[0]
    axs_center[0].quiver(
        xlim,
        np.asarray([0.8, -0.8]),
        np.asarray([arrow_len, -arrow_len]),
        np.zeros(2),
        angles="xy",
        scale_units="xy",
        scale=1,
        headwidth=4,
        headlength=4,
        headaxislength=4,
        color="k",
    )
    axs_center[0].text(
        (xlim[0] + xlim[1]) / 2, 0.9, "forward strand", va="bottom", ha="center", size=8
    )
    axs_center[0].text(
        (xlim[0] + xlim[1]) / 2, -0.9, "reverse strand", va="top", ha="center", size=8
    )
    axs_center[0].legend(
        [
            (
                mpl.patches.Rectangle(
                    (0, 0),
                    2,
                    1,
                    facecolor=mpl.colormaps["tab20"].colors[i],
                    edgecolor=None,
                ),
                mpl.patches.Rectangle(
                    (0, 0),
                    2,
                    1,
                    facecolor=mpl.colormaps["tab20"].colors[i + 2],
                    edgecolor=None,
                ),
            )
            for i in [0, 1]
        ],
        ["CDS", "UTR"],
        handler_map={tuple: mpl.legend_handler.HandlerTuple(ndivide=None, pad=0.3)},
        loc="upper center",
        bbox_to_anchor=(0.5, -0.1),
        frameon=False,
        ncols=2,
        handlelength=2,
        handleheight=1.0,
        handletextpad=0.6,
    )
    axs_center[0].axis("off")

    overlapping_genes.loc[overlapping_genes.dataset == "fu", "dataset"] = (
        overlapping_genes.loc[overlapping_genes.dataset == "fu", "dataset"].astype(str)
        + "_"
        + overlapping_genes.loc[
            overlapping_genes.dataset == "fu", "alignment_type"
        ].str.lower()
    )
    overlapping_genes["overlap"] = pd.cut(
        overlapping_genes.overlap_relative_gene_x,
        [0, 0.5, 0.75, 1],
        right=True,
        include_lowest=True,
        labels=[
            r"overlap$\leq50\%$",
            r"$50\%<$overlap$\leq75\%$",
            r"overlap$>75\%$",
        ],
    )
    sns.boxplot(
        overlapping_genes,
        x="dataset",
        y="pearsonr",
        hue="overlap",
        palette="Dark2",
        order=datasets_order,
        showfliers=False,
        ax=axs_center[1],
        **BOXPLOT_PROPS,
    )
    axs_center[1].legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
    axs_center[1].set_xlabel("")
    axs_center[1].set_ylabel("Pearson corr.")
    axs_center[1].set_xticks(
        ticks=axs_center[1].get_xticks(),
        labels=[
            DATASETS_LABELS[lab.get_text()] for lab in axs_center[1].get_xticklabels()
        ],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
    )
    sub_figs[1].text(
        left_label_x, 1, "E", transform=sub_figs[1].transSubfigure, **AXLAB_KWS
    )
    axs_center[1].text(
        -0.2,
        1,
        "F",
        transform=mpl.transforms.blended_transform_factory(
            axs_center[1].transAxes, sub_figs[1].transSubfigure
        ),
        **AXLAB_KWS,
    )

    # Bottom plots
    axs_bottom = sub_figs[2].subplots(2, 4)
    for i, row in overlapping_genes_scatter.head(8).iterrows():
        gene_x = row["Gene_x"]
        gene_y = row["Gene_y"]
        rho = row["pearsonr"]
        ax = axs_bottom.flat[i]
        ax.scatter(
            adata[:, gene_x].X.toarray().ravel(),
            adata[:, gene_y].layers["unspliced"].toarray().ravel(),
            s=0.5,
            c="k",
            rasterized=True,
        )
        ax.text(
            0.05,
            0.95,
            f"$\\rho={rho:.2f}$",
            size=8,
            ha="left",
            va="top",
            transform=ax.transAxes,
        )
        ax.set_xlabel(gene_x)
        ax.set_ylabel(gene_y)
    sub_figs[2].text(
        left_label_x, 1.03, "G", transform=sub_figs[2].transSubfigure, **AXLAB_KWS
    )
    sub_figs[2].supxlabel("Cell Ranger counts")
    sub_figs[2].supylabel("velocyto unspliced counts")

    fig.savefig(os.path.join(FIG_DIR, "fig1"))
    plt.close()


def supplementary_figure_1():
    # Load all DataFrames
    spliced_unspliced_diff = read_dataframes("spliced_unspliced_diff")
    spliced_unspliced_corr = read_dataframes("spliced_unspliced_corr")

    # Make figure
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_WIDTH / 1.5), layout="constrained")
    sub_figs = fig.subfigures(2, 3).ravel()

    # Difference in spliced/unspliced/ambiguous counts per cell
    for sub_fig, splice_state in zip(sub_figs[:3], SPLICE_STATES):
        split_boxplot(
            sub_fig,
            spliced_unspliced_diff[
                (spliced_unspliced_diff["genes"] == "all")
                & (spliced_unspliced_diff["splice_state"] == splice_state)
            ],
            "diff",
            y_label="difference",
        )
        sub_fig.suptitle(splice_state)
        if splice_state != SPLICE_STATES[0]:
            sub_fig.get_axes()[0].set_ylabel("")
        for ax in sub_fig.get_axes():
            ax.axhline(0, color="grey", linestyle="--", zorder=0)

    # Pearson correlation between methods
    for sub_fig, splice_state in zip(sub_figs[3:], SPLICE_STATES):
        split_boxplot(
            sub_fig,
            spliced_unspliced_corr[
                spliced_unspliced_corr["splice_state"] == splice_state
            ],
            "pearsonr",
            y_label="Pearson corr.",
            hue="comparison",
            palette="Accent",
        )
        if splice_state != SPLICE_STATES[0]:
            sub_fig.get_axes()[0].set_ylabel("")

    handles_top, labels_top = sub_figs[0].get_axes()[0].get_legend_handles_labels()
    handles_bottom, labels_bottom = (
        sub_figs[3].get_axes()[0].get_legend_handles_labels()
    )
    fig.legend(
        handles_top + handles_bottom,
        [lab.replace("_", " vs.\n") for lab in labels_top + labels_bottom],
        ncols=5,
        loc="outside lower center",
    )
    for i, label in enumerate(string.ascii_uppercase[:6]):
        sub_figs[i].text(
            0.02, 0.98, label, transform=sub_figs[i].transSubfigure, **AXLAB_KWS
        )
    fig.savefig(os.path.join(FIG_DIR, "suppfig1.pdf"))
    plt.close()


def supplementary_figure_2():
    # Make figure
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_WIDTH), layout="constrained")
    sub_figs = fig.subfigures(3, 1)

    # Load and plot AnnData objects for retina dataset
    adatas_retina = {
        method: sc.read_h5ad(f"../out/retina/adata_{method}.h5ad")
        for method in METHODS_ORDER[:3]
    }
    axs_top = sub_figs[0].subplots(1, 3)
    for i, (method, adata) in enumerate(adatas_retina.items()):
        scv.pl.velocity_embedding_grid(
            adata,
            basis="umap",
            color="Annotation",
            size=20,
            alpha=0.5,
            arrow_length=2,
            arrow_size=(10, 18, 8),
            arrow_color="black",
            density=0.8,
            show=False,
            ax=axs_top[i],
            title=method,
            legend_loc="right margin" if i == len(adatas_retina) - 1 else "none",
        )

    # Load and plot AnnData objects for Stewart dataset
    adatas_stewart = {
        method: sc.read_h5ad(f"../out/stewart/adata_{method}.h5ad")
        for method in METHODS_ORDER[:3]
    }
    axs_mid = sub_figs[1].subplots(1, 3)
    for i, (method, adata) in enumerate(adatas_stewart.items()):
        scv.pl.velocity_embedding_grid(
            adata,
            basis="umap",
            color="cell_type",
            size=20,
            alpha=0.5,
            arrow_length=2,
            arrow_size=(10, 18, 8),
            arrow_color="black",
            density=0.8,
            show=False,
            ax=axs_mid[i],
            title=method,
            legend_loc="right margin" if i == len(adatas_stewart) - 1 else "none",
        )

    # Load and plot AnnData objects for Fu dataset
    adatas_fu = {
        method: sc.read_h5ad(f"../out/fu/adata_{method}.h5ad")
        for method in METHODS_ORDER[:3]
    }
    axs_bottom = sub_figs[2].subplots(1, 3)
    for i, (method, adata) in enumerate(adatas_fu.items()):
        scv.pl.velocity_embedding_grid(
            adata,
            basis="umap",
            color="groups",
            size=20,
            alpha=0.5,
            arrow_length=2,
            arrow_size=(10, 18, 8),
            arrow_color="black",
            density=0.8,
            show=False,
            ax=axs_bottom[i],
            title=method,
            legend_loc="right margin" if i == len(adatas_fu) - 1 else "none",
        )

    for i, label in enumerate(string.ascii_uppercase[:3]):
        sub_figs[i].text(
            0.01, 0.98, label, transform=sub_figs[i].transSubfigure, **AXLAB_KWS
        )
    fig.savefig(os.path.join(FIG_DIR, "suppfig2.pdf"))
    plt.close()


def figure_3():
    # Load all DataFrames
    total_counts = read_dataframes("total_counts")
    counts_diff = read_dataframes("counts_diff")
    counts_corr = read_dataframes("counts_corr_cellranger")

    # Make figure
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_WIDTH / 3), layout="constrained")
    sub_figs = fig.subfigures(1, 3)
    split_boxplot(
        sub_figs[0],
        total_counts[total_counts["genes"] == "all"],
        "counts",
        palette=PALETTE,
    )
    split_boxplot(
        sub_figs[1],
        counts_diff[counts_diff["genes"] == "all"],
        "difference",
        palette=PALETTE,
    )
    for ax in sub_figs[1].get_axes():
        ax.axhline(0, color="grey", linestyle="--", zorder=0)
    split_boxplot(
        sub_figs[2],
        counts_corr[~counts_corr["alignment_type"].isin(["PE", "SE"])],
        "pearsonr",
        y_label="Pearson corr.",
        palette=PALETTE,
    )

    handles, labels = sub_figs[0].get_axes()[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        ncols=5,
        loc="outside lower center",
    )
    for i, label in enumerate(string.ascii_uppercase[: len(sub_figs)]):
        sub_figs[i].text(
            0.02, 0.98, label, transform=sub_figs[i].transSubfigure, **AXLAB_KWS
        )
    fig.savefig(os.path.join(FIG_DIR, "fig3"))
    plt.close()


def figure_4():
    # Load all DataFrames
    spliced_unspliced_counts = read_dataframes("spliced_unspliced_counts")
    spliced_unspliced_ratios = read_dataframes("spliced_unspliced_ratios")
    spliced_unspliced_cosine = read_dataframes("spliced_unspliced_cosine")

    # Make figure
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT_MAX), layout="constrained")
    sub_figs = fig.subfigures(3, 1)

    # Total spliced/unspliced/ambiguous counts per cell
    sub_figs_top = sub_figs[0].subfigures(1, 3)
    for sub_fig, splice_state in zip(sub_figs_top, SPLICE_STATES):
        split_boxplot(
            sub_fig,
            spliced_unspliced_counts[
                (spliced_unspliced_counts["genes"] == "all")
                & (spliced_unspliced_counts["splice_state"] == splice_state)
            ],
            "counts",
            palette=PALETTE,
        )
        sub_fig.suptitle(splice_state, size=10)
        if splice_state != SPLICE_STATES[0]:
            sub_fig.get_axes()[0].set_ylabel("")
    handles, labels = sub_figs_top[0].get_axes()[0].get_legend_handles_labels()
    sub_figs[0].legend(
        handles,
        [METHODS_DICT[lab] for lab in labels],
        handlelength=1.5,
        ncols=4,
        loc="outside lower center",
    )

    # Total spliced/unspliced/ambiguous counts ratios per cell
    sub_figs_center = sub_figs[1].subfigures(1, 3)
    for sub_fig, splice_state in zip(sub_figs_center, SPLICE_STATES):
        split_boxplot(
            sub_fig,
            spliced_unspliced_ratios[
                (spliced_unspliced_ratios["genes"] == "all")
                & (spliced_unspliced_ratios["splice_state"] == splice_state)
            ],
            "ratio",
            palette=PALETTE,
        )
        sub_fig.suptitle(splice_state, size=10)
        if splice_state != SPLICE_STATES[0]:
            sub_fig.get_axes()[0].set_ylabel("")
    handles, labels = sub_figs_center[0].get_axes()[0].get_legend_handles_labels()
    sub_figs[1].legend(
        handles,
        [METHODS_DICT[lab] for lab in labels],
        handlelength=1.5,
        ncols=4,
        loc="outside lower center",
    )

    # Cosine similarity between methods
    sub_figs_bottom = sub_figs[2].subfigures(1, 3)
    for sub_fig, splice_state in zip(sub_figs_bottom, SPLICE_STATES):
        split_boxplot(
            sub_fig,
            spliced_unspliced_cosine[
                spliced_unspliced_cosine["splice_state"] == splice_state
            ],
            "cosine_similarity",
            y_label="cosine similarity",
            hue="comparison",
            palette=PALETTE_COMP,
        )
        sub_fig.suptitle(splice_state, size=10)
        if splice_state != SPLICE_STATES[0]:
            sub_fig.get_axes()[0].set_ylabel("")
    handles, labels = sub_figs_bottom[0].get_axes()[0].get_legend_handles_labels()
    sub_figs[2].legend(
        handles,
        [
            " vs.\n".join(
                [METHODS_DICT[lab.split("_")[0]], METHODS_DICT[lab.split("_")[1]]]
            )
            for lab in labels
        ],
        handlelength=1.5,
        ncols=6,
        loc="outside lower center",
    )
    for i, sub_fig in enumerate(
        sub_figs_top.tolist() + sub_figs_center.tolist() + sub_figs_bottom.tolist()
    ):
        sub_fig.text(
            0.02,
            0.98,
            string.ascii_uppercase[i],
            transform=sub_fig.transSubfigure,
            **AXLAB_KWS,
        )
    fig.savefig(os.path.join(FIG_DIR, "fig4"))
    plt.close()


def figure_5():
    # Load all DataFrames
    velocities_cosine = read_dataframes("velocities_velo_cosine")
    comp_mapping = {
        "velocyto_tidesurf": "tidesurf_velocyto",
        "alevin-fry_tidesurf": "tidesurf_alevin-fry",
        "starsolo_tidesurf": "tidesurf_starsolo",
        "velocyto_alevin-fry": "velocyto_alevin-fry",
        "velocyto_starsolo": "velocyto_starsolo",
        "alevin-fry_starsolo": "alevin-fry_starsolo",
    }
    velocities_cosine["comparison"] = velocities_cosine["comparison"].map(comp_mapping)
    velocities_corr = read_dataframes("velocities_velo_corr")
    velocities_corr["comparison"] = velocities_corr["comparison"].map(comp_mapping)

    # Load AnnData objects
    adatas_pancreas = {
        method: sc.read_h5ad(f"../out/pancreas/adata_{method}.h5ad")
        for method in METHODS_ORDER[:4]
    }
    adatas_stewart = {
        method: sc.read_h5ad(f"../out/stewart/adata_{method}.h5ad")
        for method in METHODS_ORDER[:4]
    }

    # Make figure
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT_MAX), layout="constrained")
    sub_figs = fig.subfigures(1, 2, width_ratios=[2, 1])

    axs_top = sub_figs[0].subplots(4, 2)
    for i, (method, adata) in enumerate(adatas_pancreas.items()):
        scv.pl.velocity_embedding_grid(
            adata,
            basis="umap",
            color="clusters",
            size=20,
            alpha=0.5,
            arrow_length=2,
            arrow_size=(10, 18, 8),
            arrow_color="black",
            density=0.8,
            show=False,
            ax=axs_top[i, 0],
            title=method,
        )
        axs_top[i, 0].axis("equal")
    # handles, labels = axs_top[0].get_legend_handles_labels()
    axs_top[-1, 0].legend(ncols=2, loc="upper center", bbox_to_anchor=(0.5, 0.0))

    # axs_center = sub_figs[1].subplots(4, 1).flat
    for i, (method, adata) in enumerate(adatas_stewart.items()):
        scv.pl.velocity_embedding_grid(
            adata,
            basis="umap",
            color="cell_type",
            size=20,
            alpha=0.5,
            arrow_length=2,
            arrow_size=(10, 18, 8),
            arrow_color="black",
            density=0.8,
            show=False,
            ax=axs_top[i, 1],
            title=method,
        )
        axs_top[i, 1].axis("equal")
    # handles, labels = axs_center[0].get_legend_handles_labels()
    axs_top[-1, 1].legend(ncols=2, loc="upper center", bbox_to_anchor=(0.5, 0.0))

    sub_figs_bottom = sub_figs[1].subfigures(3, 1, height_ratios=[2, 2, 3])
    split_boxplot(
        sub_figs_bottom[0],
        velocities_cosine,
        "cosine_similarity",
        y_label="cosine similarity",
        hue="comparison",
        palette=PALETTE_COMP,
    )
    split_boxplot(
        sub_figs_bottom[1],
        velocities_corr,
        "pearsonr",
        y_label="Pearson corr.",
        hue="comparison",
        palette=PALETTE_COMP,
    )
    handles, labels = sub_figs_bottom[0].get_axes()[0].get_legend_handles_labels()
    sub_figs_bottom[2].legend(
        handles,
        [lab.replace("_", " vs. ") for lab in labels],
        loc="upper center",
    )

    labx, laby = 0.02, 0.995
    for i, ax in enumerate(axs_top[0, :]):
        fig.text(
            labx,
            laby,
            string.ascii_uppercase[i],
            transform=mpl.transforms.blended_transform_factory(
                ax.transAxes, fig.transFigure
            ),
            **AXLAB_KWS,
        )
    sub_figs_bottom[0].text(
        labx,
        laby,
        "C",
        transform=mpl.transforms.blended_transform_factory(
            sub_figs[1].transSubfigure, fig.transFigure
        ),
        **AXLAB_KWS,
    )
    sub_figs_bottom[1].text(
        labx,
        laby,
        "D",
        transform=mpl.transforms.blended_transform_factory(
            sub_figs[1].transSubfigure, sub_figs_bottom[1].transSubfigure
        ),
        **AXLAB_KWS,
    )
    fig.savefig(os.path.join(FIG_DIR, "fig5"))
    plt.close()


def supplementary_figure_3():
    # Load all DataFrames
    velocities_cosine = read_dataframes("velocities_cosine")
    comp_mapping = {
        "velocyto_tidesurf": "tidesurf_velocyto",
        "alevin-fry_tidesurf": "tidesurf_alevin-fry",
        "velocyto_alevin-fry": "velocyto_alevin-fry",
    }
    velocities_cosine["comparison"] = velocities_cosine["comparison"].map(comp_mapping)
    velocities_corr = read_dataframes("velocities_corr")
    velocities_corr["comparison"] = velocities_corr["comparison"].map(comp_mapping)

    # Make figure
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_WIDTH_SINGLE / 1.5), layout="constrained")
    sub_figs = fig.subfigures(1, 3)

    split_boxplot(
        sub_figs[0],
        velocities_cosine,
        "cosine_similarity",
        y_label="cosine similarity",
        hue="comparison",
        palette="Accent",
    )
    split_boxplot(
        sub_figs[1],
        velocities_corr,
        "pearsonr",
        y_label="Pearson corr.",
        hue="comparison",
        palette="Accent",
    )
    handles, labels = sub_figs[0].get_axes()[0].get_legend_handles_labels()
    sub_figs[2].legend(
        handles,
        [lab.replace("_", " vs.\n") for lab in labels],
        loc="center",
    )
    for i, label in enumerate(string.ascii_uppercase[:2]):
        sub_figs[i].text(
            0.02, 0.98, label, transform=sub_figs[i].transSubfigure, **AXLAB_KWS
        )
    fig.savefig(os.path.join(FIG_DIR, "suppfig3.pdf"))
    plt.close()


def supplementary_figure_4():
    # Load runtimes DataFrame
    run_metrics = pd.read_csv("../out/runtimes_cluster.csv", index_col=0)
    run_metrics["runtime"] = run_metrics["runtime"].astype("timedelta64[s]")
    run_metrics["runtime_hours"] = run_metrics["runtime"].dt.total_seconds() / 3600

    # Make figure
    fig = plt.figure(figsize=(FIG_WIDTH / 1.5, FIG_WIDTH / 3), layout="constrained")
    sub_figs = fig.subfigures(1, 2)
    axs = [sub_fig.subplots(1, 1) for sub_fig in sub_figs]

    sns.lineplot(
        run_metrics,
        x="n_reads",
        y="runtime_hours",
        hue="method",
        hue_order=["velocyto", "tidesurf"],
        palette=[
            mpl.cm.get_cmap("Dark2").colors[0],
            mpl.cm.get_cmap("Dark2").colors[2],
        ],
        err_style="bars",
        errorbar="sd",
        ax=axs[0],
    )
    axs[0].set_xlabel("number of reads [millions]")
    axs[0].xaxis.set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, _: f"{x / 1e6:.0f}")
    )
    axs[0].set_ylabel("run time [hours]")
    axs[0].legend(loc="upper left", markerscale=1.2)

    sns.lineplot(
        run_metrics,
        x="n_reads",
        y="memory",
        hue="method",
        hue_order=["velocyto", "tidesurf"],
        palette=[
            mpl.cm.get_cmap("Dark2").colors[0],
            mpl.cm.get_cmap("Dark2").colors[2],
        ],
        err_style="bars",
        errorbar="sd",
        ax=axs[1],
    )

    axs[1].get_legend().remove()
    axs[1].set_xlabel("number of reads [millions]")
    axs[1].xaxis.set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, _: f"{x / 1e6:.0f}")
    )
    axs[1].set_ylabel("memory [GB]")
    for sub_fig, lab in zip(sub_figs, string.ascii_uppercase[: len(sub_figs)]):
        sub_fig.text(0.015, 0.98, lab, **AXLAB_KWS)

    fig.savefig(os.path.join(FIG_DIR, "suppfig4.pdf"))
    plt.close()


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    figure_1()
    supplementary_figure_1()
    supplementary_figure_2()
    figure_3()
    figure_4()
    supplementary_figure_3()
    figure_5()
    supplementary_figure_4()


if __name__ == "__main__":
    main()
