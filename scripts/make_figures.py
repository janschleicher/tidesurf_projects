import os
import string
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns
import scanpy as sc
import scvelo as scv
from typing import Optional


# Set up plotting parameters
sc.set_figure_params(vector_friendly=True, dpi_save=300, scanpy=False, fontsize=10)
scv.set_figure_params(vector_friendly=True, dpi_save=300, fontsize=10)
plt.style.use("matplotlibrc.txt")
AXLAB_KWS = {"fontsize": 12, "fontweight": "bold", "va": "top", "ha": "left"}
FIG_WIDTH = 17.8 / 2.54  # in inches
FIG_WIDTH_SINGLE = 8.6 / 2.54  # in inches
FIG_DIR = "../figures"

DATASETS_ORDER = ["pancreas", "retina", "stewart", "fu"]
DATASETS_LABELS = {
    "pancreas": "Pancreas",
    "retina": "Retina",
    "stewart": "B cells",
    "fu": "T cells",
}
METHODS_ORDER = [
    "velocyto",
    "alevin-fry",
    "tidesurf",
    "cellranger",
    "tidesurf_velocyto",
    "tidesurf_alevin-fry",
    "velocyto_alevin-fry",
]
SEQ_TYPES = ["3'", "5'"]
SPLICE_STATES = ["spliced", "unspliced", "ambiguous"]


def read_dataframes(name: str):
    df = (
        pd.concat(
            {
                dataset: pd.read_parquet(f"../out/{dataset}/{name}.parquet.gz")
                for dataset in DATASETS_ORDER
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


def figure_2():
    # Load all DataFrames
    total_counts = read_dataframes("total_counts")
    counts_diff = read_dataframes("counts_diff")
    counts_corr = read_dataframes("counts_corr_cellranger")
    counts_cosine = read_dataframes("counts_cosine_cellranger")

    # Make figure
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_WIDTH / 3), layout="constrained")
    sub_figs = fig.subfigures(1, 4)
    split_boxplot(sub_figs[0], total_counts[total_counts["genes"] == "all"], "counts")
    split_boxplot(sub_figs[1], counts_diff[counts_diff["genes"] == "all"], "difference")
    for ax in sub_figs[1].get_axes():
        ax.axhline(0, color="grey", linestyle="--", zorder=0)
    split_boxplot(sub_figs[2], counts_corr, "pearsonr", y_label="Pearson corr.")
    split_boxplot(
        sub_figs[3], counts_cosine, "cosine_similarity", y_label="cosine similarity"
    )
    handles, labels = sub_figs[0].get_axes()[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        ncol=4,
        loc="outside lower center",
    )
    for i, label in enumerate(string.ascii_uppercase[:4]):
        sub_figs[i].text(
            0.02, 0.98, label, transform=sub_figs[i].transSubfigure, **AXLAB_KWS
        )
    fig.savefig(os.path.join(FIG_DIR, "fig2.pdf"))
    plt.close()


def figure_3():
    # Load all DataFrames
    spliced_unspliced_counts = read_dataframes("spliced_unspliced_counts")
    spliced_unspliced_ratios = read_dataframes("spliced_unspliced_ratios")
    spliced_unspliced_cosine = read_dataframes("spliced_unspliced_cosine")

    # Make figure
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_WIDTH), layout="constrained")
    sub_figs = fig.subfigures(3, 3).ravel()

    # Total spliced/unspliced/ambiguous counts per cell
    for sub_fig, splice_state in zip(sub_figs[:3], SPLICE_STATES):
        split_boxplot(
            sub_fig,
            spliced_unspliced_counts[
                (spliced_unspliced_counts["genes"] == "all")
                & (spliced_unspliced_counts["splice_state"] == splice_state)
            ],
            "counts",
        )
        sub_fig.suptitle(splice_state)
        if splice_state != SPLICE_STATES[0]:
            sub_fig.get_axes()[0].set_ylabel("")

    # Total spliced/unspliced/ambiguous counts ratios per cell
    for sub_fig, splice_state in zip(sub_figs[3:6], SPLICE_STATES):
        split_boxplot(
            sub_fig,
            spliced_unspliced_ratios[
                (spliced_unspliced_ratios["genes"] == "all")
                & (spliced_unspliced_ratios["splice_state"] == splice_state)
            ],
            "ratio",
        )
        if splice_state != SPLICE_STATES[0]:
            sub_fig.get_axes()[0].set_ylabel("")

    # Cosine similarity between methods
    for sub_fig, splice_state in zip(sub_figs[6:], SPLICE_STATES):
        split_boxplot(
            sub_fig,
            spliced_unspliced_cosine[
                spliced_unspliced_cosine["splice_state"] == splice_state
            ],
            "cosine_similarity",
            y_label="cosine similarity",
            hue="comparison",
            palette="Accent",
        )
        if splice_state != SPLICE_STATES[0]:
            sub_fig.get_axes()[0].set_ylabel("")

    handles_top, labels_top = sub_figs[0].get_axes()[0].get_legend_handles_labels()
    handles_bottom, labels_bottom = (
        sub_figs[6].get_axes()[0].get_legend_handles_labels()
    )
    fig.legend(
        handles_top + handles_bottom,
        [lab.replace("_", " vs.\n") for lab in labels_top + labels_bottom],
        ncol=6,
        loc="outside lower center",
    )
    for i, label in enumerate(string.ascii_uppercase[:9]):
        sub_figs[i].text(
            0.02, 0.99, label, transform=sub_figs[i].transSubfigure, **AXLAB_KWS
        )
    fig.savefig(os.path.join(FIG_DIR, "fig3.pdf"))
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
        ncol=5,
        loc="outside lower center",
    )
    for i, label in enumerate(string.ascii_uppercase[:6]):
        sub_figs[i].text(
            0.02, 0.98, label, transform=sub_figs[i].transSubfigure, **AXLAB_KWS
        )
    fig.savefig(os.path.join(FIG_DIR, "suppfig1.pdf"))
    plt.close()


def figure_4():
    # Load all DataFrames
    velocities_cosine = read_dataframes("velocities_velo_cosine")
    comp_mapping = {
        "velocyto_tidesurf": "tidesurf_velocyto",
        "alevin-fry_tidesurf": "tidesurf_alevin-fry",
        "velocyto_alevin-fry": "velocyto_alevin-fry",
    }
    velocities_cosine["comparison"] = velocities_cosine["comparison"].map(comp_mapping)
    velocities_corr = read_dataframes("velocities_velo_corr")
    velocities_corr["comparison"] = velocities_corr["comparison"].map(comp_mapping)

    # Load AnnData objects
    adatas = {
        method: sc.read_h5ad(f"../out/pancreas/adata_{method}.h5ad")
        for method in METHODS_ORDER[:3]
    }

    # Make figure
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_WIDTH / 1.5), layout="constrained")
    sub_figs = fig.subfigures(2, 1, height_ratios=[1.2, 1])
    axs_top = sub_figs[0].subplots(1, 3)
    for i, (method, adata) in enumerate(adatas.items()):
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
            ax=axs_top[i],
            title=method,
        )
    handles, labels = axs_top[0].get_legend_handles_labels()
    sub_figs[0].legend(handles, labels, ncol=4, loc="outside lower center")
    sub_figs_bottom = sub_figs[1].subfigures(1, 3)
    split_boxplot(
        sub_figs_bottom[0],
        velocities_cosine,
        "cosine_similarity",
        y_label="cosine similarity",
        hue="comparison",
        palette="Accent",
    )
    split_boxplot(
        sub_figs_bottom[1],
        velocities_corr,
        "pearsonr",
        y_label="Pearson corr.",
        hue="comparison",
        palette="Accent",
    )
    handles, labels = sub_figs_bottom[0].get_axes()[0].get_legend_handles_labels()
    sub_figs_bottom[2].legend(
        handles,
        [lab.replace("_", " vs.\n") for lab in labels],
        loc="center",
    )
    for sub_fig, label in zip(
        [sub_figs[0], sub_figs_bottom[0], sub_figs_bottom[1]],
        string.ascii_uppercase[:3],
    ):
        sub_fig.text(0.01, 0.99, label, transform=sub_fig.transSubfigure, **AXLAB_KWS)
    fig.savefig(os.path.join(FIG_DIR, "fig4.pdf"))
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
    figure_2()
    figure_3()
    supplementary_figure_1()
    figure_4()
    supplementary_figure_2()
    supplementary_figure_3()
    supplementary_figure_4()


if __name__ == "__main__":
    main()
