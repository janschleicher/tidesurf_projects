import pandas as pd

from make_figures import read_dataframes

if __name__ == "__main__":
    pd.set_option("display.max_columns", None)

    # Total counts and counts diff
    for df_name in ["total_counts", "counts_diff"]:
        df = read_dataframes(df_name)
        print(f"==== {df_name.replace('_', ' ').capitalize()} ====")
        print(df[df["genes"] == "all"].groupby(["method", "dataset"]).describe())
        print("\n")

    # Correlation with Cell Ranger
    counts_corr = read_dataframes("counts_corr_cellranger")
    print("==== Corr. with Cell Ranger ====")
    print(counts_corr.groupby(["method", "dataset"]).describe())
    print("\n")

    # Cosine similarity with Cell Ranger
    counts_cosine = read_dataframes("counts_cosine_cellranger")
    print("==== Cosine similarity with Cell Ranger ====")
    print(counts_cosine.groupby(["method", "dataset"]).describe())
    print("\n")

    # Spliced and unspliced counts and ratios
    for df_name in [
        "spliced_unspliced_counts",
        "spliced_unspliced_ratios",
        "spliced_unspliced_diff",
    ]:
        df = read_dataframes(df_name)
        print(f"==== {df_name.replace('_', ' ').capitalize()} ====")
        print(
            df[df["genes"] == "all"]
            .groupby(["method", "dataset", "splice_state"])
            .describe()
        )
        print("\n")

    # Spliced and unspliced cosine similarity and correlation
    for df_name in ["spliced_unspliced_cosine", "spliced_unspliced_corr"]:
        df = read_dataframes(df_name)
        print(f"==== {df_name.replace('_', ' ').capitalize()} ====")
        print(df.groupby(["comparison", "dataset", "splice_state"]).describe())
        print("\n")

    # Velocities
    for df_name in [
        "velocities_cosine",
        "velocities_corr",
        "velocities_velo_cosine",
        "velocities_velo_corr",
    ]:
        df = read_dataframes(df_name)
        print(f"==== {df_name.replace('_', ' ').capitalize()} ====")
        print(df.groupby(["comparison", "dataset"]).describe())
        print("\n")
