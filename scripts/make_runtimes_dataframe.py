import os
from glob import glob
import re
import pandas as pd

DATASETS = ["pancreas", "retina", "stewart", "fu"]

if __name__ == "__main__":
    metrics = []
    for dataset in DATASETS:
        files = glob(
            f"../data/{dataset}/cellranger/**/metrics_summary.csv", recursive=True
        )
        for file in files:
            sample = file.split("/")[-3]
            if sample == "cellranger":
                sample = dataset

            # Load tidesurf log
            ts_log_file = (
                f"../data/{dataset}/tidesurf/{sample}/tidesurf.log"
                if sample != dataset
                else f"../data/{dataset}/tidesurf/tidesurf.log"
            )
            if os.path.isfile(ts_log_file):
                with open(ts_log_file, "r") as f:
                    content = f.read()
                    time_ts = re.search(
                        r"Finished in (\d+:\d+:\d+\.\d+)", content
                    ).group(1)
            else:
                time_ts = None

            # Load velocyto time
            velocyto_time_file = (
                f"../data/{dataset}/cellranger/{sample}_time.txt"
                if sample != dataset
                else f"../data/{dataset}/cellranger/time.txt"
            )
            if os.path.isfile(velocyto_time_file):
                with open(velocyto_time_file, "r") as f:
                    content = f.read()
                    time_velo = re.search(r"real\t(\d+m\d+.\d+s)", content).group(1)
            else:
                time_velo = None

            # Load Cell Ranger metrics
            metrics_summary = pd.read_csv(file)
            metrics.append(
                {
                    "dataset": dataset,
                    "sample": sample,
                    "n_reads": metrics_summary["Number of Reads"].values[0],
                    "n_cells": metrics_summary["Estimated Number of Cells"].values[0],
                    "tidesurf": time_ts,
                    "velocyto": time_velo,
                }
            )
    df = pd.DataFrame(metrics)
    df["n_reads"] = df["n_reads"].str.replace(",", "").astype(int)
    df["n_cells"] = df["n_cells"].astype("str").str.replace(",", "").astype(int)
    df["tidesurf"] = df["tidesurf"].astype("timedelta64[s]")
    df["velocyto"] = df["velocyto"].astype("timedelta64[s]")
    df = df.melt(
        id_vars=["dataset", "sample", "n_reads", "n_cells"],
        value_vars=["tidesurf", "velocyto"],
        var_name="method",
        value_name="runtime",
    )
    pd.set_option("display.max_columns", None)
    df.to_csv("../out/runtimes.csv")
