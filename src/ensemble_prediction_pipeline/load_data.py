import pandas as pd
import argparse
import subprocess


def __main__():

    args = sys.argv[1:]
    data_source = pd.read_csv(args[0])

    df_dict = dict()
    for _, df_row in data_source.iterrows():
        file_name = "data/raw/" + df_row["dataset"] + ".csv"
        command = (
            "wget "
            + df_row["Figshare_ID"]
            + " -O "
            + "data/raw/"
            + df_row["dataset"]
            + ".csv"
        )
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        df_dict[df_row["dataset"]] = file_name

    df = pd.DataFrame(df_dict, index=[0]).transpose().reset_index()
    df.columns = ["dataset", "filename"]

    df.to_csv(args[2])
