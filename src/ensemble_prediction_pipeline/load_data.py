import pandas as pd
import argparse
import subprocess
import sys
from tqdm import tqdm
import urllib.request


def main():

    args = sys.argv[1:]
    data_source = pd.read_csv(args[0], sep="\t")

    print(data_source.head())

    df_dict = dict()
    for _, df_row in tqdm(data_source.iterrows()):
        file_name = "data/raw/" + df_row["dataset"] + ".csv"
        urllib.request.urlretrieve(df_row["Figshare_link"], file_name)
        # process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        df_dict[df_row["dataset"]] = file_name

    df = pd.DataFrame(df_dict, index=[0]).transpose().reset_index()
    df.columns = ["dataset", "filename"]

    df.to_csv(args[2])


if __name__ == "__main__":
    main()
