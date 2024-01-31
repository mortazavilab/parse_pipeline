# This script makes knee plots based on the pipeline output with 200umi cutoff

import scanpy as sc
import pandas as pd
import anndata as ad
import argparse
import numpy as np 
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description='Process some files.')
    parser.add_argument('-f', '--file', type=str, help='input file path')

    args = parser.parse_args()

    if args.file:
        process_file(args.file)
    else:
        print("Please provide a file using -f option.")

def process_file(file_path):
    
    adata = sc.read_h5ad('../preprocessed_tissues/' + file_path)


    print("Processing file:", file_path)
    knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.loglog(knee, range(len(knee)), linewidth=5, color="g")
    # ax.axvline(x=knee[expected_num_cells], linewidth=3, color="k")
    # ax.axhline(y=expected_num_cells, linewidth=3, color="k")

    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Set of Barcodes")

    plt.grid(True, which="both")
    

    filename = file_path
    filename2 = filename.replace("preprocessed.h5ad", "")
    print(filename2)
    fig.savefig('../kneeplots/' + filename2 + '_knee_plot.png')

if __name__ == "__main__":
    main()


