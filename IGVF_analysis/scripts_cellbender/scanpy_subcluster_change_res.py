import argparse
import scanpy as sc

def main():
    parser = argparse.ArgumentParser(description='Process tissue adata from pipeline with cellbender.')
    parser.add_argument('-t', '--tissue', type=str, help='Tissue', required=True)
    parser.add_argument('-c', '--clusters', type=str, nargs='+', help='List of clusters to sub-cluster', required=True)
    parser.add_argument('-r', '--resolutions', type=float, nargs='+', help='List of Leiden clustering resolutions', required=True)

    args = parser.parse_args()
        
    if args.tissue and args.clusters and args.resolutions:
        if len(args.clusters) != len(args.resolutions):
            print("Error: The number of clusters and resolutions must be the same.")
            return
        
        process_file(args.tissue, args.clusters, args.resolutions)
    else:
        print("Please provide a tissue using -t option, clusters using -c option, and resolutions using -r option.")


def process_file(tissue, clusters, resolutions):
    print("Processing tissue:", tissue)
    
    # Read the .h5ad file
    adata = sc.read_h5ad(f'../cellbender_tissues/{tissue}_processed.h5ad')
    
    for cluster, res in zip(clusters, resolutions):
        print(f"Sub-clustering for cluster: {cluster} with resolution: {res}")

        # Determine the clustering base column
        base_column = 'leiden' if 'leiden_R' not in adata.obs.columns else 'leiden_R'

        # Perform Leiden clustering on the current cluster with the specified resolution
        sc.tl.leiden(adata, restrict_to=(base_column, [cluster]), resolution=res)

    # Save the processed file
    output_path = f'../cellbender_tissues/{tissue}_processed_subclustered.h5ad'
    adata.obs['leiden_R'] = adata.obs['leiden_R'].str.replace(',', '_')
    adata.write_h5ad(output_path)
    
    print(f"Processed data saved to: {output_path}")

if __name__ == "__main__":
    main()