import argparse
import scanpy as sc

def main():
    parser = argparse.ArgumentParser(description='Process tissue adata from pipeline with cellbender.')
    parser.add_argument('-t', '--tissue', type=str, help='Tissue', required=True)
    parser.add_argument('-c', '--clusters', type=str, nargs='+', help='List of clusters to sub-cluster', required=True)
    parser.add_argument('-r', '--resolution', type=float, help='Leiden clustering resolution', required=True)

    args = parser.parse_args()
        
    print(f"Resolution: {args.resolution}")

    if args.tissue and args.clusters:
        process_file(args.tissue, args.clusters, args.resolution)
    else:
        print("Please provide a tissue using -t option and clusters using -c option.")

def process_file(tissue, clusters, res):
    print("Processing tissue:", tissue)
    
    # Read the .h5ad file
    adata = sc.read_h5ad(f'../cellbender_tissues/processed/{tissue}_processed.h5ad')
    
    for cluster in clusters:
        print(f"Sub-clustering for cluster: {cluster}")

        # Determine the clustering base column
        base_column = 'leiden' if 'leiden_R' not in adata.obs.columns else 'leiden_R'

        # Perform Leiden clustering on the current cluster
        sc.tl.leiden(adata, restrict_to=(base_column, [cluster]), resolution=res)

    # Save the processed file
    output_path = f'../cellbender_tissues/processed/{tissue}_processed_subclustered_res{res}.h5ad'
    adata.obs['leiden_R'] = adata.obs['leiden_R'].str.replace(',', '_')
    adata.write_h5ad(output_path)
    
    print(f"Processed data saved to: {output_path}")

if __name__ == "__main__":
    main()
