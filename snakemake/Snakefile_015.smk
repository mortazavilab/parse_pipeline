import pandas as pd

p = os.path.dirname(os.getcwd())
sys.path.append(p)

from utils import *

######## Only need to edit this part ########
config_tsv = 'configs/igvf_015_config.tsv'
sample_csv = 'configs/sample_metadata.csv'

kit = 'WT_mega'  # either WT (48 wells), WT_mini (12 wells), or WT_mega (96 wells)
chemistry = 'v2'  # all IGVF and ModelAD experiments are v2 so far, v3 coming soon

######## Do not change anything past this point ########

configfile: 'configs/config.yml'

# read in config / analysis spec
df = parse_config(config_tsv)
bc_df = get_bc1_matches(kit, chemistry)
sample_df = parse_sample_df(sample_csv)

mult_genotype_1s = sample_df.loc[sample_df.mult_genotype_1.notnull(), 'mult_genotype_1'].unique().tolist()
mult_genotype_2s = sample_df.loc[sample_df.mult_genotype_2.notnull(), 'mult_genotype_2'].unique().tolist()

wildcard_constraints:
    plate='|'.join([re.escape(x) for x in df.plate.tolist()]),
    subpool='|'.join([re.escape(x) for x in df.subpool.tolist()]),
    lane='|'.join([re.escape(x) for x in df.lane.tolist()]),
    sample='|'.join([re.escape(x) for x in sample_df.Mouse_Tissue_ID.tolist()]),
    tissue='|'.join([re.escape(x) for x in sample_df.Tissue.tolist()]),
    genotype='|'.join([re.escape(x) for x in sample_df.Genotype.unique().tolist()]),
    mult_genotype_1='|'.join([re.escape(x) for x in mult_genotype_1s]),
    mult_genotype_2='|'.join([re.escape(x) for x in mult_genotype_2s])

def get_subset_tissues(df, sample_df):
    temp = df.merge(sample_df, on='plate', how='inner')
    tissues = temp.Tissue.unique().tolist()
    return tissues
    
rule all:
    input:
        expand(config['tissue']['adata'],
               plate=df.plate.tolist(),
               tissue=get_subset_tissues(df, sample_df)),
        expand(config['tissue']['adata_raw_counts'],
               plate=df.plate.tolist(),
               tissue=get_subset_tissues(df, sample_df)),
        expand(config['tissue']['adata_cb_counts'],
               plate=df.plate.tolist(),
               tissue=get_subset_tissues(df, sample_df))
                        
################################################################################
################################## cellbender ##################################
################################################################################
rule cellbender:
    input:
        unfilt_adata = config['kallisto']['unfilt_adata']
    params:
        total_drops = lambda wildcards: df[df['subpool'] == wildcards.subpool]['droplets_included'].values[0],
        learning_rate = lambda wildcards: df[df['subpool'] == wildcards.subpool]['learning_rate'].values[0],
        expected_cells = lambda wildcards: df[df['subpool'] == wildcards.subpool]['expected_cells'].values[0],
    resources:
        mem_gb = 250,
        threads = 12
    output:
        filt_h5 = config['cellbender']['filt_h5'],
        unfilt_h5 = config['cellbender']['unfilt_h5'],
    conda:
        "cellbender"
    shell:
        """
        mkdir -p $(dirname {output.filt_h5})
        cd $(dirname {output.filt_h5})
        
        #SUBPOOL_ID=$(echo {wildcards.subpool} | grep -o '[0-9]*')
        #GPU_ID=$(( SUBPOOL_ID % 2 ))
        #export CUDA_VISIBLE_DEVICES=$GPU_ID

        # Run cb in target directory
        cellbender remove-background \
            --input {input.unfilt_adata} \
            --output {output.unfilt_h5} \
            --total-droplets-included {params.total_drops} \
            --learning-rate {params.learning_rate} \
            --expected-cells {params.expected_cells} \
            --cuda
        """ 

################################################################################
##################### Merge klue results and run scrublet ######################
################################################################################

rule make_filt_adata:
    resources:
        mem_gb = 128,
        threads = 4
    input:
        filt_h5 = config['cellbender']['filt_h5'],
        unfilt_adata = config['kallisto']['unfilt_adata'],
        genotype_counts = config['klue']['genotype_counts']
    output:
        filt_adata = config['cellbender']['filt_adata']
    run:
        add_meta_filter(input.filt_h5,
                        input.unfilt_adata,
                        input.genotype_counts,
                        wildcards,
                        bc_df,
                        kit,
                        chemistry,
                        sample_df,
                        output.filt_adata)

################################################################################
################################ Combine adatas ################################
################################################################################

def get_adatas(df, cfg_entry):
    # Extract unique plates and subpools from the DataFrame
    plates = df.plate.unique().tolist()
    subpools = df.subpool.unique().tolist()
    
    # Generate file paths based on plates and subpools
    files = expand(cfg_entry,
                   plate=plates,
                   subpool=subpools)
    
    # Remove duplicate file paths
    files = list(set(files))
    
    return files

######

rule make_plate_adata_raw_counts:
    input:
        adatas = lambda wc: get_adatas(df, config['cellbender']['filt_adata'])
    resources:
        mem_gb = 350,
        threads = 2
    output:
        adata = config['plate']['adata_raw_counts']
    run:
        concat_adatas_raw_counts(input.adatas, output.adata)


rule make_plate_adata_cb_counts:
    input:
        adatas = lambda wc: get_adatas(df, config['cellbender']['filt_adata'])
    resources:
        mem_gb = 350,
        threads = 2
    output:
        adata = config['plate']['adata_cb_counts']
    run:
        concat_adatas_cb_counts(input.adatas, output.adata)
        
######

rule make_plate_adata:
    input:
        adata_raw_counts = config['plate']['adata_raw_counts'],
        adata_cb_counts = config['plate']['adata_cb_counts']
    resources:
        mem_gb = 350,
        threads = 2
    output:
        adata=config['plate']['adata']
    run:
        import scanpy as sc

        adata = sc.read_h5ad(input.adata_raw_counts)
        adata_cb = sc.read_h5ad(input.adata_cb_counts)

        # Add the 'cellbender_counts' layer from adata_cb to adata
        adata.layers['cellbender_counts'] = adata_cb.X.copy()     

        # Save the combined adata object to the output file
        adata.write_h5ad(output.adata)
        
######

rule make_tissue_adata:
    input:
        plate_adata = config['plate']['adata']
    resources:
        mem_gb = 350,
        threads = 2
    output:
        tissue_adata =config['tissue']['adata']
    run:
        os.makedirs(os.path.dirname(output.tissue_adata), exist_ok=True)

        adata = sc.read_h5ad(input.plate_adata)
        adata_tissue = adata[adata.obs['Tissue'] == wildcards.tissue].copy()
        adata_tissue.write_h5ad(output.tissue_adata)
        
        
rule make_tissue_adata_raw_counts:
    input:
        plate_adata = config['plate']['adata_raw_counts']
    resources:
        mem_gb = 350,
        threads = 2
    output:
        tissue_adata = config['tissue']['adata_raw_counts']
    run:
        os.makedirs(os.path.dirname(output.tissue_adata), exist_ok=True)

        adata = sc.read_h5ad(input.plate_adata)
        adata_tissue = adata[adata.obs['Tissue'] == wildcards.tissue].copy()
        adata_tissue.write_h5ad(output.tissue_adata)
        
rule make_tissue_adata_cellbender_counts:
    input:
        plate_adata = config['plate']['adata_cb_counts']
    resources:
        mem_gb = 350,
        threads = 2
    output:
        tissue_adata = config['tissue']['adata_cb_counts']
    run:
        os.makedirs(os.path.dirname(output.tissue_adata), exist_ok=True)

        adata = sc.read_h5ad(input.plate_adata)
        adata_tissue = adata[adata.obs['Tissue'] == wildcards.tissue].copy()
        adata_tissue.write_h5ad(output.tissue_adata)