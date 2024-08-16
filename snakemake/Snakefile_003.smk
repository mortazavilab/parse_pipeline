import pandas as pd

p = os.path.dirname(os.getcwd())
sys.path.append(p)

from utils import *

######## Only need to edit this part ########
config_tsv = 'configs/igvf_003_config.tsv'
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
    mult_genotype_2='|'.join([re.escape(x) for x in mult_genotype_2s]),
    total_drops='|'.join([re.escape(str(x)) for x in df.droplets_included.tolist()])

def get_subset_tissues(df, sample_df):
    temp = df.merge(sample_df, on='plate', how='inner')
    tissues = temp.Tissue.unique().tolist()
    return tissues

rule all:
    input:
        expand(config['tissue']['adata'],
                plate=df.plate.tolist(),
                tissue=get_subset_tissues(df, sample_df))
                        
################################################################################
################################### cellbender ###################################
################################################################################
rule cellbender:
    input:
        unfilt_adata = config['kallisto']['unfilt_adata']
    params:
        total_drops = lambda wildcards: df[df['subpool'] == wildcards.subpool]['droplets_included'].values[0],
        learning_rate = lambda wildcards: df[df['subpool'] == wildcards.subpool]['learning_rate'].values[0],
    resources:
        mem_gb = 250,
        threads = 12
    output:
        filt_h5 = config['cellbender']['filt_h5'],
    conda:
        "cellbender"
    shell:
        """
        mkdir -p $(dirname {output.filt_h5})
        cd $(dirname {output.filt_h5})
        
        export CUDA_VISIBLE_DEVICES=0,1  # jaz tip!

        # Run cb in target directory
        # trying to make sure the checkpoint isn't overwritten when multiple CB run in parallel
        cellbender remove-background \
            --input {input.unfilt_adata} \
            --output {output.filt_h5} \
            --total-droplets-included {params.total_drops} \
            --learning-rate {params.learning_rate} \
            --cuda
        """

rule make_filt_adata:
    resources:
        mem_gb = 128,
        threads = 4
    input:
        filt_h5 = config['cellbender']['filt_h5'],
        unfilt_adata = config['kallisto']['unfilt_adata'],
    output:
        filt_adata = config['cellbender']['filt_adata']
    run:
        add_meta_filter(input.filt_h5,
                        input.unfilt_adata,
                        wildcards,
                        bc_df,
                        kit,
                        chemistry,
                        sample_df,
                        output.filt_adata)

################################################################################
############################### Merge klue output ##############################
################################################################################

rule klue_merge_genotype:
    input:
        genotype_counts = config['klue']['genotype_counts'],
        adata = config['cellbender']['filt_adata']
    resources:
        mem_gb = 128,
        threads = 2
    output:
        adata = config['cellbender']['genotype_adata']
    run:
        merge_kallisto_klue(input.adata,
                            input.genotype_counts,
                            output.adata)

####################
###### Scrublet ####
####################
rule make_subpool_sample_adata:
    input:
        adata = config['cellbender']['genotype_adata']
    resources:
        mem_gb = 32,
        threads = 2
    output:
        adata = config['scrublet']['adata']
    run:
        make_subpool_sample_adata(input.adata,
                                  wildcards,
                                  output.adata)

rule scrublet:
    input:
        adata = config['scrublet']['adata']
    params:
        n_pcs = 30,
        min_cells = 1,
        min_gene_variability_pctl = 85
    resources:
        mem_gb = 256,
        threads = 8
    output:
        adata = config['scrublet']['scrub_adata']
    run:
        run_scrublet(input.adata,
                     params.n_pcs,
                     params.min_cells,
                     params.min_gene_variability_pctl,
                     output.adata)


####################
### Combine adatas
###################

def get_tissue_adatas(df, sample_df, wc, cfg_entry):

    # limit to input tissue
    temp_sample = sample_df.copy(deep=True)
    temp_sample = temp_sample.loc[temp_sample.Tissue==wc.tissue]

    # merge this stuff in with the fastq df
    fastq_df = df.copy(deep=True)
    temp = fastq_df.merge(temp_sample, on='plate', how='inner')

    # get the plate / subpool / sample info for this tissue
    plates = temp.plate.tolist()
    subpools = temp.subpool.tolist()
    samples = temp.Mouse_Tissue_ID.tolist()

    files = expand(cfg_entry,
           zip,
           plate=plates,
           sample=samples,
           subpool=subpools)
    files = list(set(files))

    return files


# TODO - make one of the concatenation rules for plate+tissue
# and make another for just tissue
rule make_tissue_adata:
    input:
        adatas = lambda wc:get_tissue_adatas(df, sample_df, wc, config['scrublet']['scrub_adata'])
    resources:
        mem_gb = 256,
        threads = 2
    output:
        adata = config['tissue']['adata']
    run:
        concat_adatas(input.adatas, output.adata)
