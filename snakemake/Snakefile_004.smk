import pandas as pd

p = os.path.dirname(os.getcwd())
sys.path.append(p)

from utils import *

######## Only need to edit this part ########
config_tsv = 'configs/igvf_004_config.tsv'
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

include: "klue.smk"

rule all:
    input:
        expand(config['ref']['klue']['ind'],
               zip,
               mult_genotype_1=[g for g in mult_genotype_1s if g in get_founder_genotypes()],
               mult_genotype_2=[g for g in mult_genotype_2s if g in get_founder_genotypes()]),
        #expand(config['ref']['genome']['fa'],
        #        genotype=get_founder_genotypes()),
        expand(config['klue']['genotype_counts'],
                zip,
                plate=df.plate.tolist(),
                subpool=df.subpool.tolist()),
        expand(config['kallisto']['unfilt_adata'],
                zip,
                plate=df.plate.tolist(),
                subpool=df.subpool.tolist()),
        


################################################################################
########################## Ref download and generation #########################
################################################################################
def get_fa_link(wc, config):
    genotype = wc.genotype
    
    if genotype in get_f1_genotypes():
        d = get_f1_founder_genotype_dict()
        genotype = d[genotype]
        
    link = config['ref']['genome']['link'][genotype]
    
    return link

rule curl_fa:
    resources:
        mem_gb = 4,
        threads = 1
    params:
        link = lambda wc: get_fa_link(wc, config),
    output:
        zip = temporary("{genotype}.zip")
    shell:
        """
        if [ "{wildcards.genotype}" == "B6J" ]; then
            wget -O {wildcards.genotype}.fa.gz {params.link}
            mkdir -p {wildcards.genotype}/ncbi_dataset/data/temp/
            gunzip {wildcards.genotype}.fa.gz
            mv {wildcards.genotype}.fa {wildcards.genotype}/ncbi_dataset/data/temp/temp.fna
            zip -r {wildcards.genotype}.zip {wildcards.genotype}
        else
            curl -OJX GET "{params.link}{output.zip}"
        fi
        """

rule fa_ref_fmt:
    input:
        zip = "{genotype}.zip"
    resources:
        threads = 4,
        mem_gb =16
    output:
        fa = "ref/genomes/{genotype}.fa.gz"
    shell:
        """
        unzip {input.zip} -d {wildcards.genotype}
        gzip -cvf {wildcards.genotype}/ncbi_dataset/data/*/*fna > {output.fa}
        rm -r {wildcards.genotype}
        """

rule dl:
   resources:
       mem_gb = 4,
       threads = 1
   shell:
       "wget -O {output.out} {params.link}"

use rule dl as dl_annot with:
    params:
        link = config['ref']['annot_link']
    output:
        out = config['ref']['annot']

use rule dl as dl_fa with:
    params:
        link = config['ref']['fa_link']
    output:
        out = config['ref']['fa']

rule kallisto_ind:
    input:
        annot = config['ref']['annot'],
        fa = config['ref']['fa']
    resources:
        mem_gb = 64,
        threads = 12
    output:
        t2g = config['ref']['kallisto']['t2g'],
        ind = config['ref']['kallisto']['ind'],
        fa = config['ref']['kallisto']['fa'],
        na = config['ref']['kallisto']['na'],
        c1 = config['ref']['kallisto']['c1'],
        c2 = config['ref']['kallisto']['c2']
    shell:
        """
        kb ref \
            --workflow=nac \
            -i {output.ind} \
            -g {output.t2g} \
            -c1 {output.c1} \
            -c2 {output.c2} \
            -f1 {output.fa} \
            -f2 {output.na} \
            {input.fa} \
            {input.annot}
        """

################################################################################
########################### Symlink fastqs ###########################
################################################################################

rule symlink_fastq_r1:
    params:
        fastq = lambda wc:get_df_info(wc, df, 'fastq')
    resources:
        mem_gb = 4,
        threads = 1
    output:
        fastq = config['raw']['fastq_r1']
    shell:
        """
        ln -s {params.fastq} {output.fastq}
        """

rule symlink_fastq_r2:
    params:
        fastq = lambda wc:get_df_info(wc, df, 'fastq_r2')
    resources:
        mem_gb = 4,
        threads = 1
    output:
        fastq = config['raw']['fastq_r2']
    shell:
        """
        ln -s {params.fastq} {output.fastq}
        """

################################################################################
################################### kallisto ###################################
################################################################################
rule kallisto:
    input:
        fastq_r1 = lambda wc:get_subpool_fastqs(wc, df, config, how='list', read='R1'),
        fastq_r2 = lambda wc:get_subpool_fastqs(wc, df, config, how='list', read='R2'),
        t2g = config['ref']['kallisto']['t2g'],
        ind = config['ref']['kallisto']['ind']
    params:
        # TODO bc1 map, barcodes, c1, c2 should be output from sth, seqspec
        bc1_map = config['ref']['bc1_map'],
        barcodes = config['ref']['barcodes'],
        c1 = config['ref']['kallisto']['c1'],
        c2 = config['ref']['kallisto']['c2'],
        fastq_str = lambda wc:get_subpool_fastqs(wc, df, config, how='str'),
        odir = config['kallisto']['cgb'].split('counts_unfiltered_modified/')[0]
    resources:
        mem_gb = 250,
        threads = 12
    output:
        config['kallisto']['cgb'],
        config['kallisto']['cggn'],
        config['kallisto']['cgg'],
        config['kallisto']['mtx'],
        temporary(config['kallisto']['bus']),
        temporary(config['kallisto']['bus_modified_unfilt']),
        temporary(config['kallisto']['bus_unfilt']) 
    shell:
        """
        kb count \
            --h5ad \
            --gene-names \
            --sum=total \
            --strand=forward \
            -r {params.bc1_map} \
            -w {params.barcodes} \
            --workflow=nac \
            -g {input.t2g} \
            -x SPLIT-SEQ \
            -i {input.ind} \
            -t {resources.threads} \
            -o {params.odir} \
            -c1 {params.c1} \
            -c2 {params.c2} \
            --verbose \
            {params.fastq_str}
        """

rule make_unfilt_adata:
    resources:
        mem_gb = 128,
        threads = 4
    input:
        mtx = config['kallisto']['mtx'],
        cgb = config['kallisto']['cgb'],
        cggn = config['kallisto']['cggn'],
        cgg = config['kallisto']['cgg']
    output:
        unfilt_adata = config['kallisto']['unfilt_adata']
    run:
        make_adata_from_kallisto(input.mtx,
                        input.cgb,
                        input.cggn,
                        input.cgg,
                        wildcards,
                        output.unfilt_adata)
                       

################################################################################
##################################### klue #####################################
################################################################################

# add metadata and perform basic filtering
rule make_adata_klue:
    resources:
        mem_gb = 128,
        threads = 4
    input:
        mtx = config['klue']['mtx'],
        cgb = config['klue']['cgb'],
        cggn = config['klue']['cggn'],
        cgg = config['klue']['cgg']
    output:
        adata = config['klue']['unfilt_adata']
    run:
        add_meta_klue(input.mtx,
                        input.cgb,
                        input.cggn,
                        input.cgg,
                        wildcards,
                        bc_df,
                        kit,
                        chemistry,
                        sample_df,
                        output.adata)
        
def get_subpool_adatas(df, sample_df, wc, cfg_entry):
    """
    Get adatas that belong to the same subpool across the
    different pairs of genotypes ({mult_genotype_1},{mult_genotype_2})
    that klue was run on
    """
    # restrict to the plates in our input sample set
    temp = df.merge(sample_df, on='plate', how='inner')

    # restrict to this subpool / plate
    temp = temp.loc[temp.plate==wc.plate]
    temp = temp.loc[temp.subpool==wc.subpool]

    # restrict to multiplexed wells
    temp = temp.loc[temp.mult_genotype.notnull()]

    files = expand(cfg_entry,
                   zip,
                   plate=temp.plate.tolist(),
                   subpool=temp.subpool.tolist(),
                   mult_genotype_1=temp.mult_genotype_1.tolist(),
                   mult_genotype_2=temp.mult_genotype_2.tolist())
    files = list(set(files))

    return files

# get counts for each cell across all genotypes
rule klue_get_genotype_counts:
    input:
        adatas = lambda wc:get_subpool_adatas(df, sample_df, wc, config['klue']['unfilt_adata'])
    resources:
        mem_gb = 128,
        threads = 2
    output:
        ofile = config['klue']['genotype_counts']
    run:
        get_genotype_counts(input.adatas, output.ofile)
