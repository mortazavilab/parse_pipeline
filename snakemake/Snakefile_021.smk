import pandas as pd

p = os.path.dirname(os.getcwd())
sys.path.append(p)

from utils import *

######## Only need to edit this part ########
config_tsv = 'configs/igvf_021_config.tsv'
sample_csv = 'configs/sample_metadata.csv'

kit = 'WT_mega'  # either WT (48 wells), WT_mini (12 wells), or WT_mega (96 wells)
chemistry = 'v2'  # v1, v2, and v3 are supported

######## Do not change anything past this point ########

configfile: 'configs/config.yml'

# read in config / analysis spec
df = parse_config(config_tsv)
bc_df = get_bc1_matches(kit, chemistry)
sample_df = parse_sample_df(sample_csv)

r1_RT_filename = f"ref/r1_RT_replace_{kit}_{chemistry}.txt"
create_r1_RT_replace(kit, chemistry, r1_RT_filename)

r2r2r3_filename = f"ref/r1r2r3_{kit}_{chemistry}.txt"
create_r1r2r3_file(kit, chemistry, r2r2r3_filename)

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


rule all:
    input:
        expand(config['ref']['klue']['ind'],
               zip,
               mult_genotype_1=[g for g in mult_genotype_1s if g in get_founder_genotypes()],
               mult_genotype_2=[g for g in mult_genotype_2s if g in get_founder_genotypes()]),
        expand(config['klue']['genotype_counts'],
                zip,
                plate=df.plate.tolist(),
                subpool=df.subpool.tolist()),
        expand(config['tissue']['adata'],
                plate=df.plate.tolist(),
                tissue=get_subset_tissues(df, sample_df)),
        expand(config['cellbender']['metrics_copy'],
               plate=df.plate.tolist(),
               subpool=df.subpool.tolist())

################################################################################
########################## Ref download and generation #########################
################################################################################

rule curl_fa:
    resources:
        mem_gb = 4,
        threads = 1,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '4:00:00'
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
        mem_gb =16,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '1:00:00'
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
        mem_gb = 16,
        threads = 1,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '4:00:00'
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
        threads = 12,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '4:00:00'
    conda:
        'envs/kb_env.yaml' 
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
################################ Symlink fastqs ################################
################################################################################

rule symlink_fastq_r1:
    params:
        fastq = lambda wc:get_df_info(wc, df, 'fastq')
    resources:
        mem_gb = 4,
        threads = 1,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '1:00:00'
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
        threads = 1,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '1:00:00'
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
        c1 = config['ref']['kallisto']['c1'],
        c2 = config['ref']['kallisto']['c2'],
        fastq_str = lambda wc:get_subpool_fastqs(wc, df, config, how='str'),
        odir = config['kallisto']['cgb'].split('counts_unfiltered_modified/')[0]
    conda:
        'envs/kb_env.yaml'    
    resources:
        mem_gb = 64,
        threads = 12,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '4:00:00'
    output:
        config['kallisto']['cgb'],
        config['kallisto']['cggn'],
        config['kallisto']['cgg'],
        config['kallisto']['mtx'],
        config['kallisto']['run_info'],
        config['kallisto']['bc_info'],
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
            -r {r1_RT_filename} \
            -w {r2r2r3_filename} \
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
        mem_gb = 64,
        threads = 4,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '1:00:00'
    input:
        mtx = config['kallisto']['mtx'],
        cgb = config['kallisto']['cgb'],
        cggn = config['kallisto']['cggn'],
        cgg = config['kallisto']['cgg'],
        run_info = config['kallisto']['run_info'],
        bc_info = config['kallisto']['bc_info']
    output:
        unfilt_adata = config['kallisto']['unfilt_adata'],
        knee_df = config['report']['knee_df'],
        report_df = config['report']['report_df'],
        cb_df = config['report']['cb_df']
    run:
        make_adata_from_kallisto(input.mtx,
                                 input.cgb,
                                 input.cggn,
                                 input.cgg,
                                 bc_df,
                                 kit,
                                 chemistry,
                                 input.run_info,
                                 input.bc_info,
                                 config_tsv,
                                 wildcards,
                                 output.unfilt_adata,
                                 output.knee_df,
                                 output.report_df,
                                 output.cb_df)


################################################################################
##################################### klue #####################################
################################################################################

rule klue_fa:
    input:
        fa_g1 = lambda wc: expand(config['ref']['genome']['fa'],
                    genotype=wc.mult_genotype_1)[0],
        fa_g2 = lambda wc: expand(config['ref']['genome']['fa'],
                    genotype=wc.mult_genotype_2)[0]
    params:
        temp = config['ref']['klue']['temp']
    conda:
        'envs/kb_env.yaml'  
    resources:
        threads = 32,
        mem_gb = 64,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '4:00:00'
    output:
        fa = config['ref']['klue']['fa'],
        t2g = config['ref']['klue']['t2g']
    shell:
        """
        klue distinguish \
            -o {output.fa} \
            -g {output.t2g} \
            -t 24 \
            -r 61 \
            --all-but-one \
            --tmp {params.temp}/ \
            {input.fa_g1} {input.fa_g2}
        """

rule klue_ind:
    input:
        fa = config['ref']['klue']['fa']
    params:
        temp = config['ref']['klue']['temp']
    conda:
        'envs/kb_env.yaml'  
    resources:
        threads = 32,
        mem_gb = 64,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '4:00:00'
    output:
        ind = config['ref']['klue']['ind']
    shell:
        """
        kb ref \
            --workflow=custom \
            --distinguish \
            -i {output.ind} \
            -t {resources.threads} \
            --tmp {params.temp}/ \
            {input.fa}
        """

rule klue:
    input:
        r1_fastq = lambda wc:get_subpool_fastqs(wc, df, config, how='list', read='R1'),
        r2_fastq = lambda wc:get_subpool_fastqs(wc, df, config, how='list', read='R2'),
        t2g = config['ref']['klue']['t2g'],
        ind = config['ref']['klue']['ind']
    params:
        fastq_str = lambda wc:get_subpool_fastqs(wc, df, config, how='str'),
        odir = config['klue']['cgb'].split('counts_unfiltered_modified/')[0]
    resources:
        mem_gb = 64,
        threads = 24,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '4:00:00'
    conda:
        'envs/kb_env.yaml'  
    output:
        config['klue']['cgb'],
        config['klue']['cggn'],
        config['klue']['cgg'],
        config['klue']['adata'],
        config['klue']['mtx'],
        temporary(config['klue']['bus']),
        temporary(config['klue']['bus_modified_unfilt']),
        temporary(config['klue']['bus_unfilt'])
    shell:
        """
        kb count \
            --h5ad \
            --gene-names \
            --strand=forward \
            --mm \
            -r {r1_RT_filename} \
            -w {r2r2r3_filename} \
            --workflow=standard \
            -g {input.t2g} \
            -x SPLIT-SEQ \
            -i {input.ind} \
            -t {resources.threads} \
            -o {params.odir} \
            --tmp {params.odir}/temp/ \
            --verbose \
            {params.fastq_str}
        """

# add metadata
rule make_adata_klue:
    resources:
        mem_gb = 64,
        threads = 4,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '1:00:00'
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

# get counts for each cell across all genotypes
rule klue_get_genotype_counts:
    input:
        adatas = lambda wc:get_subpool_adatas(df, sample_df, wc, config['klue']['unfilt_adata'])
    resources:
        mem_gb = 64,
        threads = 2,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '1:00:00'
    output:
        ofile = config['klue']['genotype_counts']
    run:
        get_genotype_counts(input.adatas, output.ofile)


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
        mem_gb = 128,
        threads = 8,
        partition = 'gpu',
        account = 'seyedam_lab_gpu',
        gres = 'gpu:1',
        time = '8:00:00'
    output:
        filt_h5 = config['cellbender']['filt_h5'],
        unfilt_h5 = config['cellbender']['unfilt_h5'],
        metrics = config['cellbender']['metrics'],
    shell:
        """
        mkdir -p $(dirname {output.filt_h5})
        cd $(dirname {output.filt_h5})
        
        source ~/miniconda3/bin/activate cellbender

        # Conditionally run the command based on the value of wildcards.subpool
        if [[ "{wildcards.plate}" == "igvf_b01" || "{wildcards.plate}" == "igvf_003" ]]; then
            cellbender remove-background \
                --input ../../../{input.unfilt_adata} \
                --output ../../../{output.unfilt_h5} \
                --total-droplets-included {params.total_drops} \
                --learning-rate {params.learning_rate} \
                --cuda
        else
            cellbender remove-background \
                --input ../../../{input.unfilt_adata} \
                --output ../../../{output.unfilt_h5} \
                --total-droplets-included {params.total_drops} \
                --learning-rate {params.learning_rate} \
                --expected-cells {params.expected_cells} \
                --cuda
        fi
        """

rule copy_cellbender_metrics:
    input:
        metrics = config['cellbender']['metrics']
    resources:
        mem_gb = 4,
        threads = 1,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '1:00:00'
    output:
        metrics_copy = config['cellbender']['metrics_copy']
    shell:
        """
        mkdir -p $(dirname {output.metrics_copy})

        cp {input.metrics} {output.metrics_copy}
        """

#####################################################################################################
##################### Merge klue results with cellbender adata and run scrublet #####################
#####################################################################################################

rule make_filt_adata:
    resources:
        mem_gb = 64,
        threads = 4,
        partition = 'standard',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '1:00:00'
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

rule make_tissue_adata:
    input:
        adatas = lambda wc:get_tissue_adatas(df, sample_df, wc, config['cellbender']['filt_adata'])
    resources:
        mem_gb = 256,
        threads = 2,
        partition = 'highmem',
        account = 'seyedam_lab',
        gres = 'gpu:0',
        time = '2:00:00'
    output:
        adata = config['tissue']['adata']
    run:
        concat_adatas(input.adatas, output.adata)
