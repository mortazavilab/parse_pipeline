rule klue_ind:
    


rule klue:
    input:
        r1_fastq = lambda wc:get_subpool_fastqs(wc, df, config, how='list', read='R1'),
        r2_fastq = lambda wc:get_subpool_fastqs(wc, df, config, how='list', read='R2'),
        t2g = config['ref']['klue']['t2g'],
        ind = config['ref']['klue']['ind']
    # conda:
    #     "hpc3sc"
    params:
        # TODO bc1 map, barcodes, c1, c2 should be output from sth, seqspec
        bc1_map = config['ref']['bc1_map'],
        barcodes = config['ref']['barcodes'],
        fastq_str = lambda wc:get_subpool_fastqs(wc, df, config, how='str'),
        odir = config['klue']['cgb'].split('counts_unfiltered_modified/')[0]
    resources:
        mem_gb = 64,
        threads = 24
    output:
        config['klue']['cgb'],
        config['klue']['cggn'],
        config['klue']['cgg'],
        config['klue']['adata'],
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
            -r {params.bc1_map} \
        	-w {params.barcodes} \
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
