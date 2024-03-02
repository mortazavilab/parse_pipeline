rule klue_fa:
    input:
        fa_g1 = lambda wc: expand(config['ref']['genome']['fa'],
                    genotype=wc.mult_genotype_1)[0],
        fa_g2 = lambda wc: expand(config['ref']['genome']['fa'],
                    genotype=wc.mult_genotype_2)[0]
    params:
        temp = config['ref']['klue']['temp']
    resources:
        threads = 32,
        mem_gb = 128
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
    resources:
        threads = 32,
        mem_gb = 128
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
