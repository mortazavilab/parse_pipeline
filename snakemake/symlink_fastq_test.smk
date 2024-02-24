rule symlink_fastq:
    input:
        "/share/crsp/lab/seyedam/share/igvf_splitseq/igvf_015/fastq/015_67C_S1_L001_R1_001.fastq.gz"
    output:
        "fastq/015_67C_S1_L001_R1_001.fastq.gz"
    shell:
        "ln -s {input} {output}"