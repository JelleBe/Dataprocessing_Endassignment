"""
Snakefile to sort all .bam files
"""
SAMPLES = ['WT', 'KO1', 'KO2', 'KO3']

rule sort_bam_files:
    input: '../../unsorted_bam/{sample}.bam'
    output: '../../sorted_bam/{sample}.bam'
    shell:
        'samtools sort {input} > {output}'

rule all:
    input: expand('../../sorted_bam/{sample}.bam', sample=SAMPLES)