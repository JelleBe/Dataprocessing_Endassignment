"""
Snakefile to sort all .bam files
"""

rule sort_bam_files:
    input: 'unsorted_bam/{wildcards.sample}.bam'
    output: 'sorted_bam/{wildcards.sample}.bam'
    log: 'sorted_bam/{wildcards.sample}.bam.log'
    threads: 8
    message: "Sorting {wildcards.sample}.bam"
    shell:
        'samtools sort --threads {threads} {input} > {output}'
