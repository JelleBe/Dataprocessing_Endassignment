SAMPLES = ['WT', 'KO1', 'KO2', 'KO3']
rule all:
    input: expand('../../vcf_input/{sample}.vcf', sample=SAMPLES)

rule call_variants:
    input:
        referenceFasta = '../../reference_genome/mm10.fa',
        bams = '../../sorted_bam/{sample}.bam'
    output: '../../vcf_input/{sample}.raw.bcf'
    shell:
        'bcftools mpileup -f {input.referenceFasta} {input.bams} | bcftools call -mv -Ou -o {output}'

rule bcfToVcf:
    input:'../../vcf_input/{sample}.raw.bcf'
    output:'../../vcf_input/{sample}.vcf'
    shell:
         'bcftools view {input} > {output}'
