SAMPLES = ['WT', 'KO1', 'KO2', 'KO3']
rule all:
    input: '../../vcf_input/merged_KO.qd.vcf'

rule call_variants:
    input:
        referenceFasta = '../../reference_genome/mm10.fa',
        bams = '../../sorted_bam/{sample}.bam'
    output: '../../vcf_input/{sample}.raw.bcf'
    shell:
        'bcftools mpileup -f {input.referenceFasta} {input.bams} | bcftools call -mv -Ou -o {output}'

rule bcf_to_vcf:
    input:'../../vcf_input/{sample}.raw.bcf'
    output:'../../vcf_input/{sample}.vcf'
    shell:
         'bcftools view {input} > {output}'
rule filter_quality_of_vcf:
    input: '../../vcf_input/{sample}.vcf'
    output: '../../vcf_input/{sample}.q.vcf'
    shell: 'python3 scripts/vcf_filter.py quality --infile {input} -q 30 --outfile {output}'

rule filter_depth_of_vcf:
    input: '../../vcf_input/{sample}.q.vcf'
    output: '../../vcf_input/{sample}.q.d.vcf'
    shell: 'python3 scripts/vcf_filter.py depth --infile {input} -d 10 --outfile {output}'

rule zip_vcf:
    input: '../../vcf_input/{sample}.q.d.vcf'
    output: '../../vcf_input/{sample}.q.d.vcf.gz'
    shell: 'bgzip -i {input}'

rule index_vcf_file:
    input: '../../vcf_input/{sample}.q.d.vcf.gz'
    output: '../../vcf_input/{sample}.q.d.vcf.gz.tbi'
    shell: 'tabix {input}'

rule merge_knockout_files:
    input:
        ko1 = '../../vcf_input/KO1.q.d.vcf.gz',
        ko2 = '../../vcf_input/KO2.q.d.vcf.gz',
        ko3 = '../../vcf_input/KO3.q.d.vcf.gz',
        ko1_indexed = '../../vcf_input/KO1.q.d.vcf.gz.tbi',
        ko2_indexed = '../../vcf_input/KO2.q.d.vcf.gz.tbi',
        ko3_indexed = '../../vcf_input/KO3.q.d.vcf.gz.tbi'
    output: '../../vcf_input/merged_KO.qd.vcf'
    shell: 'bcftools merge --force-samples {input.ko1} {input.ko2} {input.ko3} > {output}'
