rule call_variants:
    input:
        reference_fasta = 'reference_genome/mm10.fa',
        bam = 'sorted_bam/{sample}.bam'
    output: 'vcf_input/{sample}.raw.bcf'
    log: 'log_files/{sample}.raw.bcf.log'
    threads: 8
    message: "Calling variants of {input.bam} into output file {output} using reference genome {input.reference_fasta} with {threads} threads"
    shell:
        'bcftools mpileup --threads {threads} -f {input.reference_fasta} {input.bam} | bcftools call -mv -Ou -o {output}'

rule bcf_to_vcf:
    input:'vcf_input/{sample}.raw.bcf'
    output:'vcf_input/{sample}.vcf'
    log: 'log_files/{sample}.vcf.log'
    message: "Converting raw binary variant call file {input} to variant call file {output}"
    shell:
         'bcftools view {input} > {output}'

rule filter_quality_of_vcf:
    input: 'vcf_input/{sample}.vcf'
    output: 'vcf_input/{sample}.q.vcf'
    log: 'log_files/{sample}.q.vcf.log'
    message: "Filtering variant call file {input} for quality with threshold being 30 into quality filtered variant call file {output}"
    shell: 'python3 scripts/vcf_filter.py quality --infile {input} -q 30 --outfile {output}'

rule filter_depth_of_vcf:
    input: 'vcf_input/{sample}.q.vcf'
    output: 'vcf_input/{sample}.q.d.vcf'
    log: 'log_files/{sample}.q.d.vcf.log'
    message: "Filtering quality filtered variant call file {input} for depth with threshold being 10 into quality filtered, depth filtere variant call file {output}"
    shell: 'python3 scripts/vcf_filter.py depth --infile {input} -d 10 --outfile {output}'

rule zip_vcf:
    input: expand('vcf_input/{sample}.q.d.vcf', sample=config["samples"])
    output: 'vcf_input/{sample}.q.d.vcf.gz'
    log: 'log_files/{sample}.q.d.vcf.gz.log'
    message: "Zipping quality and depth filtered variant call file {input} to zipped {output}"
    shell: 'bgzip -i {input} > {output}'

rule index_vcf_file:
    input: 'vcf_input/{sample}.q.d.vcf.gz'
    output: 'vcf_input/{sample}.q.d.vcf.gz.tbi'
    log: 'log_files/{sample}.q.d.vcf.gz.tbi.log'
    message: "Indexing zipped, quality and depth filtered variant call file {input} to zipped and indexed {output}"
    shell: 'tabix {input}'

rule merge_knockout_files:
    input:
        ko1 = 'vcf_input/KO1.q.d.vcf.gz',
        ko2 = 'vcf_input/KO2.q.d.vcf.gz',
        ko3 = 'vcf_input/KO3.q.d.vcf.gz',
        ko1_indexed = 'vcf_input/KO1.q.d.vcf.gz.tbi',
        ko2_indexed = 'vcf_input/KO2.q.d.vcf.gz.tbi',
        ko3_indexed = 'vcf_input/KO3.q.d.vcf.gz.tbi'
    output: 'vcf_input/merged_KO.qd.vcf'
    log: 'log_files/merged_KO.log'
    message: "Merging variant call files of all knockouts into {output}"
    threads: 8
    shell: 'bcftools merge --threads {threads} --force-samples {input.ko1} {input.ko2} {input.ko3} > {output}'
