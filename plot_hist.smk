configfile: "config.yaml"
print(config)

rule plot_histogram:
    input:
        wild_type = '../../vcf_input/WT.q.d.vcf',
        merged_KO = '../../vcf_input/merged_KO.qd.vcf'
    output: 'results/histograms.pdf'
    shell: 'Rscript scripts/plot.R {input.wild_type} {input.merged_KO} {output}'

rule all:
    input: 'results/histograms.pdf'
