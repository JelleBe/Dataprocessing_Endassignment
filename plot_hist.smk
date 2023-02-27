configfile: "config.yaml"

rule plot_histogram:
    input:
        wild_type = 'vcf_input/WT.q.d.vcf',
        merged_KO = 'vcf_input/merged_KO.qd.vcf'
    output: 'results/histograms.pdf'
    log: 'log_files/histograms.log'
    script: 'scripts/plot.R'