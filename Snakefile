"""
Main script.
Includes all the Snakefiles necessary
for executing the pipeline.
"""
configfile: "config.yaml"

include: "sort_bams.smk"
include: 'vcf.smk'
include: "plot_hist.smk"

rule all:
    input:
        'results/histograms.pdf'

onsuccess:
    print("Pipeline ran succesfully, no errors occurred.")       

onerror:
    print("An error occurred, check the log files for additional information on what went wrong.")       
