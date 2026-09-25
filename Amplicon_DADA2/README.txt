16S amplicon analysis:
make_trainset_final.R: Takes bac120_ssu_reps_r226.fna and ar53_ssu_reps_r226.fna, which contains the data for labeling amplicons, and turns them into train set capable of being used by the DADA2 pipeline.

run_Cutadapt.bat: Runs cutadapt on the raw amplicon reads to remove adapters and primers.

DADA2_amplicon_analysis.Rmd: Processes and visualizes the amplicon data post-cutadapt. Credit goes to Justus Nwese for the original script.
