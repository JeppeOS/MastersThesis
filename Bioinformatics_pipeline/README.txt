Comparative genomics pipeline used for Jeppes masters thesis

00: Genome extraction.
00_extract_methylococcales_from_globdb.sh: Extracts the completete Methylococcales genome set from the globdb r226 dataset. 

01: FeGenie annotation
01_run_fegenie_methylococcales.sbatch: Runs the FeGenie tool on the extracted Methylococcales genomes. Prodigal provides .faa outputs in addition to FeGenie's predictions.

02: Signal peptide prediction
02_run_signalp6_array.sbatch: Runs SignalP6 on the FeGenie derived proteomes.

03: CXXCH-prediction
03_run_findmehemes.sbatch: Runs FindMeHemes on the derived proteomes.

04: Protein topology analysis
04.1_prepare_deeptmhmm_candidates.py: Constructs the protein set used as DeepTMHMM input. These are drawn from FeGenie- and/or FindMeHemes-positive proteins.

04.2_run_deeptmhmm_candidates.sbatch: Runs DeepTMHMM on the candidate proteins.

05: Integrated protein annotation
05_build_integrated_annotations.py: Combines the output from FeGenie, SignalP6, FindMeHemes, and DeepTMHMM. The constructed table keep annotations for downsteam analysis.

06: Protein-family clustering
06.1_prepare_clustering_input.py: Selects proteins with heme-motifs and evidence for export for clustering.

06.2_run_mmseqs2_clustering.sbatch: Runs MMseqs2 on the selected proteins with parameters 40% minimum sequence identity and 80% minimum coverage.

06.3_build_cluster_tables.py: Converts raw clustering output to tables for downsteam use.

07: Protein-family co-occurence network
07.1_build_cooccurrence_diagnostics.py: Evaluates the distribution and co-occurrence of MMseqs2 protein families across the 631 genomes. Clusters present in fewer than five genomes are excluded from the main network analysis. Pairwise genome-distribution similarity is calculated using the Jaccard similarity. 

07.2_build_final_cooccurrence_network.py: Constructs protein family co-occurence network for Cytoscape. In cytoscape community detection was performed using Markov Clustering (MCL) though the clusterMaker2 cytoscape plugin. MCL-transformed modules are exported.

07.3_build_module_membership.py: Reads the MCL-transformed modules for future use.

08: Protein family protein characterization
08.1_build_fegenie_cluster_composition.py: Summarizes FeGenie annotation composition within the MMseqs2 clusters.

08.2_analyze_cluster_00035_sequence_coherence.py: Performs detailed analysis of cluster_00035.

08.3_summarize_focal_cluster_sizes_lengths.py: Produces deltailed analysis of all other clusters of interest.

08.4_build_cluster00063_MCA0421_tables.py: Produces deltailed analysis of cluster_00063.

08.5_screen_methylobacter_smmo.py: Screens Methylobacter-genomes with cluster_00063 for soluble methane monooxygenase.

09: Module occurence
09_build_module_occurrence.py: Combines MMseq2 clusters, MCL module assignment, and genome membership.

10: Genomic coordinates
10_build_genomic_coordinates.py: Maps module-associated proteins back to their genomic coordinates.

11: Annotated gene catalog
11_build_annotated_gene_catalog.py: Constructs a unified genomic gene catalog for the Methylococcales genomes. The gene catalog integrates genomic coordinates with available functional annotations and protein-family assignments. Annotations comes from FeGenie, COG, and corresponding Globdb .gff files.

12: Neighborhood extraction
12.1_build_context_observability.py: Determines how much genomic context is observable around each protein-family protein.

12.2_extract_observed_neighborhoods.py: Extracts genomic information based on 12.1.

12.3_extract_cluster00063_context.py: Extracts the Stage-12 genomic-context information specifically associated with Cluster_00063. This was done at a later point hence not a part of 12.2.

13: Neighborhood conservation (tries to build a consensus neighborhood around genes of interest. This ended up being exploratory however as the output is better at individual genome-gene-level.)
13A_summarize_module_colocalization.py: Tests if protein families in the same module are also physically co-located within member genomes.

13B.1_summarize_mmseqs_neighbor_conservation.py: Measures reccurence of MMseq2 protein families around the focal protein family.

13B.2_summarize_mmseqs_gene_window_conservation.py: Performs a complementary neighbor-conservation analysis using gene-count windows rather than only physical base-pair distances.

13C.1_summarize_cog_context_fast.py: Examines conserved functional annotation around focal proteins using COG assignments.

13C.2_summarize_intervening_gene_context.py: Examines the actual genes located between pairs of module-associated proteins.

13C.3_summarize_local_intervening_gene_context.py: Summarizes recurrent local arrangements and identifies conserved gene patterns occurring between physically associated module proteins.

13D_build_consensus_neighborhoods.py: Integrates the Stage-13 neighborhood analyses into consensus neighborhood for the diffrent modules.

14: Module report (This stage collects all the results so far for each module. These are also only explorary.)
14.1_build_module_report_data.py: Collects relevant data.

14.2_build_module_architecture.py: Uses the report data to describe recurrent genomic architectures associated with each module based on 13.

14.3_render_module_reports.py: Renders readable reports.

14.4_build_refined_module_reports.py: Produces a final, refined module reports.

15: Gene-level neighborhood analysis
15.1_build_gene_level_resolution_recovery.py: Reconstructs gene-level neighborhood information for selected focal modules while retaining the full set of actual genes surrounding focal proteins.

15.2_build_gene_level_associations.py: Summarizes recurring gene-level associations in the reconstructed neighborhoods.




