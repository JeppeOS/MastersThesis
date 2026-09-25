Scripts used for visualizing result and decisions made in regard to Globdb-data.

01: Plot jaccard sensitivity
01_plot_jaccard_sensitivity.R: Visualizes the sensitivity of the protein-family co-occurrence network to alternative Jaccard similarity thresholds, including changes in retained nodes, edges, connected components, and largest-component size. Supports selection of the final Jaccard threshold of 0.25. (Figure X in thesis).

02: MCL-inflation
02_plot_mcl_inflation_sensitivity: Visualizes the sensitivity of MCL community detection to alternative inflation values. Compares community number and assignment stability relative to the selected inflation value of 2.0. (Figure X in thesis)

03: Module occurrence vs. completeness
03_plot_module_occurrence_and_completeness.R: Visualizes the distribution and median completeness of selected MCL modules across the Methylococcales genome collection. (Figure X in thesis)

04: Overview figure
04_plot_comparative_genomics_pipeline_vertical.R: Generates a schematic overview of the comparative-genomics workflow, from genome annotation and exported-heme protein identification through MMseqs2 clustering, co-occurrence analysis, MCL clustering, and genomic-context analysis.

05: Module06 gene map
05.1_prepare_module06_gene_map_inputs.py: Builds the genomic-context dataset for Module 06 and its focal Cluster_00063 MCA0421-like proteins, combining protein membership, coordinates, neighborhood information, and relevant annotations.

05.2_select_module06_representatives.py: Selects representative Module-06/Cluster_00063 loci for visualization using taxonomy, module occurrence, genomic-context observability, and neighborhood composition. Produces the final selected gene-map dataset.

06: Module10 gene map
06.1_audit_cluster00050_module_architectures_v2.py: Audits recurrent genomic architectures surrounding Cluster_00050 Cyc2 proteins in Module 10 and resolves duplicated focal observations before representative-locus selection.

06.2_prepare_module10_gene_map_data.py: Constructs gene-resolved genomic neighborhoods centred on Module-10 Cluster_00050 proteins for subsequent local-family and architecture analysis.

06-7.3-8.2_run_local_architecture_pipeline.py: Generalized local-neighborhood analysis used for Modules 10, 20, and 35. Locally clusters neighborhood proteins with MMseqs2, assigns module-specific local-family IDs, and summarizes family prevalence, gene order, spacing, orientation, and recurrent genomic architectures.

06.4_inspect_module10_local_architectures.py: Summarizes recurrent local protein families and neighborhood architectures around Module-10/Cyc2 loci to guide interpretation and representative selection.

06.5_prepare_module10_selected_gene_maps.py: Selects representative Module-10 architectures and produces the final oriented and annotated gene-level table used for the Module-10 gene-map figure.

07: Module20 gene map
07.1_inspect_module20_inputs.py: Audits the Module-20/Cluster_00035 input data, including focal-protein assignments, duplicate observations, annotations, and neighborhood information before architecture analysis.

07.2_prepare_module20_gene_map_data.py: Constructs gene-resolved genomic neighborhoods surrounding Module-20 Cluster_00035 MtoA/MtrA-family proteins.

07.4_build_module20_architecture_matrix.py: Converts Module-20 local-family occurrences into region-level architecture profiles and identifies recurrent combinations of focal proteins, neighboring families, and accessory genes.

07.5_prepare_module20_selected_gene_maps.py: Selects representative Module-20 genomic architectures and generates the final oriented and annotated gene-map dataset used for visualization.

08: Module 35 gene map
08.1_prepare_module35_gene_map_data.py: Constructs the gene-resolved neighborhood dataset for Module 35 and its focal Cyc2-family proteins prior to local-family and architecture analysis.

09: Find contig boundaries
09_check_selected_contig_boundaries.py: Checks the amount of assembled sequence available around all selected Module-06, -10, -20, and -35 loci. Produces cross-figure QC used to distinguish true contig boundaries from the chosen plotting window.

10: Plot gene maps
10.1_plot_module06_selected_gene_maps_streamlined.R: Generates the final streamlined gene-map figure for the selected Module-06/Cluster_00063 neighborhoods. (Figure X)

10.2_plot_module10_selected_gene_maps_streamlined.R: Generates the final streamlined gene-map figure for representative Module-10/Cluster_00050 Cyc2 neighborhoods.

10.3_plot_module20_selected_gene_maps_streamlined.R: Generates the final streamlined gene-map figure for representative Module-20/Cluster_00035 MtoA/MtrA-family neighborhoods.

10.4_plot_module35_gene_maps_streamlined.R: Generates the final streamlined gene-map figure for Module-35 Cyc2-centered neighborhoods, including recurrent and lineage-specific local genomic architectures.











