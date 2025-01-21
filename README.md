# Asymmetric gene flow maintains range edges in a marine invertebrate

ACCESS INFORMATION
1. Licenses/restrictions: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication
2. Recommended citation for this data/code archive: https://github.com/zachbhancock/asym_gene_flow_proj

DATA & CODE FILE OVERVIEW
This data repository consist of 12 data files, 13 code scripts, 9 figure files, and this README document, with the following data and code filenames and variables.

Data files and variables in data folder
1. Filtered VCF file (feems_filt.recode.vcf)
2. BED files (feems_bed.bed, feems_bed.bim, feems_bed.bed)
3. FST and geographic distance data frame (fst_dist_data.txt)
4. Location files (haust_locations.coord, haust_locations.outer)
5. List of individual locations (new_indv_file.txt)
6. Population locations (new_pop_file.txt)
7. Empirical psi results (psi_results)
8. Summary statistics (shovel_bugs_new_results.txt)
9. Simulated psi results (sim_psi_results.txt)

Code scripts and workflow in script folder
1. Statistical analysis in R (analysis_plots.R)
2. SLiM simulation recipes (center_asym.slim, center_sym.slim, equal_asym.slim, equal_sym.slim)
3. Perform the feems analysis in Python (feems_script.py)
4. Bash scripts for de novo assembly and filtering (de_novo_map.sh, populations_stacks.sh)
5. Custom stacks script for handling 2RAD data, provided by J. Catchen (stacks-2.67.tar.gz)
6. Parse the simulation psi results in Python (psi_sim_parse.py)
7. Prepare data and perform conStruct analysis (prep_amphipod_conStruct.R, run_amphipod_conStruct.R)
8. Filter data for pixy (filtering_pixy.sh)

SOFTWARE VERSIONS
1. STACKS v.2.67
2. vcftools v.1.16
3. pixy v.1.2.7.beta1
4. VCF2PCACluster v.1.41
5. fastSTRUCTURE v.1.0
6. conStruct v.1.0.3
7. triangulaR v.0.0.0.9000
8. feems v.1.0.0
9. rangeExpansion (no official release version)
10. SLiM v.4.0.1
11. pyslim v.1.0.1
12. tskit v.0.5.3
13. msprime v.1.2.0
14. plink v.1.9
15. R v.4.4.2
