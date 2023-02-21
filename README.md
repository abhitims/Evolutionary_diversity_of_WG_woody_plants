# Evolutionary diveristy of Western Ghats woody plants

The code examines the evolutionary diversity of woody plants in the Western Ghats. It takes in the output of the SDMs run for 348 species and occurence locations for 122 species, which is stacked at 10x10 km resolution. The phylogenetic tree is built using V. phylomaker. For more details refer to the Methods section in the manuscript. The code for the SDMs is taken from https://github.com/bhartidk/centipede_diversity_endemism

# GENERAL INFORMATION

**Title: Range restricted old and young lineages show the southern Western Ghats to be both museum and cradle of diversity for woody plants**  

Authors: Abhishek Gopal, D. K. Bharti,  View Navendu Page,  Kyle G. Dexter, Ramanathan Krishnamani, Ajith Kumar, Jahnavi Joshi

**Abstract**:
The Western Ghats (WG) mountain chain is a global biodiversity hotspot with high diversity and endemicity of woody plants. The latitudinal breadth of the WG offers an opportunity to determine the evolutionary drivers of latitudinal diversity patterns. We examined the spatial patterns of evolutionary diversity using complementary phylogenetic diversity and endemism measures. To examine if different regions of the WG serve as a museum or cradle of evolutionary diversity, we examined the distribution of 470 species based on distribution modelling and occurrence locations across the entire region. In accordance with the expectation, we found that the southern WG is both a museum and cradle of woody plant evolutionary diversity, as both old and young evolutionary lineages are restricted to the southern WG. The diversity gradient is likely driven by high geo-climatic stability in the south and phylogenetic niche conservatism for moist and aseasonal sites. This is corroborated by persistent lineage nestedness at almost all evolutionary depths (10–135 million years), and a strong correlation of evolutionary diversity with drought seasonality, precipitation and elevation. Our results highlight the global value of the WG, demonstrating, in particular, the importance of protecting the southern WG – an engine of plant diversification and persistence.

Keywords: Cradles; Museums; Nestedness; Phylogenetic diversity; Phylogenetic endemism; Western Ghats 


# SHARING/ACCESS INFORMATION

The data and the code can be used for non-commercial purposes (The MIT Licence).  

**Recommended citation for this dataset:**  
Gopal, A., Bharti, D. K., Page, N., Dexter, K. G., Krishnamani, R., Kumar, A., & Joshi, J. (2023). Range restricted old and young lineages show the southern Western Ghats to be both museum and cradle of diversity for woody plants. *Proceedings of the Royal Society B*.

# DATA OVERVIEW

**CSV files**

1. biodiv_final_all_470sp_6Dec22:  
2. biodiv_final_all_subset_470sp_6Dec22:  
3. biodiv.mat.clip.SDM_only_348_6Dec22:  
4. All_sp_checklist_node_final_6Dec22:  
5. WG_pred_470_CWD:  

**Tree file as created by V. phylomaker**  
1. Tree_470_6Dec22:  Phylogenetic tree created using V.phylomaker


**Shapefiles**  

1. ne_10m_coastline:  
2. WG_boundary:  


# Code overview   

1. Loading all packages and reading and cleaning the data    
  
     1.a. Loading all the packages    
     1.b. Input shape files of WG    
     1.c. Read input data regarding occ location    
     1.d. Read input tree    
     1.e. Data cleaning for making all the names compatible  
     1.f. Custom theme for ggplot  
     
2. Calculating all diversity indices  
 
     2.a. Measuring TILD  
     2.b. Measuring SR and PD  
     2.c. Measuring PE  
     2.d. Plotting the evolutionary indices  
     2.e. Diversity measure wrt Biog zone  
     2.f. Correlations between diversity indices  
     
3. PD at diff time period, LTT at each lat and lineage age trends  
 
     3.a. PD at diff time period  
     3.a.1. Tree slices  
     3.b. LTT plots for each lat bin  
     3.c. lineage age trends  
     3.d. lineage turnover  
     
4. Order family genus contribution to PD  
     
     4.a. Contribution of Superorder to PD at each time period  
     4.b. Key order family genus contributing to PD per biog  
     4.c. Key families contributing to PD per lat bins  
     4.d. Summarizing Family level and species information     
     
5. Examining the correlates of evolutionary diveristy  
     
     5.a. Examining correlations among WorldClim the predictors  
     5.b. Key predictor distribution and relationship with lat  
     5.c. Correlates of evolutionary diveristy  
     
6. Plotting all the figures in the manuscript and the supplementary  
 

