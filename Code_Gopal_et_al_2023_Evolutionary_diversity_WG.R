# --------------------------------General Information--------------------------------------------------#
#Range restricted old and young lineages show the southern Western Ghats to be both museum and cradle of diversity for woody plants
#Authors: Abhishek Gopal,  D. K. Bharti,  Navendu Page,  Kyle G. Dexter, Ramanathan Krishnamani, Ajith Kumar, Jahnavi Joshi

#Recommended citation for this script:
#Gopal, A., Bharti, D. K., Page, N., Dexter, K. G., Krishnamani, R., Kumar, A., & Joshi, J. (2023). Range restricted old and young lineages show the southern Western Ghats 
#to be both museum and cradle of diversity for woody plants. Proceedings of the Royal Society B.

#Correspondence: abhishekg@csirccmb.org; abhishekgopal1993@gmail.com
# ------------------------------------------------------------------------------------------------------#

#rm(list = ls())

#Collapse all for better readability and to get the structure of the code.
#Edit -> Folding -> Collapse_all


final_sub_path <- paste0("C:/Users/Downloads/Evolutionary_diveristy_of_Western_Ghats_woody_plants-main/Evolutionary_diveristy_of_Western_Ghats_woody_plants-main")

# ----------------------------Structure of the code----------------------------------------------#
# 1.Loading all packages and reading and cleaning the data
#     1.a.Loading all the packages
#     1.b.Input shape files of WG
#     1.c.Read input data regarding occ location
#     1.d.Read input tree
#     1.e.Data cleaning for making all the names compatible
#     1.f.Custom theme for ggplot
# 2.Calculating all diversity indices
#     2.a.Measuring TILD
#     2.b.Measuring SR and PD
#     2.c.Measuring PE
#     2.d.Plotting the evolutionary indices
#     2.e.Diversity measure wrt Biog zone
#     2.f.Correlations between diversity indices
# 3.PD at diff time period, LTT at each lat and lineage age trends
#     3.a. PD at diff time period
#     3.a.1.Tree slices
#     3.b. LTT plots for each lat bin
#     3.c. lineage age trends
#     3.d. lineage turnover
# 4.Order family genus contribution to PD
#     4.a. Contribution of Superorder to PD at each time period
#     4.b. Key order family genus contributing to PD per biog
#     4.c. Key families contributing to PD per lat bins
#     4.d. Summarizing Family level and species information   
# 5.Examining the correlates of evolutionary diveristy
#     5.a.Examining correlations among WorldClim the predictors
#     5.b.Key predictor distribution and relationship with lat
#     5.c.Correlates of evolutionary diveristy
# 6. Plotting all the figures in the manuscript and the supplementary
# -----------------------------------------------------------------------------------------#


# -----------------------------------------------------------------------------------------#
# --------------------1.Loading all packages and reading and cleaning the data-------------#
# -----------------------------------------------------------------------------------------#

#-------------------------1.a.Loading all the packages---------------------------------------
#-------------Run this to install all the packages---------------#
# Code modified from https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them

packages_install_load <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}
#--------------Loading/installing all packages------------------------------#
#Data wrangling
packages_install_load(c("tidyverse","reshape2","stringr", "stringi", "readxl","Taxonstand" ))
#SDM maxent
packages_install_load(c("raster","sp","rgdal","maps","rgeos","dismo","dplyr",
                        "Hmisc","devtools","digest","rJava","geosphere","ncdf4","sf","MASS"))
#Tree plotting, pruning, and wrangling
packages_install_load(c("phylobase","ggtree","V.PhyloMaker"))
#Regression analysis
packages_install_load(c("WRTDStidal","quantreg","lmtest","corplot","jtools"))
#Phylogenetic indices
packages_install_load(c("ape","picante","zoo","phytools","phyloregion"))

#lineage age graphs
packages_install_load(c("dendextend","vegan","abind"))
#Nestedness
packages_install_load("metacom")
#Plotting
packages_install_load(c("paletteer","patchwork","gridExtra","gt","ggrepel","corrplot","cowplot","psych"))#,"phyloch"))

#-----------------------------------------------------------------------------------------#
#-------------------------1.b.Input shape files of WG---------------------------------------

wgs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
wg <- readOGR(dsn=paste0(final_sub_path,'/Appendix_S3_processed_data/processed_data/Shapefile/WG_boundary'), layer='wg_boundary')

#-------------------------1.c.Read input data regarding occ location----------------------------
biodiv_final_all <- read.csv(paste0(final_sub_path,"/Appendix_S3_processed_data/processed_data/biodiv_final_all_470sp_6Dec22.csv"),strip.white=TRUE)
biodiv_final_all.subset <- read.csv(paste0(final_sub_path,"/Appendix_S3_processed_data/processed_data/biodiv_final_all_subset_470sp_6Dec22.csv"),strip.white=TRUE)
biodiv.mat.clip <- read.csv(paste0(final_sub_path,"/Appendix_S3_processed_data/processed_data/biodiv.mat.clip.SDM_only_348_6Dec22.csv"),strip.white=TRUE)

WG_pred <- read.csv(paste0(final_sub_path,"/Appendix_S3_processed_data/processed_data/WG_pred_470_CWD.csv"),strip.white=TRUE)


#-------------------------1.d.Read input tree -----------------------------------------------
phylo_tree_final <- read.tree(paste0(final_sub_path,"/Appendix_S3_processed_data/processed_data/Tree_470_6Dec22.tre"))

# The data modified from the input tree to include order, family name and node number
All_sp_checklist.node <- read.csv(paste0(final_sub_path,"/Appendix_S3_processed_data/processed_data/All_sp_checklist_node_final_6Dec22.csv"),strip.white=TRUE)

#----------------------------------------------------------------------------------------------#
#-------------------------1.e.Data cleaning for making all the names compatible-----------------
#As csv changes - to . while saving/reading
biodiv_final_all.subset
setdiff(colnames(biodiv_final_all.subset)[-c(1,2)],phylo_tree_final$tip.label )
colnames(biodiv_final_all.subset)[69] <- "Chionanthus_mala-elengi"
colnames(biodiv_final_all.subset)[241] <- "Melicope_lunu-ankenda"
colnames(biodiv_final_all.subset)[426] <- "Breynia_vitis-idaea"
colnames(biodiv_final_all.subset)[157] <- "Garcinia_gummi-gutta"


biodiv_final_all
setdiff(colnames(biodiv_final_all)[-c(1,2)],phylo_tree_final$tip.label )
colnames(biodiv_final_all)[68] <- "Chionanthus_linocieroides"
colnames(biodiv_final_all)[69] <- "Chionanthus_mala-elengi"
colnames(biodiv_final_all)[241] <- "Melicope_lunu-ankenda"
colnames(biodiv_final_all)[426] <- "Breynia_vitis-idaea"
colnames(biodiv_final_all)[157] <- "Garcinia_gummi-gutta"


biodiv.mat.clip
setdiff(colnames(biodiv.mat.clip)[-c(1,2)],phylo_tree_final$tip.label )
colnames(biodiv.mat.clip)[69] <- "Chionanthus_mala-elengi"
colnames(biodiv.mat.clip)[68] <- "Chionanthus_linocieroides"
colnames(biodiv.mat.clip)[241] <- "Melicope_lunu-ankenda"
colnames(biodiv.mat.clip)[157] <- "Garcinia_gummi-gutta"

#Creating a matrix for diversity calculations
WG_mat <-  as.matrix(biodiv_final_all.subset[,-c(1,2)])
colnames(WG_mat)
rownames(WG_mat) <- paste0(biodiv_final_all.subset$x,":",biodiv_final_all.subset$y)

#-------------------------1.f.Custom theme for ggplot---------------------------
theme_custom <- 
  theme_classic(base_size = 10)+# labels = c("<0.01", "<0.025", "Not significant", ">0.975",">0.99"))
  theme( axis.title = element_text(size=14, face="bold", family="serif"),
         axis.text = element_text(size=10,colour = "black",family = "serif" ),
         legend.title = element_text(size=10, face = "bold",family = "serif"),
         legend.text = element_text(size=10,family = "serif"))

theme_custom2<- 
  theme_bw(base_size = 15)+# labels = c("<0.01", "<0.025", "Not significant", ">0.975",">0.99"))
  theme( axis.title = element_text(size=14,face="bold", family="serif"),
         axis.text = element_text(size=14,colour = "black",family = "serif" ),
         legend.title = element_text(size=10, face = "bold",family = "serif"),
         legend.text = element_text(size=10,family = "serif"))
############################################################################################


#------------------------------------------------------------------------------------------#
#----------------------2.Calculating all diversity indices---------------------------------#
#------------------------------------------------------------------------------------------#


#-------------------------2.a. Measuring TILD------------------------------------------------
#------------------------------------------------------------------------------#
# This code is modified from Dexter et al. 2019 and Griffiths et al. 2021

#Function to calculate the number of lineages in a grid at specified time slice
#Input is commat: A community matrix in a wide format; phy: a phylogeny; bin_size: time slice

no_lineages <- function(commat,phy,bin_size) {
  
  tmp_seq <- rev(seq(0,max(branching.times(phy),by=bin_size)))
  matrix_out <- matrix(NA,nrow(commat),length(tmp_seq))
  rownames(matrix_out) <- rownames(commat)
  colnames(matrix_out) <- tmp_seq
  
  for (i in 1:nrow(commat)) {
    
    tmp_spp <- colnames(commat)[which(commat[i,]>0)]
    tmp_phy <- drop.tip(phy,which(!phy$tip.label%in%tmp_spp))
    tmp_btimes <- branching.times(tmp_phy)
    
    for (j in 1:ncol(matrix_out)) {
      
      matrix_out[i,j] <- length(which(tmp_btimes > tmp_seq[j]))
      
    }
    
  }
  
  matrix_out <- matrix_out + 1
  matrix_out[,1] <- 1
  
  matrix_out[,ncol(matrix_out)] <- apply(commat,1,sum)
  
  return(matrix_out)
  
}

rownames(biodiv_final_all.subset) <- paste0(biodiv_final_all.subset$x,":",biodiv_final_all.subset$y)

#Calling the function to calculate the number of lineages in each 10 x 10 grid cell
# at a time slice for 1 Mya

no_lineages_sim <- no_lineages( biodiv_final_all.subset,phylo_tree_final,1)

# plot this to make sure it looks reasonable
plot(-as.numeric(colnames(no_lineages_sim)),no_lineages_sim[1,],type="l")
plot(-as.numeric(colnames(no_lineages_sim)),no_lineages_sim[105,],type="l")

# "This function calculates TILD, both 'raw' and log-transformed version
#  lineages = a lineage composition matrix derived from the above function"

#  "In Dexter et al. 2019. Forests, were referred to the log-transformed version
#   as the TILD metric
#   The raw (non-log-tranformed) version is in fact mathematically identical to raw 
#   phylogenetic diversity, and is included here for testing"

TILD <- function(lineages) {
  
  lineages <- lineages[,ncol(lineages):1]
  
  PD <- vector("numeric",nrow(lineages))
  TILD <- vector("numeric",nrow(lineages))
  tmpx <- as.numeric(colnames(lineages))
  
  for (i in 1:nrow(lineages)) {
    
    tmpy <- lineages[i,]
    tmplogy <- log(lineages[i,])
    PD[i] <- sum(diff(tmpx)*rollmean(tmpy,2))
    TILD[i] <- sum(diff(tmpx)*rollmean(tmplogy,2))
    
  }
  
  output <- cbind(PD,TILD)
  rownames(output) <- rownames(lineages)
  
  return(output)
  
}

TILD_full <- TILD(no_lineages_sim)


str(TILD_full)
TILD_full <- as.data.frame(TILD_full)

TILD_full$site <- rownames(TILD_full)

#check values, TILD should approximately equal pd
plot(TILD_full[,1],pd(WG_mat,phylo_tree_final)[,1])

TILD.df<- TILD_full%>% data.frame() %>% 
  mutate(site=rownames(TILD_full)) %>% 
  separate(site, c("x", "y"), sep = ":") %>% 
  mutate(across(c(x, y), parse_number)) %>% dplyr::select(-PD)

nrow(TILD.df)#1248

#-------------------------2.b. Measuring SR and PD-------------------------------------------
# converting to sparse matrix for faster calculation
WG_mat_sparse<- dense2sparse(WG_mat)
colnames(WG_mat_sparse)#470

#1 Calcualatig SR and PD values

SR_PD <- pd(WG_mat, phylo_tree_final, include.root=T) 

#-------------------------2.c. Measuring PE--------------------------------------------------
#Calculating PE; Refer to the main text for more info; We are using only those species 
#for which we have run the SDM for, 348 species. For PE Calculations we will be using this 
#subset of data.

biodiv.mat.clip

names(biodiv.mat.clip)#350

biodiv.mat.clip.348 <- biodiv.mat.clip
names(biodiv.mat.clip.348)

#Pruning  the tree
tree.348 <- phylo_tree_final
tree.348 <- prune.sample(biodiv.mat.clip.348[,-c(1:2)],tree.348)

#Checking if all the names match
setdiff(names(biodiv.mat.clip.348)[-c(1,2)],tree.348$tip.label )

nrow(biodiv.mat.clip.348) #1756
nrow(biodiv_final_all.subset) #1248

#Grid cells with less than 9 species have been removed so left_join this df with
#revised dataframe so that number of rows match

biodiv.mat.clip.348.subset <- left_join(biodiv_final_all.subset[,c(1,2)],
                                        biodiv.mat.clip.348)

WG_mat.348<- as.matrix(biodiv.mat.clip.348.subset[,-c(1,2)])
colnames(WG_mat.348)
rownames(WG_mat.348) <- paste0(biodiv.mat.clip.348.subset$x,":",biodiv.mat.clip.348.subset$y)


#This PE is not scaled, not scaling it so that it is easy to compare with other measures in terms of magnitude
#for scaled PE can use canaper.

PE_348_sp  <- phylo_endemism(dense2sparse(WG_mat.348), tree.348,weighted = TRUE)
PE_348_sp <- data.frame(PE_348_sp)

names(PE_348_sp)[1] <- "PE"


#Merging all the richness indices into a dataframe

WG_all_indices <-  data.frame(cbind(TILD.df,SR_PD,PE_348_sp))
str(WG_all_indices)
summary(WG_all_indices)
nrow(WG_all_indices)#1248
head(WG_all_indices)

WG_all_indices <- WG_all_indices %>% relocate(x,y,TILD,PD,PE,SR)



#-------------------------2.d.Plotting the evolutionary indices-------------------------------
plot_SR_WG <- WG_all_indices %>% ggplot (aes(x=x , y=y,fill=SR )) + 
  geom_raster()+
  geom_polygon(data=wg, aes(x=long, y=lat, group=group),size=0.7, color="black", fill=NA) +
  theme_bw(base_size = 15)+ coord_fixed()+ylab("Latitude")+xlab("Longitude")+
  scale_fill_paletteer_c("grDevices::YlOrRd", direction = -1,
                         breaks=c(min(WG_all_indices$SR),max(WG_all_indices$SR)))+theme_custom

plot_PD_WG<-    WG_all_indices %>%  ggplot (aes(x=x , y=y,fill=PD )) +
  geom_raster()+ 
  geom_polygon(data=wg, aes(x=long, y=lat, group=group), size=0.7, color="black", fill=NA) +
  theme_bw(base_size = 15)+coord_fixed()+
  ylab("Latitude")+xlab("Longitude") +
  scale_fill_paletteer_c("grDevices::YlOrRd", direction = -1,
                         breaks=c(round(min(WG_all_indices$PD),2), max(WG_all_indices$PD)))+theme_custom

plot_PE_350sp_WG <-WG_all_indices %>% mutate(PE=round(PE,2)) %>% ggplot (aes(x=x , y=y,fill=PE )) +
  geom_raster()+
  geom_polygon(data=wg, aes(x=long, y=lat, group=group), size=0.7, color="black", fill=NA) +
  theme_bw(base_size = 15)+
  coord_fixed()+ylab("Latitude")+xlab("Longitude")+
  scale_fill_paletteer_c("grDevices::YlOrRd", direction = -1,
                         breaks=c( round(min(WG_all_indices$PE),2)+0.01, #adding because rounding takes the min value which is not in the dataset
                                   round(max(WG_all_indices$PE),2)))+theme_custom

plot_TILD_WG <- WG_all_indices%>%  ggplot (aes(x=x , y=y,fill=TILD )) +
  geom_raster()+ 
  geom_polygon(data=wg, aes(x=long, y=lat, group=group), size=0.7, color="black", fill=NA) +
  theme_bw(base_size = 15)+
  coord_fixed()+ylab("Latitude")+xlab("Longitude")+
  scale_fill_paletteer_c("grDevices::YlOrRd", direction = -1,
                         breaks=c(round(min(WG_all_indices$TILD),2)+0.01,
                                  round(max(WG_all_indices$TILD),2)))+theme_custom

plot_SR_WG+plot_PD_WG+plot_PE_350sp_WG+plot_TILD_WG

#-------------------------2.e. Diversity measure wrt Biog zone--------------------------------------------

WG_all_indices
nrow(WG_all_indices)

WG_all_indices$biog <- "CWG"

WG_all_indices <- WG_all_indices %>% mutate(biog = ifelse(y <= 11,
                                                          yes =  "SWG",
                                                          no= ifelse(y >= 15.8, "NWG","CWG")))

biog_summ <- WG_all_indices %>% group_by(biog) %>% 
            summarise(mean(SR), mean(TILD), mean(PE), mean(PD),
            median(SR), median(TILD), median(PE), median(PD)) %>% as.data.frame()



#SWG vs NWG
biog_summ$`mean(SR)`[3]/biog_summ$`mean(SR)`[2] # 5.5
biog_summ$`mean(PD)`[3]/biog_summ$`mean(PD)`[2] #3.4
biog_summ$`mean(TILD)`[3]/biog_summ$`mean(TILD)`[2] #1.4
biog_summ$`mean(PE)`[3]/biog_summ$`mean(PE)`[2] # 6.3

#SWG vs CWG
biog_summ$`mean(SR)`[3]/biog_summ$`mean(SR)`[1] #1.5
biog_summ$`mean(PD)`[3]/biog_summ$`mean(PD)`[1] #1.4
biog_summ$`mean(TILD)`[3]/biog_summ$`mean(TILD)`[1] #1.1
biog_summ$`mean(PE)`[3]/biog_summ$`mean(PE)`[1] # 1.7

#CWG vs NWG
biog_summ$`mean(SR)`[1]/biog_summ$`mean(SR)`[2] #3.7
biog_summ$`mean(PD)`[1]/biog_summ$`mean(PD)`[2] #2.4
biog_summ$`mean(TILD)`[1]/biog_summ$`mean(TILD)`[2] #1.26
biog_summ$`mean(PE)`[1]/biog_summ$`mean(PE)`[2] # 3.6


WG_all_indices %>% mutate( biog=factor(biog , levels=c("SWG", "CWG", "NWG"))) %>% ggplot(aes(y=TILD, x=biog))+geom_boxplot()+theme_bw()+xlab("")+
  WG_all_indices %>% ggplot(aes(y=PD, x=biog))+geom_boxplot()+theme_bw()+xlab("")+
  WG_all_indices %>% ggplot(aes(y=PE, x=biog))+geom_boxplot()+theme_bw()+xlab("")+
  WG_all_indices %>% ggplot(aes(y=SR, x=biog))+geom_boxplot()+theme_bw()+xlab("")



#-------------------------2.f. Correlations between diversity indices----------------------

# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)}


pairs(WG_all_indices[,c(-1,-2,-7)], lower.panel = panel.cor)#removing x,y,biog

pairs.panels( WG_all_indices[,c(-1,-2,-7)] ,
              method = "pearson", # correlation method
              hist.col = "#00AFBB",
              density = TRUE,  # show density plots
              ellipses = F # show correlation ellipses
)

dev.off()

cor.test(WG_all_indices$SR,WG_all_indices$PD) 
cor.test(WG_all_indices$SR,WG_all_indices$PE)
cor.test(WG_all_indices$SR,WG_all_indices$TILD)
cor.test(WG_all_indices$PD,WG_all_indices$PE)
cor.test(WG_all_indices$PD,WG_all_indices$TILD)
cor.test(WG_all_indices$TILD,WG_all_indices$PE)

############################################################################################

# -----------------------------------------------------------------------------------------#
# ---------------------3.PD at diff time period, LTT at each lat and lineage age trends----#
# -----------------------------------------------------------------------------------------#

#-------------------------3.a. PD at diff time period---------------------------------------
#-----------------------------3.a.1.Tree slices---------------------------------------------
# Labelling all the key orders in the tree-
order_names <-  matrix(nrow=35,ncol=2)
order_names[,1] <- All_sp_checklist.node %>% dplyr::select(order) %>%distinct() %>% unlist(use.names = F)

for(i in 1: 35)#modify this for family/genus level labelling
{
  order_names[i,2] <-ifelse(is.null(getMRCA(phylo_tree_final,All_sp_checklist.node %>%
                                              filter(order==order_names[i]& !is.na(node)) %>%
                                              dplyr::select(node) %>% unlist(use.names = F))), yes="NA",
                            no=getMRCA(phylo_tree_final,All_sp_checklist.node %>%
                                         filter(order==order_names[i]) %>%
                                         dplyr::select(node) %>% unlist(use.names = F)))
}
order_names
order_names <- as.data.frame(order_names)
names(order_names) <- c("order", "node")


order_names$node <-  as.numeric(order_names$node)

order_names.sub <- order_names %>% filter(node!=is.na(node))

revts(ggtree(phylo_tree_final) %<+%   order_names+ #%<+% node_label +geom_label(aes(label = revised_lab)) +
        ggtitle(paste0( "0 Mya"))+ theme(plot.title = element_text(hjust = 0.5))+
        geom_label(aes(label = order),colour="red")+
        theme_tree2()+geom_tiplab(size=1))



#finding the node labels of the deeper nodes of the tree
node_label <-phylo_tree_final%>% as_tibble() %>%
  filter(!label %in% names(biodiv_final_all.subset)[c(-1,-2)]) %>%
  filter(!grepl("mrcaott", label)) %>% data.frame() %>%
  distinct(label)
node_label

node_label <- data.frame(label=node_label$label[-5])# this is the oldest for the current phylogeny not labelling this
node_label <-phylo_tree_final%>% as_tibble() %>% filter(label %in% node_label$label)
node_label$revised_lab <-node_label$label

#Renaming for standardisation
node_label$revised_lab[c(2,4,5:10)] <- c( "Eudicotyledons", "Campanulids" ,
                                          "Lamiids", "Gentianales",
                                          "Santalales",  "Rosales",
                                          "Myrtales",   "Magnoliales")


node_label <- node_label %>% filter(label!="Mesangiospermae")

new.tree <- list()
plot_tree <- list()
tree <- phylo_tree_final


age_list <- c(seq(10,130, by=10),135)

# The code for slicing the tree is modified from Daru et al. 2018.
# Here the tree is collapsed from the tips to the root at different time periods.
# The terminals are sticky and as such each  branch length has an identity of all the
# descendants, the branch length of the terminal nodes is
# 0 at at time period older than their divergence times.
# This outputs a tree specified at each time period

for( i in age_list ){
  new.tree[[i]]<-tree

  #specify here threshold age for collapsing nodes in the tree
  node.age<-i

  print(node.age)

  goal.length<-max(cophenetic(tree))/2-node.age

  while (round(max(cophenetic(new.tree[[i]]))/2, 10)>round(goal.length,10)){
    k<-nodeHeights(new.tree[[i]])

    w<-which(k[,2]==max(k[,2]))[1]

    if(k[w,2]>goal.length){
      z<-k[w,2]-goal.length
      new.tree[[i]]$edge.length[w]<-new.tree[[i]]$edge.length[w]-z
      if(new.tree[[i]]$edge.length[w]<0){
        new.tree[[i]]$edge.length[w]<-0
      }#end if
    }#end if


  }#end while

  plot(new.tree[[i]], show.tip=F, main=paste0(i, " Mya"))
  plot_tree[[i+1]] <- revts(ggtree(new.tree[[i]]) %<+% order_names %<+%node_label +
                              geom_label(aes(label = order)) +geom_label(aes(label = revised_lab)) +
                              ggtitle(paste0(i, " Mya"))+ theme(plot.title = element_text(hjust = 0.5)))

}

# Inputting the tree at current time period
plot_tree[[1]] <- ggtree(phylo_tree_final) %<+% order_names %<+% node_label +#%<+% family_names2+
  geom_label(aes(label = order)) +geom_label(aes(label = revised_lab)) +
  ggtitle(paste0( "0 Mya"))+ theme(plot.title = element_text(hjust = 0.5))+
  theme_tree2()

plot_tree <- plot_tree%>% discard(is.null)


#
rm(phylo_tree)

#plotting the PD per  tree and saving it in a list
PD.depth <- list()
for(i in age_list) { PD.depth[[i]] <- PD(WG_mat_sparse, new.tree[[i]])}

PD.depth<- PD.depth %>% discard(is.null)


PD.depth<- do.call(cbind,PD.depth) %>% as.data.frame()
PD.depth$site <- rownames(PD.depth)

PD.depth <- PD.depth %>%
  separate(site, c("long", "lat"), sep = ":") %>%
  mutate(across(c(long, lat), parse_number))

names(PD.depth)[1:14] <- c(10,20,30, 40,50,60,70,80,90,100,110,120,130,135)

plot_PD.depth <- list()

for(i in 1:14)
{
  temp_data <- data.frame(cbind(PD.depth[,c("long", "lat")],age= PD.depth[,i])) %>% mutate(age=round(age,2))

  temp_plot  <- ggplot (data=temp_data,aes(x=long , y=lat,fill=age )) +
    geom_raster()+ ggtitle(paste0("PD ",names(PD.depth)[i], " Mya")) +
    geom_polygon(data=wg, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) +
    theme_bw(base_size = 15)+
    coord_fixed()+theme_void()+
    scale_fill_paletteer_c("grDevices::YlOrRd", direction = -1,
                           breaks=c(round(min(temp_data$age),2),max(temp_data$age)),"PD")


  plot_PD.depth[[i+1]] <- temp_plot
}

rm(temp_data)
plot_PD.depth[[1]]  <- WG_all_indices %>%
  ggplot (aes(x=x , y=y,fill=PD )) +
  geom_raster()+ ggtitle("0 Mya")+
  geom_polygon(data=wg, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) +
  theme_classic(base_size = 15) +scale_fill_gradientn(colours=c("#FFFFFF",  "#FEC44F", "#FE9929", "#EC7014", "#CC4C02","#993404", "#662506"),
                                                      breaks=c( round(min(WG_all_indices$PD),2) ,max(WG_all_indices$PD)))+
  coord_fixed()+theme_void()


plot_PD.depth<- plot_PD.depth %>% discard(is.null)
#this is imp to remove the null values

plot_PD.depth[[1]]+ plot_PD.depth[[9]]+ plot_PD.depth[[5]]+ plot_PD.depth[[8]]


#
#----------for select age 135 120  90 60 30 10---------------------#

wrap_plots(
  plot_PD.depth[[15]]+ylab("Latitude")+xlab("Longitude")+theme_custom+
    theme(legend.key.size = unit(0.3, "cm"),legend.position=c(0.3,0.25),
          axis.title = element_text( size=20,face="bold",   family="serif"),
          axis.text = element_text( size=20,colour = "black", family = "serif" ),
          legend.title = element_text(size=16,face = "bold", family = "serif"),
          legend.text = element_text( size=16, family = "serif")),

  plot_PD.depth[[13]]+theme(legend.key.size = unit(0.3, "cm"),legend.position=c(0.3,0.25),legend.title = element_text(size=16,face = "bold", family = "serif"),legend.text = element_text( size=16, family = "serif")),
  plot_PD.depth[[10]]+theme(legend.key.size = unit(0.3, "cm"),legend.position=c(0.3,0.25),legend.title = element_text(size=16,face = "bold", family = "serif"),legend.text = element_text( size=16, family = "serif")),
  plot_PD.depth[[7]]+theme(legend.key.size = unit(0.3, "cm"),legend.position=c(0.3,0.25),legend.title = element_text(size=16,face = "bold", family = "serif"),legend.text = element_text( size=16, family = "serif")),
  plot_PD.depth[[4]]+theme(legend.key.size = unit(0.3, "cm"),legend.position=c(0.3,0.25),legend.title = element_text(size=16,face = "bold", family = "serif"),legend.text = element_text( size=16, family = "serif")),
  plot_PD.depth[[2]]+theme(legend.key.size = unit(0.3, "cm"),legend.position=c(0.3,0.25),legend.title = element_text(size=16,face = "bold", family = "serif"),legend.text = element_text( size=16, family = "serif")),
  nrow=1)


#------------------------------------------------------------------------------#
#-----------------------plotting PD VS Latitude
#------------------------------------------------------------------------------#
# plot_PD.depth.lat <- list()
# 
# for(i in 1:14)
# {
#   temp_data <- data.frame(cbind(PD.depth[,c("long", "lat")],age= PD.depth[,i]))
#   
#   temp_plot  <- ggplot (data=temp_data,aes(y=lat , x=age )) + geom_point() +theme_bw() + 
#     theme_bw(base_size = 15)+xlab(paste0("PD at ",names(PD.depth)[i], " Mya"))+
#     ylab("")
#   
#   plot_PD.depth.lat[[i+1]] <- temp_plot
# }
# 
# rm(temp_data)
# rm(temp_plot)
# 
# plot_PD.depth.lat[[1]]  <- ggplot (data=WG_all_indices,aes(y=y , x=PD ))+ geom_point() +theme_bw() + 
#   theme_bw(base_size = 15)+xlab(paste0("PD at 0 Mya"))+ylab("Latitude")
# 


#-----------------------------3.b. LTT plots for each lat bin---------------------------------------
names(biodiv_final_all.subset)

#Converting wide to long and summing up occurrence wrt to each lat bin
LTT    <- biodiv_final_all.subset %>% 
          pivot_longer(cols = Acrocarpus_fraxinifolius:Agasthiyamalai_pauciflora,
          names_to="species", values_to="occ") %>% filter(occ!=0) %>% 
          mutate(lat_bin=round(y)) %>% dplyr::select(-x,-y) %>% 
          group_by(lat_bin,species) %>% summarise(occ=sum(occ))

min(LTT$occ)
LTT$occ[LTT$occ>1] <- 1   #Converting into presence absence maps      
max(LTT$occ)

#converting back to wide format
LTT.final <- LTT %>% pivot_wider(names_from = "species", values_from = "occ",values_fill=0)
rownames(LTT.final) <- LTT.final$lat_bin

#Similar function used for TILD, this counts the number of lineages in each bin for each
#time period

no_lineages.ltt <- function(commat,phy,bin_size) {
  
  tmp_seq <- rev(seq(0,max(branching.times(phy),by=bin_size)))
  matrix_out <- matrix(NA,nrow(commat),length(tmp_seq))
  rownames(matrix_out) <- rownames(commat)
  colnames(matrix_out) <- tmp_seq
  
  for (i in 1:nrow(commat)) {
    
    tmp_spp <- colnames(commat)[which(commat[i,]>0)]
    tmp_phy <- drop.tip(phy,which(!phy$tip.label%in%tmp_spp))
    tmp_btimes <- branching.times(tmp_phy)
    
    for (j in 1:ncol(matrix_out)) {
      
      matrix_out[i,j] <- length(which(tmp_btimes > tmp_seq[j]))
      
    }
    
  }
  
  matrix_out <- matrix_out + 1
  matrix_out[,1] <- 1
  
  matrix_out[,ncol(matrix_out)] <- apply(commat,1,sum)
  
  return(matrix_out)
  
}

LTT_lat <- as.data.frame(no_lineages.ltt(LTT.final,phylo_tree_final,1))

LTT_lat$lat.bin <- rownames(LTT_lat)

LTT_lat.long <- LTT_lat %>% pivot_longer(cols="0":"135",names_to = "age", values_to = "num")
str(LTT_lat.long)
LTT_lat.long$age<- as.character(LTT_lat.long$age)
LTT_lat.long$age<- factor(LTT_lat.long$age,levels=unique(LTT_lat.long$age))  

LTT_lat.long$lat.bin<- as.numeric(LTT_lat.long$lat.bin)

str(LTT_lat.long)

LTT_lat_plot <- LTT_lat.long %>% 
                ggplot(aes(x = reorder(age, desc(age)),y=log10(num), colour=lat.bin, group=lat.bin))+
                geom_line()+xlab("Age")+ylab(expression(paste(log[10],"(Num of Lineages)")))+
                scale_color_paletteer_c("grDevices::YlOrRd", direction = -1,name = "Latitudnal \nbins") +
                geom_point( data= LTT_lat.long %>% filter(age==0),aes(x=age, y=log10(num)))+
                geom_text_repel(data=LTT_lat.long %>% filter(age==0),aes(x=age, y=log10(num),label=lat.bin),direction = c("both", "y", "x"),size=8,max.overlaps = getOption("ggrepel.max.overlaps", default = 14))+
                scale_x_discrete(breaks=c(135,seq(130,10, -10),0))+  theme_classic(base_size = 30)+
                theme(legend.key.size = unit(1, "cm"),
                      axis.title = element_text(size=20,face="bold", family="serif"),
                      axis.text = element_text(size=20,colour = "black",family = "serif" ),
                      legend.title = element_text( size=24,face = "bold", family = "serif"),
                      legend.text = element_text(size=24,family = "serif"),
                      legend.position=c(0.9,0.3))
              
LTT_lat_plot



#-----------------------------3.c. lineage age trends-----------------------------------------
# This code is modified from Griffiths et al. 2021
# there are two functions, 1. age_lat_distn: to create the latitudinal distribution and 
# 2. age_profile to plot the peak richness


# function to create latitudinal distribution of lineages at different time slices
# assumes that tip labels in phylogeny are in same order as columns in commat
# metadata object must have column named "latitudinal" and column named "grid" 

age_lat_distn <- function(phylog,commat,metadata,age) {
  
  groups <- distconnected(as.dist(cophenetic(phylog)),toolong = age*2)
  
  unique_groups <- unique(groups)
  
  min_lat <- vector("numeric",length(unique_groups))
  max_lat <- vector("numeric",length(unique_groups))
  mean_lat <- vector("numeric",length(unique_groups))
  
  for (i in 1:length(unique_groups)) {
    
    tmp_commat <- as.matrix(commat[,which(groups==unique_groups[i])])
    grid <- rownames(tmp_commat)[which(apply(tmp_commat,1,sum)>0)]
    Lat <- metadata$lat[which(metadata$grid%in%grid)]
    Lat <- as.numeric( Lat)
    
    min_lat[i] <- min(Lat)
    max_lat[i] <- max(Lat)
    mean_lat[i] <- mean(Lat)
    
  }
  
  out_table <- as.data.frame(cbind(unique_groups,mean_lat,min_lat,max_lat))
  out_table <- out_table[order(out_table$mean_lat),]
  
  #Uncomment to see the plots at each time period
  plot(1:length(unique_groups),mean_lat,data=out_table,type="n",ylab="Latitude",xlab="")
  segments(1:length(unique_groups),out_table$min_lat,1:length(unique_groups),out_table$max_lat)
  
  return(out_table)
  
}

#function using output table from above age_lat_distn function to count number of 
#lineages at given latitudinal intervals for diff evolutionary age

age_profile <- function(in_mat,lat_interval=1) {
  
  lat_slices <- seq(min(in_mat$min_lat),max(in_mat$max_lat),lat_interval)
  
  number_of_lineages <- vector("numeric",length(lat_slices))
  
  for (i in 1:length(lat_slices)) {
    
    number_of_lineages[i] <- length(which(lat_slices[i] >= in_mat$min_lat & lat_slices[i] <= in_mat$max_lat))
    
  }
  
  out_table <- as.data.frame(cbind(lat_slices,number_of_lineages))
  
  return(out_table)
  
}

#Setting up the data
#assumes that tip labels in phylogeny are in same order as columns in commat
setdiff(colnames(WG_mat),phylo_tree_final$tip.label)

WG_mat <- as.matrix(WG_mat)
#Commat = WGmat.reorder has to be a matrix..

WG.mat.reorder <-WG_mat[,match(phylo_tree_final$tip.label,colnames(WG_mat))]
str(WG.mat.reorder)
#Checking
setdiff(colnames(WG.mat.reorder),phylo_tree_final$tip.label)

max(branching.times(phylo_tree_final)) #135.91

#time slices
age_list <- c(seq(10,130, by=10),135)
time_slices <- c(0,age_list)
time_slices <- replace(time_slices, time_slices==0, 1)

#empty vectors to fill
peak_mean <- vector("numeric",length(time_slices))
peak_max <- vector("numeric",length(time_slices))
peak_min <- vector("numeric",length(time_slices))

# this loop plots for latitudinal range of lineages at diff time slices with a plot of richness trend

for (i in 1:length(time_slices)) {
  
  ages_tmp <- age_lat_distn(phylo_tree_final,WG.mat.reorder,
                            as.data.frame(cbind(grid=paste0(WG_pred$x,":",WG_pred$y),
                                                lat=as.numeric(WG_pred$y))),
                            time_slices[i])
  
  ages_tmp_out <- age_profile(ages_tmp,1)
  #plot(ages_tmp_out, type="l")
  peak_ages_tmp <- (ages_tmp_out$lat_slices[which(ages_tmp_out$number_of_lineages==
                                                    max(ages_tmp_out$number_of_lineages))])
  peak_ages_tmp
  peak_mean[i] <- mean(peak_ages_tmp)
  peak_max[i] <- max(peak_ages_tmp)
  peak_min[i] <- min(peak_ages_tmp)
  
}

#saving the output as df
ages <- as.data.frame(cbind(time_slices,peak_mean,peak_max,peak_min))


# This function uses the above two function, age_lat_distn and age_profile and formats it for
# the table

Age_lat_table<- function(Age){
  age=Age
  
  temp_age <- age_lat_distn(phylo_tree_final,WG.mat.reorder,
                            as.data.frame(cbind(grid=paste0(WG_pred$x,":",WG_pred$y),
                                                lat=as.numeric(WG_pred$y))),  age)
  temp_age_out <- age_profile(temp_age,1)
  
  
  temp <- as.data.frame(cbind(Age=age,
                              temp_age %>% filter(unique_groups==max(unique_groups)) %>% dplyr::select(unique_groups),
                              temp_age %>% filter(max_lat <= 13) %>% summarise(lt_13=n()),
                              temp_age %>% filter(min_lat > 13) %>% summarise(gt_13=n()),
                              #ages %>% filter(time_slices==135) %>% dplyr::select(peak_max)
                              temp_age_out %>% filter(lat_slices %in% (ages %>% filter(time_slices==age) %>% 
                                                                         dplyr::select(peak_max))) %>% 
                                rename(peak_lat=lat_slices, lineage_peak_lat=number_of_lineages)))
  return(temp)
  
}

Age_lineage_table_full <- lapply(rev(age_list),Age_lat_table)
Age_lineage_table_full <- as.data.frame(Reduce(rbind,Age_lineage_table_full))


#Outputting the table at all time periods
Age_lineage_table_full %>%   gt()%>% 
  fmt_number(
    columns = peak_lat,decimal=1,
    suffixing = TRUE) %>% cols_label(unique_groups=md("Total lineages"),
                                     Age=md("Age (Mya)"), lt_13=md("No. lineages only < 13 deg"),
                                     gt_13=md("No. lineages only > 13 deg"),peak_lat=md("Latitude of peak richness"),
                                     lineage_peak_lat=md("Number of lineages at peak richness"))


#Outputting the table at select time periods
temp_age_list <- c(10,30,60,90,120,135)

Age_lineage_table_select<- lapply(rev(temp_age_list),Age_lat_table)
Age_lineage_table_select <- as.data.frame(Reduce(rbind,Age_lineage_table_select))


Age_lineage_table_select %>%   gt()%>% 
  fmt_number(
    columns = peak_lat,decimal=1,
    suffixing = TRUE) %>% cols_label(unique_groups=md("Total lineages"),
                                     Age=md("Age (Mya)"), lt_13=md("No. lineages only < 13 deg"),
                                     gt_13=md("No. lineages only > 13 deg"),peak_lat=md("Latitude of peak richness"),
                                     lineage_peak_lat=md("Number of lineages at peak richness"))

#Similar to Age_lat_table, this function takes the two functions, age_lat_distn and age_profile, and plots the output 
# plots in ggplot in the desired sequence

plot_age_lat <- function(age,label)
{  
  temp_age <- age_lat_distn(phylo_tree_final,WG.mat.reorder,
                            as.data.frame(cbind(grid=paste0(WG_pred$x,":",WG_pred$y),
                                                lat=as.numeric(WG_pred$y))),  age)
  temp_age_out <- age_profile(temp_age,1)
  
  
  
  p <- wrap_plots(temp_age %>% ggplot(aes(x=1:length(unique_groups),y=mean_lat))+ylab("")+xlab("")+theme_bw()+
                    geom_segment(aes(x=1:length(unique_groups),y=min_lat,xend=1:length(unique_groups),yend=max_lat))+theme_custom2+
                    ggtitle(paste0(label,"     ",age, "Mya"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()),
                  temp_age_out %>% ggplot(aes(number_of_lineages,lat_slices))+geom_path()+xlab("")+ylab("")+theme_bw(),ncol=1)+
    theme_custom2+xlim(0,max(temp_age_out$number_of_lineages))
  return(p)
  
}


label <- c(letters[seq( from = 1, to = 6 )])

wrap_plots( mapply(plot_age_lat,rev(temp_age_list),label=label, SIMPLIFY = F),nrow = 1)

#------------------------------------------------------------------------------------------#

#-----------------------------3.d. Turnover-----------------------------------------
#converting the matrix into lat bins and reordering them

WG_mat_lat_bin <- WG_mat %>% as.data.frame() %>% 
                  mutate(grid=rownames(WG_mat)) %>% 
                  separate(grid, c("long", "lat"), sep = ":") %>% 
                  mutate(across(c(long, lat), parse_number)) %>% 
                  dplyr::select(-long) %>% 
                  pivot_longer(cols = colnames(WG_mat)[1]:
                  colnames(WG_mat)[length(colnames(WG_mat))],
                  names_to="species", values_to="occ") %>% filter(occ!=0) %>% 
                  mutate(lat_bin=round(lat)) %>% dplyr::select(-lat) %>% 
                  group_by(lat_bin,species) %>% summarise(occ=sum(occ))

min(WG_mat_lat_bin$occ)
WG_mat_lat_bin$occ[WG_mat_lat_bin$occ>1] <- 1   #Converting into presence absence maps      
max(WG_mat_lat_bin$occ)

WG_mat_lat_bin.final <- WG_mat_lat_bin %>% pivot_wider(names_from = "species", values_from = "occ",values_fill=0) %>% as.data.frame()
rownames(WG_mat_lat_bin.final) <- WG_mat_lat_bin.final$lat_bin 


#Commat = WGmat.reorder has to be a matrix..

WG_mat_lat_bin.final.reordered <- WG_mat_lat_bin.final[,match(phylo_tree_final$tip.label,
                                                              colnames(WG_mat_lat_bin.final))]




age_lat_distn_turnover <- function(phylog,commat,metadata,age) {
  
  groups <- distconnected(as.dist(cophenetic(phylo_tree_final)),toolong = age*2)
  
  unique_groups <- unique(groups)
  
  temp_list <- list()
  
  for (i in 1:length(unique_groups)) {
    
    tmp_commat  <- as.matrix(commat[which(groups==unique_groups[i])],rownames = TRUE)
    lat <- rownames(tmp_commat)[which(apply(tmp_commat,1,sum)>0)]
    lat <- metadata$lat_bin[which(metadata$lat_bin%in%lat)]
    
    tmp_commat.tmp <- left_join(metadata %>% dplyr::select(lat_bin) %>% rename(lat=lat_bin),
                                as.data.frame(cbind(lat,
                                                    unique_groups[i])))
    
    names(tmp_commat.tmp)[2]   <-    unique_groups[i]
    tmp_commat.tmp[is.na(tmp_commat.tmp)] <- 0
    temp_list[[i]] <- tmp_commat.tmp[,2]
    
  }
  
  tmp_commat_final <-  cbind(metadata$lat_bin,
                             do.call(cbind,temp_list)) %>% as.data.frame()
  names(tmp_commat_final) <- c("Lat",unique_groups)
  
  
  return(tmp_commat_final)
  
}



age_list_select <- c(135,120,90,60,30,10)
age_lat_distn_turnover_all <- list()
age_lat_distn_turnover_all_pres_abs <- list()

for(i in 1: length(age_list_select))
{
  
  age_lat_distn_turnover_all[[i]] <- age_lat_distn_turnover(phylo_tree_final,
                                                            WG_mat_lat_bin.final.reordered,
                                                            WG_mat_lat_bin %>% dplyr::select(lat_bin) %>% distinct(lat_bin),
                                                            age_list_select[i])
  #converting to pres abs
  age_lat_distn_turnover_all_pres_abs[[i]] <- age_lat_distn_turnover_all[[i]]
  
  age_lat_distn_turnover_all_pres_abs[[i]] [age_lat_distn_turnover_all_pres_abs[[i]]>1] <- 1
  
}


lineage_turnover_lat_order <- list()

set.seed(10)

#lineage_nestedess

for (i in 1:length(age_list_select))
{
  tmp_turnover  <- Turnover(age_lat_distn_turnover_all_pres_abs[[i]][,-1], 
                            method = "r1", sims = 1000, scores = 1,
                            order = T, orderNulls = FALSE, allowEmpty = FALSE,
                            binary = TRUE, verbose = T, seed = 1, 
                            fill = TRUE)
  lineage_turnover_lat_order[[i]]   <- cbind(age=age_list_select[i],  num_lineages=length(colnames(age_lat_distn_turnover_all_pres_abs[[i]])[-1]),
                                             tmp_turnover)
  
  
  print(c(i, age_list_select[i]))
  
}

lineage_turnover_tmp <- do.call(rbind,lineage_turnover_lat_order)


lineage_turnover_final <- lineage_turnover_tmp %>% filter(name == "turnover"| name== "simMean"| name=="p") %>%
  pivot_wider(names_from = "name", values_from = "stat") %>% 
  rename(obs_turnover=turnover, mean_turnover=simMean) %>% 
  relocate(age, num_lineages, obs_turnover, mean_turnover, p) %>% 
  mutate(mean_turnover=round(mean_turnover,2), p=round(p,3))

lineage_turnover_final %>%   gt()%>% 
  cols_label(age=md("Age (Mya)"),
             num_lineages=md("Num. of lineages"),
             obs_turnover=md("Observed turnover"),
             mean_turnover=md("Expected Turnover")) %>% 
  tab_options(table.font.names = "Times New Roman")


#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
############################################################################################
# -----------------------------------------------------------------------------------------#
#----------------------4.Order family genus contribution to PD-----------------------------#
# -----------------------------------------------------------------------------------------#
#-------------------------4.a. Contribution of Superorder to PD at each time period----------------------
sort(unique(All_sp_checklist.node$order)) #35 orders

#Creating the superorder list 
Superasteridae <-    c("Apiales", "Aquifoliales","Asterales","Lamiales","Dipsacales","Ericales",
                       "Cornales","Santalales",
                       "Dilleniales","Gentianales","Icacinales","Metteniusales")

Monocotyledoneae <-  c("Arecales","Asparagales","Pandanales")
Superrosidae <-      c("Brassicales","Malvales","Myrtales","Cucurbitales","Celastrales",
                       "Malpighiales","Oxalidales",
                       "Fabales","Rosales","Sapindales","Crossosomatales")
Basal_Eudicots <-    c("Vitales","Saxifragales","Proteales","Buxales","Sabiales","Ranunculales")
Magnoliidae <-       c("Laurales","Magnoliales")

#monocot distribution none found above 15 in this dataset
LTT.final %>% dplyr::select(c(Caryota_urens,Pinanga_dicksonii,Corypha_umbraculifera,Arenga_wightii,Pandanus_furcatus,Dracaena_elliptica))

#Creating the dataframe
superorder <- All_sp_checklist.node %>% dplyr::select(label, family,order)

superorder$super_order <- ""  
superorder$super_order[superorder$order %in% Superasteridae ] <- "Superasteridae"
superorder$super_order[superorder$order %in% Superrosidae ] <- "Superrosidae"
superorder$super_order[superorder$order %in% Monocotyledoneae ] <- "Monocotyledoneae"
superorder$super_order[superorder$order %in% Magnoliidae ] <- "Magnoliidae"
superorder$super_order[superorder$order %in% Basal_Eudicots ] <- "Basal Eudicots"
superorder$super_order[superorder$super_order==""] <- "Chloranthales"


#Excluding Cholranthales as it has only one entry of Ehretia
# label       family       order   super_order
# 1 Ehretia_canarensis Boraginaceae Boraginales Chloranthales

key.super_orders <- c("Superasteridae",     "Superrosidae",   
                      "Monocotyledoneae",  "Basal Eudicots", "Magnoliidae")

#using Ltt.final as the data is summarised to lat. bins
LTT.final

#New tree1 adding the current unsliced tree
new.tree[[1]] <- phylo_tree_final
Lat_comm <- LTT.final  

# Adding the biog zones
Lat_comm$biog <- "CWG"
Lat_comm <- Lat_comm %>% mutate(biog=ifelse(lat_bin <=11, "SWG",
                                            no=ifelse(lat_bin>=15.8,"NWG","CWG")))
Lat_comm %>% dplyr::select(lat_bin,biog)
names(Lat_comm)

Lat_comm <- Lat_comm[,c(1,472,2:471)]

#creating a new dataframe by pooling lats at biog zones
Lat_biog_long <- Lat_comm %>% pivot_longer(cols = names(Lat_comm)[3]:names(Lat_comm)[472],
                                           names_to="species", values_to="occ") %>%
                filter(occ!=0) %>% 
                group_by(biog,species) %>% summarise(occ=sum(occ))

Lat_biog_long$occ[Lat_biog_long$occ>1] <- 1 

Lat_biog_wide <-Lat_biog_long %>% pivot_wider(names_from = "species", values_from = "occ",values_fill=0)


#Function which inputs biog zone and age, and calculates the PD for each superorder at each age
# can be modified to examine lat level trends

lat_age_branch_len<- function(biog_zone,age){
  
  
  #Lat_comm <- Lat_comm %>% filter(lat_bin==lat)# gives the comm for the lat bin
  #Use the above if you need for lat trends
  #Lat_comm_temp <- Lat_biog_wide
  
  biog <- biog_zone
  Lat_comm_final <- Lat_biog_wide %>% filter(biog==biog_zone)
  Lat_comm_final <- Lat_comm_final[,-c(1,2)]
  Lat_comm_final <- Lat_comm_final[, colSums(Lat_comm_final) > 0]# removing species with no occ
  
  Lat_age_tree <-  prune.sample(Lat_comm_final,new.tree[[age]])# new tree age is sliced tree at respective age
  
  #temp_tree_BL is the Tree for the given age and latitude
  
  #using the same logic as the order/family pd contribution 
  
  superorderorder.PD.contribution <- list()
  
  for( i in 1:5)#length(unique(superorder$super_order))) #
  {
    #labels which corresond to the select orders
    tmp.label <- superorder %>% filter(super_order==key.super_orders[i]) %>% distinct(label) 
    
    #selecting the label at the prt lat wrt the order 
    tmp_super.order.lat.age <- tmp.label$label[tmp.label$label %in% names(Lat_comm_final)]
    
    if(length(tmp_super.order.lat.age)<1) next # 
    select_comm <- Lat_comm_final %>%   dplyr::select(all_of(tmp_super.order.lat.age))
    
    #tmp.Ltt.sparse <- dense2sparse( select_comm)
    
    superorderorder.PD.contribution[[i]]  <- as.data.frame(cbind(biog=biog,#lat_bin=lat, 
                                                                 order=key.super_orders[i],
                                                                 sum_bl= prune.sample(select_comm,
                                                                                      Lat_age_tree) %>% 
                                                                   as_tibble() %>% 
                                                                   summarise(sum_bl=sum(branch.length,na.rm=T))))
  }
  
  superorderorder.PD.contribution <- superorderorder.PD.contribution %>% discard(is.null)
  superorderorder.PD.contribution <- do.call(rbind,superorderorder.PD.contribution)
  
  Lat_age_PD <-    cbind( Age=age,
                          superorderorder.PD.contribution,
                          pd(Lat_comm_final,Lat_age_tree)[1])  # Lat_age_tree %>% as_tibble() %>%  summarise(Tot_branch_length=sum(branch.length,na.rm=T))
  
  names(Lat_age_PD)[5] <- "Tot_pd"
  Lat_age_PD$prop_pd <- Lat_age_PD$sum_bl/Lat_age_PD$Tot_pd
  
  
  
  return(Lat_age_PD)
  
}

Order_lat_age_PD_final <- 
  rbind(do.call(rbind,mapply(lat_age_branch_len,c("NWG"),#lat=8:19,
                             age=c(1,age_list),SIMPLIFY = F)),
        do.call(rbind,mapply(lat_age_branch_len,c("CWG"),#lat=8:19,
                             age=c(1,age_list),SIMPLIFY = F)),
        do.call(rbind,mapply(lat_age_branch_len,c("SWG"),#lat=8:19,
                             age=c(1,age_list),SIMPLIFY = F)))

#adding zero for Monocots in NWG
Order_lat_age_PD_final <- rbind(Order_lat_age_PD_final,cbind(Age=as.numeric(unique(Order_lat_age_PD_final$Age)),biog="NWG",
                                                             order="Monocotyledoneae",sum_bl=0,Tot_pd=0,prop_pd=0))

Order_lat_age_PD_final$Age <- as.numeric(Order_lat_age_PD_final$Age)
Order_lat_age_PD_final$Age[Order_lat_age_PD_final$Age==1] <- 0

wrap_plots(
  Order_lat_age_PD_final %>% filter(biog=="SWG") %>% group_by(Age,order) %>% 
    mutate(Age=as.factor(Age), prop_pd=as.numeric(prop_pd))  %>% 
    ggplot(aes( Age,order, fill= prop_pd)) + ggtitle("SWG 8-11 deg")+scale_x_discrete(limits = rev)+
    geom_tile()+  scale_fill_gradient(low="#FFFFE5", high="#662506", name= "Prop. PD") +theme_bw(base_size = 15)+ylab("")+
    theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 15), 
          axis.text = element_text(colour = "black")),
  
  Order_lat_age_PD_final %>% filter(biog=="CWG") %>% group_by(Age,order) %>% 
    mutate(Age=as.factor(Age), prop_pd=as.numeric(prop_pd))  %>% 
    ggplot(aes( Age,order, fill= prop_pd)) + ggtitle("CWG 12-15 deg")+scale_x_discrete(limits = rev)+
    geom_tile()+  scale_fill_gradient(low="#FFFFE5", high="#662506", name= "Prop. PD") +theme_bw(15)+ylab("")+
    theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(),axis.text = element_text(colour = "black")),
  
  Order_lat_age_PD_final %>% filter(biog=="NWG") %>% group_by(Age,order) %>% 
    mutate(Age=factor(Age), prop_pd=as.numeric(prop_pd)) %>% 
    ggplot(aes( Age,order, fill= prop_pd)) + ggtitle("NWG 16-19 deg")+scale_x_discrete(limits = rev)+
    geom_tile()+  scale_fill_gradient(low="#FFFFE5", high="#662506", name= "Prop. PD") +theme_bw(base_size = 15)+ylab("")+
    theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
          axis.text = element_text(colour = "black")))


#------------------------------------------------------------------------------#
#-------------------------4.b. Key family  contributing to PD per biog-----------------------------
#function to calculate PD of order/family/genus per biog zone
# selecting the key families 

sing_sp_fam <-  All_sp_checklist.node %>% group_by(family,label) %>% summarise(num=n()) %>% group_by(family) %>%
  summarise(num=sum(num)) %>%  filter(num==1)

#Filtering  families, genera which are not single

key.family <-  All_sp_checklist.node %>% distinct(family) %>% filter(!family %in% sing_sp_fam$family)

tax_group_PD_lat <- function(taxon_group,key.group)  {
  
  tmp.label <- All_sp_checklist.node %>% filter(All_sp_checklist.node[[taxon_group]]==key.group) %>% distinct(label) 
  
  PD_Lat_biog <- biodiv_final_all.subset %>% pivot_longer(cols = Acrocarpus_fraxinifolius:Agasthiyamalai_pauciflora,
                                                          names_to="species", values_to="occ") %>% filter(occ!=0) %>% 
    mutate(lat_bin=round(y)) %>% dplyr::select(-x,-y) %>% 
    group_by(lat_bin,species) %>% 
    mutate(biog=ifelse(lat_bin <=11, "SWG",#comment this if you need all lat trends
                       no=ifelse(lat_bin>=15.8,"NWG","CWG"))) %>% 
    group_by(biog,species) %>% 
    summarise(occ=sum(occ))
  
  PD_Lat_biog$occ[PD_Lat_biog$occ>1] <- 1  
  
  PD_Lat_biog.final <- PD_Lat_biog  %>% pivot_wider(names_from = "species", values_from = "occ",values_fill=0)
  
  PD_Lat_biog.final.contribution  <- as.data.frame(cbind(lat_bin=PD_Lat_biog.final[,1],
                                                         taxonomic_group=key.group,
                                                         #Calculating PD only for the select labels of a family/order
                                                         pd=PD(dense2sparse(PD_Lat_biog.final[,-1]  %>% 
                                                                              dplyr::select(tmp.label$label)),
                                                               prune.sample(PD_Lat_biog.final[,-1]  %>% 
                                                                              dplyr::select(tmp.label$label),
                                                                            phylo_tree_final))))
  return(PD_Lat_biog.final.contribution)
}

family_PD_biog_lat <- mapply(tax_group_PD_lat,"family",key.family$family, SIMPLIFY = F)
family_PD_biog_lat <- (do.call(rbind,family_PD_biog_lat))




#-------------------------4.c. Key families contributing to PD per lat bins-----------------------------
family.PD.lat.df <- biodiv_final_all.subset %>% pivot_longer(cols = Acrocarpus_fraxinifolius:Agasthiyamalai_pauciflora,
                                                             names_to="species", values_to="occ") %>% filter(occ!=0) %>% 
                     mutate(lat_bin=round(y)) %>% dplyr::select(-x,-y) %>% group_by(lat_bin,species) %>% 
                     summarise(occ=sum(occ))

min(family.PD.lat.df$occ)
family.PD.lat.df$occ[family.PD.lat.df$occ>1] <- 1         
max(family.PD.lat.df$occ)

family.PD.lat.final <-family.PD.lat.df %>% pivot_wider(names_from = "species", values_from = "occ",values_fill=0)
rownames(family.PD.lat.final) <- family.PD.lat.final$lat_bin



family.PD <- list()

for( i in 1:nrow(key.family))
{
  tmp.label <- All_sp_checklist.node %>% filter(family==key.family[i,]) %>% distinct(label) 
  
  family.tree  <- prune.sample(family.PD.lat.final[,-1]  %>% dplyr::select(tmp.label$label),
                               phylo_tree_final)
  family.PD[[i]] <-   cbind(family=key.family[i,],PD=sum(family.tree$edge.length)) %>% 
                      as.data.frame() %>% mutate(PD=as.numeric(PD))
  
}

family.PD.final <- (do.call(rbind,family.PD))

#Total conntribution of PD by all families
family.PD.final %>% summarise(sum(PD))  #9663.718
#Total PD of the entire tree
sum(phylo_tree_final$edge.length) #15158.58

family_PD_biog_lat%>% arrange(desc(pd)) %>% distinct(taxonomic_group,.keep_all = T)
family.PD.final %>% arrange(desc(PD)) %>% filter(PD >350)
#        family       PD
# 1      Annonaceae 800.9482
# 2       Myrtaceae 766.5849
# 3   Euphorbiaceae 644.5336
# 4       Lauraceae 604.9198
# 5       Rubiaceae 482.4925
# 6       Meliaceae 462.0323
# 7  Phyllanthaceae 452.2757
# 8        Fabaceae 450.3357
# 9   Anacardiaceae 400.6749
# 10       Rutaceae 397.3416

#10 families which have greater than 350 pd 

family.PD.final %>% arrange(desc(PD)) %>% filter(PD >350) %>% summarise(select_pd=sum(PD))/
                    family.PD.final %>% summarise(tot=sum(PD))
#0.56 total PD of *familes

dom_family <- family.PD.final %>% arrange(desc(PD)) %>% filter(PD >300) %>% distinct(family)

Key_family_biog_PD_plot <- family_PD_biog_lat %>% group_by(biog) %>% mutate(tot_pd=sum(pd),) %>% 
  rename("family"=taxonomic_group) %>% 
  mutate(prop_pd=pd/tot_pd) %>% filter(family %in% dom_family$family) %>% 
  filter(pd!=0) %>% #Fabaceae is removed at NWG
  mutate(biog = factor(biog, levels=c("SWG","CWG","NWG"))) %>% 
  ggplot(aes(fill=family, y=prop_pd, x=biog)) +
  geom_bar(position="stack", stat="identity", colour="black")+theme_bw(base_size = 15)+coord_flip()+xlab("")+
  ylab("Prop. of PD ")+geom_text(aes(y = prop_pd, label = family),
                                 colour = "black",position=position_stack(vjust=0.5), stat="identity", size=5,
                                 fontface="bold")+theme(legend.position = "none",#"top",legend.direction = "horizontal",
                                                        axis.text = element_text(colour = "black", size=25),
                                                        axis.title.x  = element_text(colour = "black", size=25))+ 
                                 scale_fill_manual(values=c(rep("#FFFFE5",10)))+#scale_fill_manual(values=c(colorRampPalette(c("#662506","#FFFFE5"))(10)))+
                                  guides(fill = guide_legend(nrow = 2, title=NULL)) 

#Adding number of families
family_PD_biog_lat_plot <- family_PD_biog_lat %>% group_by(biog) %>% mutate(tot_pd=sum(pd),) %>% 
                          rename("family"=taxonomic_group) %>% group_by(biog) %>% filter(pd!=0) %>%
                          mutate(Num_fam=n_distinct(family))

family_PD_biog_lat_plot<- left_join(family_PD_biog_lat %>% group_by(biog) %>% mutate(tot_pd=sum(pd),) %>% 
                                      rename("family"=taxonomic_group) %>% 
                                      mutate(prop_pd=pd/tot_pd) %>% filter(family %in% dom_family$family) %>% 
                                      filter(pd!=0) %>% #Fabaceae is removed at NWG
                                      mutate(biog = factor(biog, levels=c("SWG","CWG","NWG"))),
                                    family_PD_biog_lat %>% group_by(biog) %>%   rename("family"=taxonomic_group) %>% 
                                      filter(pd!=0) %>% summarise(num=n_distinct(family))) %>% as.data.frame()

Key_family_biog_PD_plot <- Key_family_biog_PD_plot+ geom_text(data=family_PD_biog_lat_plot,
                                                              aes(label= paste("Num = ",num)),position="identity", 
                                                              stat="identity", size=5,y=0.584,fontface="bold",family="Times")


rm(family_PD_biog_lat_plot)# used only for num familes

#-------------------------4.d. Summarizing Family level and species information----


All_sp_checklist.node %>% group_by(family) %>% mutate(Endemic=as.numeric(Endemic))%>% summarise(num_sp=n_distinct(label),
                                                                                                endem=sum(Endemic,na.rm=T)) %>% arrange(desc(num_sp))


All_sp_checklist.node %>% group_by(order,family) %>% filter(family %in% dom_family$family) %>%
  mutate(Endemic=as.numeric(Endemic))%>% summarise(num_sp=n_distinct(label),
                                                   endem=sum(Endemic,na.rm=T)) %>% arrange(desc(num_sp)) %>%
  mutate(prop_end=round(endem/num_sp,2)) %>% arrange(desc(prop_end)) %>%
  as.data.frame()
#            order         family num_sp endem prop_end
# 1   Magnoliales     Annonaceae     26    23     0.88
# 2      Myrtales      Myrtaceae     29    25     0.86
# 3      Laurales      Lauraceae     43    35     0.81
# 4    Sapindales  Anacardiaceae     16    13     0.81
# 5   Gentianales      Rubiaceae     43    34     0.79
# 6       Fabales       Fabaceae     10     7     0.70
# 7  Malpighiales  Euphorbiaceae     27    15     0.56
# 8  Malpighiales Phyllanthaceae     16     8     0.50
# 9    Sapindales       Rutaceae     13     5     0.38
# 10   Sapindales      Meliaceae     19     6     0.32

fam_endem_list <- All_sp_checklist.node %>% group_by(order,family) %>% 
                  mutate(Endemic=as.numeric(Endemic))%>% 
                  summarise(num_sp=n_distinct(label),
                            endem=sum(Endemic,na.rm=T)) %>% arrange(desc(num_sp)) %>%
                  mutate(prop_end=round(endem/num_sp,2)) %>% 
                  arrange(-num_sp,-prop_end) %>% as.data.frame()


sp_endem_list <- All_sp_checklist.node %>% group_by(order,family, label) %>% arrange(order) %>% 
                 dplyr::select(order,family, genus, label, Endemic)


############################################################################################
############################################################################################
#------------------------------------------------------------------------------------------#
#----------------------5.Examining the correlates of evolutionary diversity----------------#
#------------------------------------------------------------------------------------------#
#-------------------------5.a.Examining correlations among WorldClim the predictors-------


WG_pred_indices <- left_join(WG_all_indices[,-7] %>% 
                              mutate(site=rownames(WG_all_indices)),
                              WG_pred %>% mutate(site=rownames(WG_all_indices)),by="site") %>% 
                              rename(x=x.x,y=y.x) %>% dplyr::select(-site,-x.y,-y.y)

WG_pred_indices %>% filter(is.na(elevation))
#1 grid cell with no world clim value
# Lat       Long  SR      PD   PE   TILD   PE.all
# 18.45833 72.875 16 1572.85 1.76 315.49 1.759526


summary(WG_pred_indices)
head(WG_pred_indices)#1 NA
nrow(WG_pred_indices)#1248
names(WG_pred_indices)

WG_pred %>% filter(is.na(elevation))


corrplot(cor(WG_pred_indices[,c(2,8:28)]%>% filter(!is.na(elevation)) %>% rename(Latitude=y)), method="circle", addCoef.col ='black', 
         diag=FALSE, type = "lower",order='alphabet', number.cex = 0.3)
dev.off()

corrplot(cor(WG_pred_indices %>% dplyr::select(y, elevation, ann_prec, ann_mean_temp,  
                                               prec_season,temp_season,CWD) %>% 
               filter(!is.na(elevation)) %>% 
               rename(Latitude=y)),
               method="circle", addCoef.col ='black', 
               diag=FALSE,order='alphabet', number.cex = 0.5)

dev.off()

#-------------------------5.b.correlations of Key predictor with evolutionary diversity indices----------------------------------------

#TILD
TILD_plot <-WG_pred_indices %>% mutate(TILD=round(TILD,2)) %>% 
  ggplot(aes(y=ann_prec,x=CWD, colour=TILD)) +geom_point()+
  theme_bw()+ylab("Annual precipitation")+xlab("CWD")+ scale_colour_gradient(
  breaks=c(min(WG_pred_indices$TILD),max(WG_pred_indices$TILD)),
  low="#FFFFE5", high="#662506", name= "TILD")+theme_custom2+theme(legend.position = "none")+
  xlab("")+
  
  WG_pred_indices %>% mutate(TILD=round(TILD,2)) %>% 
  ggplot(aes(x=ann_prec,y=elevation, colour=TILD)) +geom_point()+
  theme_bw()+xlab("Annual precipitation")+ylab("Elevation")+ scale_colour_gradient(
  breaks=c(min(WG_pred_indices$TILD),max(WG_pred_indices$TILD)),
  low="#FFFFE5", high="#662506", name= "TILD")+ theme_custom2+
  theme(legend.position = "none")+ scale_x_continuous(expand = c(0.1,0.1))+    xlab("")+
  
  WG_pred_indices %>%  mutate(TILD=round(TILD,2)) %>% 
  ggplot(aes(x=elevation,y=CWD, colour=TILD)) +geom_point()+
  theme_bw()+xlab("Elevation")+ylab("CWD")+ scale_colour_gradient(
  breaks=c(round(min(WG_pred_indices$TILD),2),round(max(WG_pred_indices$TILD),2)),
  low="#FFFFE5", high="#662506", name= "TILD")+theme_custom2+    xlab("")
#PD    
PD_plot <-  WG_pred_indices %>% mutate(PD=round(PD,2)) %>% 
  ggplot(aes(y=ann_prec,x=CWD, colour=PD)) +geom_point()+
  theme_bw()+ylab("Annual precipitation")+xlab("CWD")+ scale_colour_gradient(
  breaks=c(min(WG_pred_indices$PD),max(WG_pred_indices$PD)),
  low="#FFFFE5", high="#662506", name= "PD")+theme_custom2+theme(legend.position = "none")+
  xlab("")+
  
  WG_pred_indices %>% mutate(PD=round(PD,2)) %>% 
  ggplot(aes(x=ann_prec,y=elevation, colour=PD)) +geom_point()+
  theme_bw()+xlab("Annual precipitation")+ylab("Elevation")+ scale_colour_gradient(
  breaks=c(min(WG_pred_indices$PD),max(WG_pred_indices$PD)),
  low="#FFFFE5", high="#662506", name= "PD")+ theme_custom2+
  scale_x_continuous(expand = c(0.1,0.1))+theme(legend.position = "none")+    xlab("")+
  
  WG_pred_indices %>%  mutate(PD=round(PD,2)) %>% 
  ggplot(aes(x=elevation,y=CWD, colour=PD)) +geom_point()+
  theme_bw()+xlab("Elevation")+ylab("CWD")+ scale_colour_gradient(
  breaks=c(round(min(WG_pred_indices$PD),2),round(max(WG_pred_indices$PD),2)),
  low="#FFFFE5", high="#662506", name= "PD")+theme_custom2+  xlab("")
#PE
PE_plot <-  WG_pred_indices %>% mutate(PE=round(PE,2)) %>% 
           ggplot(aes(y=ann_prec,x=CWD, colour=PE)) +geom_point()+
           theme_bw()+ylab("Annual precipitation")+xlab("CWD")+ scale_colour_gradient(
          breaks=c(min(WG_pred_indices$PE),max(WG_pred_indices$PE)),
          low="#FFFFE5", high="#662506", name= "PE")+theme_custom2+theme(legend.position = "none")+
  
  
  WG_pred_indices %>% mutate(PE=round(PE,2)) %>% 
  ggplot(aes(x=ann_prec,y=elevation, colour=PE)) +geom_point()+
  theme_bw()+xlab("Annual precipitation")+ylab("Elevation")+ scale_colour_gradient(
  breaks=c(min(WG_pred_indices$PE),max(WG_pred_indices$PE)),
  low="#FFFFE5", high="#662506", name= "PE")+ theme_custom2+
  scale_x_continuous(expand = c(0.1,0.1))+theme(legend.position = "none")+  
  
  WG_pred_indices %>%  mutate(PE=round(PE,2)) %>% 
  ggplot(aes(x=elevation,y=CWD, colour=PE)) +geom_point()+
  theme_bw()+xlab("Elevation")+ylab("CWD")+ scale_colour_gradient(
  breaks=c(round(min(WG_pred_indices$PE),2),round(max(WG_pred_indices$PE),2)),
  low="#FFFFE5", high="#662506", name= "PE")+theme_custom2

wrap_plots(TILD_plot,PD_plot,PE_plot,nrow=3)



dat <- WG_pred_indices %>% filter(!is.na(elevation))  %>% 
      dplyr::select(PD, CWD)

cor.test(dat$CWD,dat$PD) 
cor.test(dat$CWD,dat$PD,method=c("spearman")) 

rm(dat)

reg_cor <- function(div_index,preds)
{
  dat <- WG_pred_indices %>% filter(!is.na(elevation))  %>% 
    dplyr::select(div_index, preds)
  
  cor_pearson <- with(dat,cor.test(dat[[preds]],dat[[div_index]]))
  
  # cor_spearman <- with(dat,cor.test(dat[[preds]],dat[[div_index]],
  #                  method="spearman", conf.level = 0.95))
  
  estimates <- 
    cbind(div_index=div_index,pred=preds,
          cor_estimate=cor_pearson$estimate,
          cor_p=cor_pearson$p.value,
          cor_conf=paste0(round(cor_pearson$conf.int[1],2)," - " ,
                          round(cor_pearson$conf.int[2],2))
          #spearman
          # spearman_estimate=cor_spearman$estimate,
          # spear_p=cor_spearman$p.value
          # # spear_conf=paste0(round(cor_spearman$conf.int[1],2)," - " ,
          #                   round(cor_spearman$conf.int[2],2))
    ) %>% 
    as.data.frame() #%>%  mutate_at(c(3,4), as.numeric) %>% #mutate_at(c(3,4,6,7), as.numeric) %>% 
  #mutate_if(is.numeric, round, 4)
  
  
  return(estimates)
  
}

div_index <- c("TILD","PD", "PE")
preds <- c("elevation","ann_prec","CWD")


cor_est <- list()
for(j in 1:3){
  for(i in 1:length(preds))
  {
    cor_est[[paste0(j,i)]] <- reg_cor(div_index[j], preds[i])
  }
}

cor_est <- (do.call(rbind,cor_est))
cor_est %>% arrange(pred)
#write.xlsx(cor_est,paste0(final_sub_path,"/Final_plots/6_Dec_22/Tables/corr_evol_div_6Dec22.xlsx"))

###############################################################################
#-------------------------------6.a.Main figures-----------------------------------
###############################################################################
#-------------------------------Figure 1: Diversity indices-------------------
wrap_plots(
  plot_SR_WG+theme_custom+ylab("Latitude")+xlab("Longitude")+ggtitle('A')+
    ylim(8, max(WG_all_indices$y))+#clipping the WG at 19 as no predictions were done over that
    theme(legend.key.size = unit(0.3, "cm"),legend.position=c(0.3,0.25),
          axis.title = element_text( size=20,face="bold",   family="serif"),
          axis.text = element_text( size=20,colour = "black", family = "serif" ),
          legend.title = element_text(size=16,face = "bold", family = "serif"),
          legend.text = element_text( size=16, family = "serif"),
          plot.title = element_text(size = 20,face = "bold",  family = "serif",hjust = 0.5,vjust = 0.3)),
  
  
  plot_PD_WG+theme(legend.title = element_text(size=10,face = "bold", family = "serif"),
                   legend.text = element_text(size=10,family = "serif"))+ggtitle("")+
    ylim(8, max(WG_all_indices$y))+#clipping the WG at 19 as no predictions were done over that
    ylab("")+theme_void()+ggtitle('B')+theme(legend.key.size = unit(0.3, "cm"),
                                             legend.position=c(0.3,0.25),
                                             legend.title = element_text(size=16,face = "bold", family = "serif"),
                                             legend.text = element_text( size=16, family = "serif"),
                                             plot.title = element_text(size = 20,face = "bold",  
                                                                       family = "serif",hjust = 0.5)),
  
  plot_PE_350sp_WG+theme(legend.title = element_text(size=10,face = "bold", family = "serif"),
                         legend.text = element_text(size=10,family = "serif"))+ggtitle("")+
    ylim(8, max(WG_all_indices$y))+#clipping the WG at 19 as no predictions were done over that
    ylab("")+theme_void()+ggtitle('C')+theme(legend.key.size = unit(0.3, "cm"),
                                             legend.position=c(0.3,0.25),
                                             legend.title = element_text(size=16,face = "bold", family = "serif"),
                                             legend.text = element_text( size=16, family = "serif"),
                                             plot.title = element_text(size = 20,face = "bold", 
                                                                       family = "serif",hjust = 0.5)),
  
  
  plot_TILD_WG+theme(legend.title = element_text(size=10,face = "bold", family = "serif"),
                     legend.text = element_text(size=10,family = "serif"))+ggtitle("")+
    ylim(8, max(WG_all_indices$y))+#clipping the WG at 19 as no predictions were done over that
    ylab("")+theme_void()+ggtitle('D')+theme(legend.key.size = unit(0.3, "cm"),
                                             legend.position=c(0.3,0.25),
                                             legend.title = element_text(size=16,face = "bold", family = "serif"),
                                             legend.text = element_text( size=16, family = "serif"),
                                             plot.title = element_text(size = 20,face = "bold",  family = "serif",hjust = 0.5)),
  
  nrow = 1, ncol=4)

#ggsave(paste0(final_sub_path,"/Final_plots/6_Dec_22/1_Diversity_indices_clipped.png"), plot=last_plot(),units="in", 
width= 12, height = 8, dpi =600,bg="white")


# -------------------------------Figure 2a: Diversity at diff time depth-------------------
# PD_depth_lat_plot <-   wrap_plots(
#   plot_PD.depth[[15]]+ylab("Latitude")+xlab("Longitude")+theme_custom+
#     ylim(8, max(WG_all_indices$y))+#clipping the WG at 19 as no predictions were done over that
#     theme(legend.key.size = unit(0.3, "cm"),legend.position=c(0.3,0.25),
#           axis.title = element_text( size=20,face="bold",   family="serif"),
#           axis.text = element_text( size=20,colour = "black", family = "serif" ),
#           legend.title = element_text(size=16,face = "bold", family = "serif"),
#           legend.text = element_text( size=16, family = "serif")),
# 
#   plot_PD.depth[[13]]+ylim(8, max(WG_all_indices$y))+#clipping the WG at 19 as no predictions were done over that
#     theme(legend.key.size = unit(0.3, "cm"),
#           legend.position=c(0.3,0.25),legend.title = element_text(size=16,face = "bold", family = "serif"),
#           legend.text = element_text( size=16, family = "serif")),
#   plot_PD.depth[[10]]+ylim(8, max(WG_all_indices$y))+#clipping the WG at 19 as no predictions were done over that
#     theme(legend.key.size = unit(0.3, "cm"),
#           legend.position=c(0.3,0.25),legend.title = element_text(size=16,face = "bold", family = "serif"),
#           legend.text = element_text( size=16, family = "serif")),
#   plot_PD.depth[[7]]+ylim(8, max(WG_all_indices$y))+#clipping the WG at 19 as no predictions were done over that
#     theme(legend.key.size = unit(0.3, "cm"),
#           legend.position=c(0.3,0.25),legend.title = element_text(size=16,face = "bold", family = "serif"),
#           legend.text = element_text( size=16, family = "serif")),
#   plot_PD.depth[[4]]+ylim(8, max(WG_all_indices$y))+#clipping the WG at 19 as no predictions were done over that
#     theme(legend.key.size = unit(0.3, "cm"),
#           legend.position=c(0.3,0.25),legend.title = element_text(size=16,face = "bold", family = "serif"),
#           legend.text = element_text( size=16, family = "serif")),
#   plot_PD.depth[[2]]+ylim(8, max(WG_all_indices$y))+#clipping the WG at 19 as no predictions were done over that
#     theme(legend.key.size = unit(0.3, "cm"),
#           legend.position=c(0.3,0.25),legend.title = element_text(size=16,face = "bold", family = "serif"),
#           legend.text = element_text( size=16, family = "serif")),nrow=1)
# 
# 
# #ggsave(paste0(final_sub_path,"/Final_plots/6_Dec_22/2a_PD_time_slice_select_clipped.png"),
# plot=last_plot(),
# device = "png",width = 16, height = 8, units = "in",dpi=600)

#-------------------------------Figure 2a: LTT -------------------
devtools::install_github("aljrico/gameofthrones")
library(gameofthrones)

LTT_lat_plot <- LTT_lat.long %>% 
  ggplot(aes(x = reorder(age, desc(age)),y=log10(num), colour=lat.bin, group=lat.bin))+
  geom_line()+xlab("Age")+ylab(expression(paste(log[10],"(Num of Lineages)")))+
  scale_color_paletteer_c("gameofthrones::baratheon2", direction = -1,name = "Latitudnal \nbins") +
  geom_point( data= LTT_lat.long %>% filter(age==0),aes(x=age, y=log10(num)))+
  geom_text_repel(data=LTT_lat.long %>% filter(age==0),aes(x=age, y=log10(num),label=lat.bin),direction = c("both", "y", "x"),
                  size=8,max.overlaps = getOption("ggrepel.max.overlaps", default = 14))+
  scale_x_discrete(breaks=c(135,seq(130,10, -10),0))+  theme_classic(base_size = 30)+
  theme(legend.key.size = unit(1, "cm"),
        axis.title = element_text(size=30,face="bold", family="Times New Roman"),
        axis.text = element_text(size=30,colour = "black",family = "Times New Roman" ),
        legend.title = element_text( size=24,face = "bold", family = "Times New Roman"),
        legend.text = element_text(size=24,family = "Times New Roman"),
        legend.position=c(0.9,0.3))

LTT_lat_plot
#ggsave(paste0(final_sub_path,"/Final_plots/6_Dec_22/2b_LTT.png"), plot=last_plot(), units = "in",
width= 12, height = 10, dpi =600)

#-------------------------------Figure 2b: Lineage age trends-------------------
plot_age_lat <- function(age)#,label)
{  
  temp_age <- age_lat_distn(phylo_tree_final,WG.mat.reorder,
                            as.data.frame(cbind(grid=paste0(WG_pred$x,":",WG_pred$y),
                                                lat=as.numeric(WG_pred$y))),  age)
  #temp_age_out <- age_profile(temp_age,1)
  
  p <-  temp_age %>% ggplot(aes(x=1:length(unique_groups),y=mean_lat))+ylab("")+xlab("")+theme_bw()+
    geom_segment(aes(x=1:length(unique_groups),y=min_lat,xend=1:length(unique_groups),yend=max_lat))+theme_custom2+
    ggtitle(paste0(age, " Mya"))+theme_custom2+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                                                     plot.title = element_text(hjust = 0.5,family = 'Times New Roman' ))
  return(p)
  
}

#label <- c(letters[seq( from = 1, to = 6 )])

Lineage_age_lat_plot <- wrap_plots( mapply(plot_age_lat,rev(temp_age_list),#label=label,
                                           SIMPLIFY = F),nrow = 1)

Lineage_age_lat_plot + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(family = 'Times New Roman'))
#ggsave(paste0(final_sub_path,"/Final_plots/6_Dec_22/2c_Lineage_age_lat_rich.png"), plot=last_plot(),
units = "in", width= 22, height = 4, dpi =600)

#-------------------------------Figure 3: Correlates of evolutionary diversity-------------------
wrap_plots(TILD_plot,PD_plot,PE_plot,nrow=3)
#ggsave(paste0(final_sub_path,"/Final_plots/6_Dec_22/3_Environmental_correlates.png"), plot=last_plot(), units="in", dpi=600,
width=12, height=8)
#-------------------------------Figure 4a: Superorder contribution at each time period-------------------
superorder_age_PD_plot <- 
  wrap_plots(
    Order_lat_age_PD_final %>% filter(biog=="SWG") %>% group_by(Age,order) %>% 
      mutate(Age=as.factor(Age), prop_pd=as.numeric(prop_pd))  %>% 
      ggplot(aes( Age,order, fill= prop_pd)) + ggtitle("SWG 8-11 deg")+scale_x_discrete(limits = rev)+
      geom_tile()+  scale_fill_gradient(low="#FFFFE5", high="#662506", name= "Prop. PD") +theme_bw(base_size = 25)+ylab("")+
      theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),
            legend.text=element_text(size=15),legend.key.size = unit(1, "cm"),
            axis.text.y = element_text(size = 30), 
            axis.text = element_text(colour = "black"),text=element_text(family="Times")),
    
    Order_lat_age_PD_final %>% filter(biog=="CWG") %>% group_by(Age,order) %>% 
      mutate(Age=as.factor(Age), prop_pd=as.numeric(prop_pd))  %>% 
      ggplot(aes( Age,order, fill= prop_pd)) + ggtitle("CWG 12-15 deg")+scale_x_discrete(limits = rev)+
      geom_tile()+  scale_fill_gradient(low="#FFFFE5", high="#662506", name= "Prop. PD") +theme_bw(base_size =25)+ylab("")+
      theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),
            legend.text=element_text(size=15),legend.key.size = unit(1, "cm"),
            axis.text.y = element_blank(),axis.text = element_text(colour = "black"),text=element_text(family="Times")),
    
    Order_lat_age_PD_final %>% filter(biog=="NWG") %>% group_by(Age,order) %>% 
      mutate(Age=factor(Age), prop_pd=as.numeric(prop_pd)) %>% 
      ggplot(aes( Age,order, fill= prop_pd)) + ggtitle("NWG 16-19 deg")+scale_x_discrete(limits = rev)+
      geom_tile()+  scale_fill_gradient(low="#FFFFE5", high="#662506", name= "Prop. PD") +theme_bw(base_size = 25)+ylab("")+
      theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),
            legend.text=element_text(size=15),legend.key.size = unit(1, "cm"),
            axis.text.y = element_blank(),
            axis.text = element_text(colour = "black"),text=element_text(family="Times")))

superorder_age_PD_plot
#ggsave(paste0(final_sub_path,"/Final_plots/6_Dec_22/4a_Superorder_PD_biog.png"), plot=last_plot(),
units = "in", width= 26, height = 8, dpi =600)


#-------------------------------Figure 4b: Family level contribution to PD at each biog-------------------
Key_family_biog_PD_plot <- family_PD_biog_lat %>% group_by(biog) %>% mutate(tot_pd=sum(pd),) %>% 
  rename("family"=taxonomic_group) %>% 
  mutate(prop_pd=pd/tot_pd) %>% filter(family %in% dom_family$family) %>% 
  filter(pd!=0) %>% #Fabaceae is removed at NWG
  mutate(biog = factor(biog, levels=c("SWG","CWG","NWG"))) %>% 
  ggplot(aes(fill=family, y=prop_pd, x=biog)) +
  geom_bar(position="stack", stat="identity", colour="black")+theme_bw(base_size = 15)+coord_flip()+xlab("")+
  ylab("Prop. of PD ")+theme_custom +geom_text(aes(y = prop_pd, label = family),
                                               colour = "black",position=position_stack(vjust=0.5), stat="identity", size=6,
                                               fontface="bold",family="Times")+
  theme( axis.text = element_text(colour = "black", size=30),
         axis.title.x  = element_text(colour = "black", size=25),
         text=element_text(family="Times"))+theme(legend.position = "none",#"top",legend.direction = "horizontal",
         )+ scale_fill_manual(values=c(rep("#FFFFE5",10)))
  #scale_fill_manual(values=c(colorRampPalette(c("#996D50","#FFFFE5"))(10)))
#guides(fill = guide_legend(nrow = 2, title=NULL))

#colorRampPalette(c("#996D50","#FFFFE5"))(10)

#Adding number of families
family_PD_biog_lat_plot <- family_PD_biog_lat %>% group_by(biog) %>% mutate(tot_pd=sum(pd),) %>% 
  rename("family"=taxonomic_group) %>% group_by(biog) %>% filter(pd!=0) %>% 
  mutate(Num_fam=n_distinct(family))

family_PD_biog_lat_plot<- left_join(family_PD_biog_lat %>% group_by(biog) %>% mutate(tot_pd=sum(pd),) %>% 
                                      rename("family"=taxonomic_group) %>% 
                                      mutate(prop_pd=pd/tot_pd) %>% filter(family %in% dom_family$family) %>% 
                                      filter(pd!=0) %>% #Fabaceae is removed at NWG
                                      mutate(biog = factor(biog, levels=c("SWG","CWG","NWG"))),
                                    family_PD_biog_lat %>% group_by(biog) %>%   rename("family"=taxonomic_group) %>% 
                                      filter(pd!=0) %>% summarise(num=n_distinct(family))) %>% as.data.frame()




Key_family_biog_PD_plot <- Key_family_biog_PD_plot+ geom_text(data=family_PD_biog_lat_plot,
                                                              aes(label= paste("Num = ",num)),position="identity", 
                                                              stat="identity", size=5,y=0.59,fontface="bold",family="Times")


Key_family_biog_PD_plot + ylim(0,0.61)

rm(family_PD_biog_lat_plot)# used only for num familes

#ggsave(paste0(final_sub_path,"/Final_plots/6_Dec_22/4b_Dominant Family_PD_Biog2.png"), plot=last_plot(),units="in",
dpi=600,width=29, height=6)


#-------------------------------Table 1: Lineage turnover----------------------------------------------
lineage_turnover_final %>%   gt()%>% 
  cols_label(age=md("Age (Mya)"),
             num_lineages=md("Num. of lineages"),
             obs_turnover=md("Observed turnover"),
             mean_turnover=md("Expected Turnover")) %>% 
  tab_options(table.font.names = "Times New Roman")

library(openxlsx)

#write.xlsx(lineage_turnover_final,paste0(final_sub_path,"/Final_plots/6_Dec_22/Tables/lineage_turnover_6Dec22.xlsx"))


#-------------------------------Table 2: Lineage age trends----------------------------------------------
#Outputting the table at select time periods
temp_age_list <- c(10,30,60,90,120,135)

Age_lineage_table_select<- lapply(rev(temp_age_list),Age_lat_table)
Age_lineage_table_select <- as.data.frame(Reduce(rbind,Age_lineage_table_select))


Age_lineage_table_select %>%   gt()%>% 
  fmt_number(
    columns = peak_lat,decimal=1,
    suffixing = TRUE) %>% cols_label(unique_groups=md("Total lineages"),
                                     Age=md("Age (Mya)"), lt_13=md("Num. lineages only < 13 deg"),
                                     gt_13=md("Num. lineages only > 13 deg"),peak_lat=md("Latitude of peak richness"),
                                     lineage_peak_lat=md("Number of lineages at peak richness"))




###############################################################################
#-------------------------------6b.Supplemetary Figures----------------------------
###############################################################################


#-------------------------------Figure S2---------------------------------------
rownames(biodiv_final_all) <- paste0(biodiv_final_all$x,":",biodiv_final_all$y)

biodiv_final_all %>% mutate(site=rownames(biodiv_final_all)) %>% 
  pivot_longer(cols=Acrocarpus_fraxinifolius:Agasthiyamalai_pauciflora, names_to = "sp", values_to = "num") %>% 
  filter(num>0) %>% group_by(site) %>% summarise(sum_num=sum(num)) %>% filter(sum_num < 10) %>%
  separate(site, c("long", "lat"), sep = ":") %>% 
  mutate(across(c(long, lat), parse_number)) %>% ggplot()+geom_tile(aes(x=long , y=lat, fill=sum_num)) +
  geom_polygon(data=wg, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + 
  theme_bw(base_size = 15)+coord_fixed()+ scale_fill_paletteer_c("grDevices::YlOrRd", direction = -1,#scale_fill_gradientn(colours=c("#FFFFFF",  "#FEC44F", "#FE9929", "#EC7014", "#CC4C02","#993404", "#662506"),
                                                                 breaks=c(1,9), "SR")+ylab("Latitude")+xlab("Longitude")+theme_custom

# ggsave(paste0(final_sub_path,"/Final_plots/6_Dec_22/Supple/S2_Grids_with_lt_9species_6Dec22.png"),
plot=last_plot(), units = "in", width= 8, height = 8, dpi =600)

#-------------------------------Figure S3---------------------------------------
#png(filename=paste0(final_sub_path,"/Final_plots/6_Dec_22/Supple/S3a_All_pred_correlations_10km_6Dec22.png"),
width=12, height = 10, units = "in",res=600)

par(cex=1.5)
corrplot(cor(WG_pred_indices[,c(2,8:28)]%>% filter(!is.na(elevation)) %>% rename(latitude=y)), method="circle", addCoef.col ='black', 
         diag=FALSE, type = "lower",order='alphabet', number.cex = 0.3)
dev.off()

par(cex=2)

#png(filename=paste0(final_sub_path,"/Final_plots/6_Dec_22/Supple/S3b_Select_pred_correlations_10km_6Dec22.png"),width=12, height = 8, units = "in",res=600)

corrplot(cor(WG_pred_indices %>% dplyr::select(y, elevation, ann_prec, ann_mean_temp,  
                                               prec_season,temp_season,CWD) %>% 
               filter(!is.na(elevation)) %>% 
               rename(Latitude=y)),
         method="circle", addCoef.col ='black', 
         diag=FALSE,order='alphabet', number.cex = 1.5, cl.cex = 1.5, tl.cex = 1.5)

dev.off()

#-------------------------------Figure S4---------------------------------------

WG_all_indices %>% mutate( biog=factor(biog , levels=c("SWG", "CWG", "NWG"))) %>%  
  ggplot(aes(y=SR, x=biog))+geom_violin()+
  geom_boxplot(width=0.08,notch=TRUE)+theme_bw(base_size = 12)+
  xlab("")+ggtitle("A")+
  WG_all_indices %>% mutate( biog=factor(biog , levels=c("SWG", "CWG", "NWG"))) %>% 
  ggplot(aes(y=PD, x=biog))+geom_violin()+
  geom_boxplot(width=0.08,notch=TRUE)+theme_bw(base_size = 12)+
  xlab("")+ggtitle("B")+
  WG_all_indices %>% mutate( biog=factor(biog , levels=c("SWG", "CWG", "NWG"))) %>% 
  ggplot(aes(y=PE, x=biog))+geom_violin()+geom_boxplot(width=0.08,notch=TRUE)+
  theme_bw(base_size = 12)+xlab("")+ggtitle("C")+
  WG_all_indices %>% mutate( biog=factor(biog , levels=c("SWG", "CWG", "NWG"))) %>%
  ggplot(aes(y=TILD, x=biog))+geom_violin()+geom_boxplot(width=0.08,notch=TRUE)+
  theme_bw(base_size = 12)+xlab("")+ggtitle("D")


# ggsave("Final_plots/6_Dec_22/Supple/S4_Indices_biog.png", plot=last_plot(), units="in", dpi=600,
width=12, height=8)

#-------------------------------Figure S5---------------------------------------
#png(filename="Final_plots/6_Dec_22/Supple/S5_Correlations_6Dec22.png",width=8, height = 8, units = "in",res=600)
pairs.panels( WG_all_indices[,c(-1,-2,-7)] ,
              method = "pearson", # correlation method
              hist.col = "#00AFBB",
              density = TRUE,  # show density plots
              ellipses = F # show correlation ellipses
)

dev.off()



#-------------------------------Table S1---------------------------------------

Age_lineage_table_full <- lapply(rev(age_list),Age_lat_table)
Age_lineage_table_full <- as.data.frame(Reduce(rbind,Age_lineage_table_full))


#Outputting the table at all time periods
Age_lineage_table_full %>%   gt()%>% 
  fmt_number(
    columns = peak_lat,decimal=1,
    suffixing = TRUE) %>% cols_label(unique_groups=md("Total lineages"),
                                     Age=md("Age (Mya)"), lt_13=md("No. lineages only < 13 deg"),
                                     gt_13=md("No. lineages only > 13 deg"),peak_lat=md("Latitude of peak richness"),
                                     lineage_peak_lat=md("Number of lineages at peak richness"))


#write.xlsx(Age_lineage_table_full,paste0(final_sub_path,"/Final_plots/6_Dec_22/Tables/lineage_age_lat_6Dec22.xlsx"))

#-------------------------------Table S2---------------------------------------
All_sp_checklist.node %>% group_by(order,family) %>% 
  mutate(Endemic=as.numeric(Endemic))%>% 
  summarise(num_sp=n_distinct(label),
            endem=sum(Endemic,na.rm=T)) %>% arrange(desc(num_sp)) %>%
  mutate(prop_end=round(endem/num_sp,2)) %>% 
  arrange(-num_sp,-prop_end) %>% as.data.frame()


#-------------------------------Table S3---------------------------------------
All_sp_checklist.node %>% group_by(order,family, label) %>% arrange(order) %>% 
  dplyr::select(order,family, genus, label, Endemic)





