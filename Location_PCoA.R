## Analyzing 2019 location data

#load libraries
library(phyloseq)
library(ggplot2)
library(ape)
library(vegan)
library(ggpubr)
library(agricolae)
library(ggsignif)
library(dplyr)
library(microbiome)
library(tibble)
library(ggforce)
library(tidyverse) 
library(ggrepel) # optional - nice to help keep labels on ggplots from running into each other


#### ====================================================================== ####

#### Load 16S data ####
#### ====================================================================== ####

#import feature table, taxonomy table, and metadata
feat_table = read.csv("mergedredo/12500/split_feat_table.csv", sep = ",", row.names =1)
feat_table = as.matrix(feat_table)

taxonomy = read.csv("mergedredo/12500/split_tax.csv", sep = ",", row.names = 1)
taxonomy <- replace(taxonomy, taxonomy == "", NA) #make any empty strings in taxonomy read NA
taxonomy = as.matrix(taxonomy)

metadata_full = read.table("mergedredo/updated_metadata.tsv", sep = "\t", row.names = 1, header = TRUE)
metadata_full$year<-as.factor(metadata_full$year)


#import as phyloseq objects
OTU = otu_table(feat_table, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy)
META = sample_data(metadata_full)
#tree was already imported as phyloseq object

#check that your OTU names are consistent across objects
head(taxa_names(TAX))
head(taxa_names(OTU))
head(taxa_names(phy_tree))

#make sure files have same sample names
sample_names(OTU)
sample_names(META)

#create phyloseq object
physeq_16s = phyloseq(OTU, TAX, META) 

# Removing contaminants found in water controls
ps<-subset_taxa(physeq_16s, Family != "FamilyNA")
ps<-subset_taxa(ps, Family != "D_4__Streptomycetaceae")
ps<-subset_taxa(ps, Family != "D_4__Corynebacteriaceae")
ps<-subset_taxa(ps, Family != "D_4__Enterobacteriaceae")
ps<-subset_taxa(ps, Family != "D_4__Bacillaceae")
ps<-subset_taxa(ps, Family != "D_4__Staphylococcaceae")
ps<-subset_taxa(ps, Family != "D_4__Acidobacteriaceae (Subgroup 1)")
ps<-subset_taxa(ps, Family != "D_4__Pseudomonadaceae")
ps<-subset_taxa(ps, Family != "D_4__0319-6G20")
ps<-subset_taxa(ps, Family != "D_4__Bacteroidaceae")
ps<-subset_taxa(ps, Family != "D_4__Intrasporangiaceae")
ps<-subset_taxa(ps, Family != "D_4__Lactobacillaceae")
ps<-subset_taxa(ps, Family != "D_4__Burkholderiaceae")
ps<-subset_taxa(ps, Family != "D_4__Family XI")
ps<-subset_taxa(ps, Class != "D_2__Chloroplast")
ps<-subset_taxa(ps, Family != "D_4__Mitochondria")


###looking at location differences - subset into just 2019 wild data
loc_ps<- subset_samples(physeq_16s, year.plot == "2019/Plot1")
loc_ps<- subset_samples(loc_ps, site_type == "wild")



#### ====================================================================== ####

#### Load ITS data ####
#### ====================================================================== ####

#import feature table, taxonomy, and metadata
feat_tableits = read.csv("QDM011/merged/its/split_feat_table.csv", sep = ",", row.names =1, check.names = TRUE)
feat_tableits = as.matrix(feat_tableits)

taxonomyits = read.csv("QDM011/merged/its/split_tax.csv", sep = ",", row.names = 1)
taxonomyits <- replace(taxonomyits, taxonomyits == "", NA)
taxits<- taxonomyits %>%
  mutate(Class = ifelse(is.na(Class), "ClassNA", Class))
taxonomyits = as.matrix(taxits)

metadata_fullits = read.table("QDM011/merged/its/merged_itsmetadata.tsv", sep = "\t", row.names =1, header = TRUE)
metadata_fullits$year.plot<-as.factor(metadata_fullits$year.plot)


#import as phyloseq objects
OTUits = otu_table(feat_tableits, taxa_are_rows = TRUE)
TAXits = tax_table(taxonomyits)
METAits = sample_data(metadata_fullits)

#create phyloseq object
physeq_its = phyloseq(OTUits, TAXits, METAits, phy_treeits)


###looking at location differences - subset by 2019 wild data
locits_ps<- subset_samples(physeq_its, year.plot == "2019/Plot1")
locits_ps<- subset_samples(locits_ps, site_type == "wild")

#create new metadata using subsetted data
loc_metadata <- data.frame(sample_data(loc_ps)) 
locits_metadata <- data.frame(sample_data(locits_ps)) 

# Extract the taxonomic information - this will be useful for plotting later
loc_tax_info <- data.frame(tax_table(loc_ps))
locits_tax_info <- data.frame(tax_table(locits_ps))

# Extract OTU table and convert it out of phyloseq format
loc_otu_df <- data.frame(otu_table(loc_ps), check.names = FALSE) 
locits_otu_df <- data.frame(otu_table(locits_ps), check.names = FALSE) 



# Rename the OTUs
# I don't like my OTUs to only be named numbers, so I'll append the word "OTU" in front of them
# this makes it easier so R doesn't mistake them for numbers

# Rename them in the otu table
rownames(loc_otu_df) <- paste0("OTU_", rownames(loc_otu_df))
rownames(locits_otu_df) <- paste0("OTU_", rownames(locits_otu_df))

# Rename them in the taxonomy table
rownames(loc_tax_info) <- paste0("OTU_", rownames(loc_tax_info))
rownames(locits_tax_info) <- paste0("OTU_", rownames(locits_tax_info))


# Transform otu table so species are columns and rows are samples - this is the default that 
# vegan and most other non-microbial ecological packages expects
loc_otu_df_t <- t(loc_otu_df) # note: this converts from dataframe into matrix format
locits_otu_df_t <- t(locits_otu_df)


#### ====================================================================== ####

#### Step 2: Compute distance matrices  ####
# different distance metrics say slightly different things about the community,

#### ====================================================================== ####

# Bray-curtis dissimilarities
loc_bc <- vegdist(x = loc_otu_df_t, method = "bray")
locits_bc <- vegdist(x = locits_otu_df_t, method = "bray")



#### ====================================================================== ####


#### Step 3: Compute ordination #####

#### ====================================================================== ####

# PCOA oridination, with bray-curtis
loc_bc.pcoa <- pcoa(loc_bc) # this is from the package ape; 
loc_bc.pct_ex <- round((loc_bc.pcoa$values$Relative_eig) * 100, 1) # this line computes the variation on each axis

locits_bc.pcoa <- pcoa(locits_bc) # this is from the package ape; 
locits_bc.pct_ex <- round((locits_bc.pcoa$values$Relative_eig) * 100, 1)


#### ====================================================================== ####

#### Step 5: Prepare plotting dataframes #####
# here we join our sample metadata with the points that were generated by the ordination calculations above
# For simplicity we'll actually just create 1 dataframe with all of the different sets of ordination points
# but you could also make them separately
#### ====================================================================== ####


loc_ord_plot.df <- loc_metadata %>%
  # First join the pcoa data
  left_join(loc_bc.pcoa$vectors[,1:2] %>% data.frame() %>%
              rownames_to_column(var = "SampleID"), # here we make a SampleID column before joining the dataframes
            by = "SampleID") %>%
  dplyr::rename(PCOA1_BC = Axis.1, PCOA2_BC = Axis.2)


locits_ord_plot.df <- locits_metadata %>%
  # First join the pcoa data
  left_join(locits_bc.pcoa$vectors[,1:2] %>% data.frame() %>%
              rownames_to_column(var = "SampleID"), # here we make a SampleID column before joining the dataframes
            by = "SampleID") %>%
  dplyr::rename(PCOA1_BC = Axis.1, PCOA2_BC = Axis.2)



#### ====================================================================== ####

#### Compute Arrows ####
#### ====================================================================== ####

# Envfit  - more robust solution works for nmds or PCOA
# this wrapper modifies OTU table to go into envfit, but the same idea can be applied to 
# any numeric variable or a factor that can be converted into a numeric variable. Example: pH, SoilTexture
env_fit_wrapper <- function(ordination, otu_tab, ordination_type = "NMDS", tax_cutoff = NULL) {
  
  # ordination = ordination output
  # ordination_type = either NMDS or PCOA
  # otu_tab = table of organisms abundances in each sample; can be full otu table or subset; columns are organisms
  # tax_cutoff = A number, the abundance cutoff for OTU table, the top xxx most abundant taxa in table will be selected
  
  # Prepare data
  if(ordination_type == "NMDS") {
    ord_axes <- ordination$points
  }
  
  if(ordination_type == "PCOA") {
    ord_axes <- ordination$vectors
  }
  
  # Filter OTU table
  if(!is.null(tax_cutoff)) {
    top_tax <- sort(colSums(otu_tab, na.rm = T), decreasing = T)[1:tax_cutoff]
    names(top_tax)
    otu_filt <- otu_tab[,names(top_tax)]
  } else {
    otu_filt <- otu_tab
  }
  
  # Run envfit
  fit <- envfit(ord = ord_axes, env = otu_filt, perm = 999, na.rm = TRUE)
  
  # adjust p-values for many tests
  pvals.adj <- p.adjust(fit$vectors$pvals, method = "fdr") #false discovery rate adjustment
  fit$vectors$pval.adj <- pvals.adj
  
  # Use scores to transforms the axes values so the legth of the arrow is propotional to R2
  # Further multiply it for plotting scale with 'OrdiArrowMul"
  scaled_arrow_heads <- data.frame(scores(fit, "vectors"))*ordiArrowMul(fit) 
  
  # Do some data cleaning up for ease of plotting
  arrow_output <- data.frame(scaled_arrow_heads,
                             fit$vectors$arrows,
                             r2 = fit$vectors$r,
                             pvals = fit$vectors$pvals, 
                             pval.adj = fit$vectors$pval.adj) %>%
    rename(Axis.1 = 1, Axis.2 = 2, 
           Axis.1.unscaled = 3, Axis.2.unscaled = 4) %>%
    rownames_to_column(var = "FactorLabel")
  return(arrow_output)
  
} 



loc_bc.pcoa.arrows <- env_fit_wrapper(ordination = loc_bc.pcoa,
                                       ordination_type = "PCOA",
                                       otu_tab = loc_otu_df_t, 
                                       tax_cutoff = 50) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but nice
  left_join(loc_tax_info %>% 
              rownames_to_column(var = "OTU_ID"),
            by = c("FactorLabel" = "OTU_ID"))


locits_bc.pcoa.arrows <- env_fit_wrapper(ordination = locits_bc.pcoa,
                                      ordination_type = "PCOA",
                                      otu_tab = locits_otu_df_t, 
                                      tax_cutoff = 50) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but nice
  left_join(locits_tax_info %>% 
              rownames_to_column(var = "OTU_ID"),
            by = c("FactorLabel" = "OTU_ID"))


#### ====================================================================== ####

#### Step 5: Plot data #####
#### ====================================================================== ####


# PCOA Bray-Curtis
ggplot(data = loc_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = location)) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", loc_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", loc_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = loc_bc.pcoa.arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text(data = loc_bc.pcoa.arrows,
            mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel)) +
  # Alternatively
  geom_text_repel(data = loc_bc.pcoa.arrows,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel))


#### ====================================================================== ####

# Plot by taxonomic group
#### ====================================================================== ####

# Summarize OTU table by class
ClassTab <- t(loc_otu_df_t) %>% 
  data.frame() %>%
  rownames_to_column(var = "OTU_ID") %>%
  left_join(loc_tax_info %>% rownames_to_column(var = "OTU_ID"),
            by = "OTU_ID") %>%
  mutate(Class = ifelse(is.na(Class), "ClasssNA", Class)) %>% # fix NA classes
  group_by(Class) %>%
  summarize(across(BYUL1B1:VTWin32, ~sum(.))) %>% # sum class abundances in each sample column
  column_to_rownames(var = "Class") %>%
  t() # retransform data


# Summarize OTU table by class - ITS
itsClassTab <- t(locits_otu_df_t) %>% 
  data.frame() %>%
  rownames_to_column(var = "OTU_ID") %>%
  left_join(locits_tax_info %>% rownames_to_column(var = "OTU_ID"),
            by = "OTU_ID") %>%
  mutate(Class = ifelse(is.na(Class), "ClasssNA", Class)) %>% # fix NA classes
  group_by(Class) %>%
  summarize(across(BYUL1B2:VTStarr2, ~sum(.))) %>% # sum class abundances in each sample column
  column_to_rownames(var = "Class") %>%
  t() # retransform data




loc_bc.pcoa.class.arrows <- env_fit_wrapper(ordination = loc_bc.pcoa,
                                             ordination_type = "PCOA",
                                             otu_tab = ClassTab, 
                                             tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but  nice
  left_join(loc_tax_info %>% 
              rownames_to_column(var = "OTU_ID") %>%
              select(Phylum:Class) %>% distinct(),
            by = c("FactorLabel" = "Class"))



locits_bc.pcoa.class.arrows <- env_fit_wrapper(ordination = locits_bc.pcoa,
                                            ordination_type = "PCOA",
                                            otu_tab = itsClassTab, 
                                            tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but  nice
  left_join(locits_tax_info %>% 
              rownames_to_column(var = "OTU_ID") %>%
              select(Phylum:Class) %>% distinct(),
            by = c("FactorLabel" = "Class"))




# Plot Arrows by Class
# PCOA Bray-Curtis
ggplot(data = loc_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = texture, shape = location)) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", loc_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", loc_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = loc_bc.pcoa.class.arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text_repel(data = loc_bc.pcoa.class.arrows, max.overlaps = 20,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel))



ggplot(data = locits_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = texture, shape = location)) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", locits_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", locits_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = locits_bc.pcoa.class.arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text_repel(data = locits_bc.pcoa.class.arrows, max.overlaps = 20,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel))



#### ====================================================================== ####

# Plot an environmental factor 
#### ====================================================================== ####


env_dat<- loc_metadata[, 15:ncol(loc_metadata)] #remove first 14 columns
env_dat1 <- subset(env_dat, select = -site) #remove the site column
env_dat1$ph<-as.numeric(env_dat1$ph)
env_dat1$calcium<-as.numeric(env_dat1$calcium)
env_dat1$magnesium<-as.numeric(env_dat1$magnesium)
env_dat1$potassium<-as.numeric(env_dat1$potassium)
env_dat1$phosphorus<-as.numeric(env_dat1$phosphorus)
env_dat1$lead<-as.numeric(env_dat1$lead)
env_dat1$om<-as.numeric(env_dat1$om)
env_dat1$elevation<-as.numeric(env_dat1$elevation)
env_dat1$av_temp<-as.numeric(env_dat1$av_temp)
env_dat1$av_precip<-as.numeric(env_dat1$av_precip)


envits_dat<- locits_metadata[, 15:ncol(locits_metadata)] #remove first 14 columns
envits_dat1 <- subset(envits_dat, select = -site) #remove the site column
envits_dat1$ph<-as.numeric(envits_dat1$ph)
envits_dat1$calcium<-as.numeric(envits_dat1$calcium)
envits_dat1$magnesium<-as.numeric(envits_dat1$magnesium)
envits_dat1$potassium<-as.numeric(envits_dat1$potassium)
envits_dat1$phosphorus<-as.numeric(envits_dat1$phosphorus)
envits_dat1$lead<-as.numeric(envits_dat1$lead)
envits_dat1$om<-as.numeric(envits_dat1$om)
envits_dat1$elevation<-as.numeric(envits_dat1$elevation)
envits_dat1$av_temp<-as.numeric(envits_dat1$av_temp)
envits_dat1$av_precip<-as.numeric(envits_dat1$av_precip)




fit <- envfit(ord = loc_bc.pcoa$vectors, # we're using the Bray-curtis pcoa in this instance 
              env = env_dat1, perm = 999, na.rm = TRUE)

# adjust p-values for many tests
pvals.adj <- p.adjust(fit$vectors$pvals, method = "fdr")
fit$vectors$pval.adj <- pvals.adj

fitits <- envfit(ord = locits_bc.pcoa$vectors, # we're using the Bray-curtis pcoa in this instance 
              env = envits_dat1, perm = 999, na.rm = TRUE)

# adjust p-values for many tests
pvals.adjits <- p.adjust(fitits$vectors$pvals, method = "fdr")
fitits$vectors$pval.adj <- pvals.adjits


# While numeric columns can be represented by arrows, factor levels should be represented by centroids

fit_cont <- (scores(fit, "vectors") * ordiArrowMul(fit)) %>% 
  data.frame() %>%
  rownames_to_column(var = "FactorLabel")

fit_cat <- (scores(fit, "factors") * ordiArrowMul(fit)) %>%
  data.frame() %>%
  rownames_to_column(var = "FactorLabel")

fitits_cont <- (scores(fitits, "vectors") * ordiArrowMul(fitits)) %>% 
  data.frame() %>%
  rownames_to_column(var = "FactorLabel")

fitits_cat <- (scores(fitits, "factors") * ordiArrowMul(fitits)) %>%
  data.frame() %>%
  rownames_to_column(var = "FactorLabel")


# PCOA Bray-Curtis
# PCOA Bray-Curtis
bac<-ggplot(data = loc_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = texture, shape = location), size = 3) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", loc_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", loc_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = fit_cont,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text_repel(data = fit_cont,
                  mapping = aes(x = Axis.1, y = Axis.2, 
                                label = FactorLabel, size = 45)) +
  theme(axis.title.x = element_text(color="black", size = 22)) +
  theme(axis.title.y = element_text(color="black", size =22)) +
  theme(axis.text.x = element_text(size =20, color = "black")) +
  theme(axis.text.y = element_text(size =20, color = "black")) +
  theme(legend.text = element_text(size =20)) +
  theme(legend.title = element_blank()) +
  ggtitle("Bacteria") +
  guides(size = "none") +
  theme(plot.title = element_text(color = "black", size = 24, face = "bold"))




fun<-ggplot(data = locits_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = texture, shape = location), size = 3) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", loc_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", loc_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = fitits_cont,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text_repel(data = fitits_cont,
                  mapping = aes(x = Axis.1, y = Axis.2, 
                                label = FactorLabel, size = 45)) +
  theme(axis.title.x = element_text(color="black", size = 22)) +
  theme(axis.title.y = element_text(color="black", size =22)) +
  theme(axis.text.x = element_text(size =20, color = "black")) +
  theme(axis.text.y = element_text(size =20, color = "black")) +
  theme(legend.text = element_text(size =20)) +
  theme(legend.title = element_blank()) +
  ggtitle("Fungi") +
  guides(size = "none") +
  theme(plot.title = element_text(color = "black", size = 24, face = "bold"))


ggarrange(bac, fun, nrow =1 , ncol =2, common.legend = TRUE, legend = "right")
  
