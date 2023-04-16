##Analyzing both 2019 and 2020 data

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
library(GUniFrac) 
library(ggrepel)


#### ====================================================================== ####

#### Load 16S data ####
#### ====================================================================== ####


#read in feature table, taxonomy, and metadata
feat_table = read.csv("mergedredo/12500/split_feat_table.csv", sep = ",", row.names =1)
feat_table = as.matrix(feat_table)

taxonomy = read.csv("mergedredo/12500/split_tax.csv", sep = ",", row.names = 1)
taxonomy <- replace(taxonomy, taxonomy == "", NA) # replace empty strings with NA, if Family is not assigned, label it "FamilyNA"
tax<- taxonomy %>%
  mutate(Family = ifelse(is.na(Family), "FamilyNA", Family))
taxonomy = as.matrix(tax)

metadata_full = read.table("mergedredo/updated_metadata.tsv", sep = "\t", row.names = 1, header = TRUE)
metadata_full$year<-as.factor(metadata_full$year)



#import as phyloseq objects
OTU = otu_table(feat_table, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy)
META = sample_data(metadata_full)

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


#Subset phyloseq object to just look at Woodman agricultural samples
wo_ps<- subset_samples(ps, location == "Durham_NH")
wo_ps<- subset_samples(wo_ps, site_type == "agricultural")

#Subset phyloseq object to just look at Vermont samples
vt_ps<- subset_samples(ps, location == "Vermont")



#### ====================================================================== ####

#### Load ITS data ####
#### ====================================================================== ####

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
physeq_its = phyloseq(OTUits, TAXits, METAits)


##looking at just woodman samples
wo_psits<- subset_samples(physeq_its, location == "Durham_NH")
wo_psits<- subset_samples(wo_psits, site_type == "agricultural")

#Making the ag vs wild its dataset
duits_ps<-subset_samples(physeq_its, location == "Durham_NH")
ryeits_ps<-subset_samples(physeq_its, location == "Rye_NH")
nhits_ps = merge_phyloseq(duits_ps, ryeits_ps)

#Splitting into 2019 and 2020
funsite19 <- subset_samples(nhits_ps, year.plot == "2019/Plot1")
funsite20 <- subset_samples(nhits_ps, year.plot == "2020/Plot2")

##Subsetting out just VT samples
vt_ps_its <- subset_samples(physeq_its, location == "Vermont")




#### ====================================================================== ####

#### PcoA's, PERMANOVA in Phyloseq ####
#### ====================================================================== ####





# 16s woodman by year
ordu = ordinate(ps, "PCoA", "bray") #other options like "NMDS" and "unifrac"
plot_ordination(ps, ordu, color = "soil_type", shape = "year") +
  geom_point(size = 4) +
  theme(axis.title.x = element_text(color="black", size = 20)) +
  theme(axis.title.y = element_text(color="black", size =20)) +
  theme(axis.text.x = element_text(size =12, color = "black")) +
  theme(axis.text.y = element_text(size =12, color = "black")) +
  theme(legend.text = element_text(size =14)) +
  theme(legend.title = element_text(size =14, face = "bold"))

# ITS woodman by year
ordu = ordinate(wo_psits, "PCoA", "bray") #other options like "NMDS" and "unifrac"
plot_ordination(wo_psits, ordu, color = "host", shape = "year") +
  geom_point(size = 4) +
  theme(axis.title.x = element_text(color="black", size = 20)) +
  theme(axis.title.y = element_text(color="black", size =20)) +
  theme(axis.text.x = element_text(size =12, color = "black")) +
  theme(axis.text.y = element_text(size =12, color = "black")) +
  theme(legend.text = element_text(size =14)) +
  theme(legend.title = element_text(size =14, face = "bold"))


###### PERMANOVA - Analyzing Woodman data separately by year
##permanova on host
ps19<-subset_samples(wo_ps, year.plot == "2019/Plot1")

ps19_meta <- as(sample_data(ps19), "data.frame")
adonis2(distance(ps19, method="bray") ~ host_type*soil_type, data = ps19_meta) 
#host_type p = 0.46
#soil_type p = 0.01
#host_type:soil_type p = 0.87




ps20<-subset_samples(ps, year.plot == "2020/Plot2")

ps20_meta <- as(sample_data(ps20), "data.frame")
adonis2(distance(ps20, method="bray") ~ host_type*soil_type, data = ps20_meta) 
#host_type p = 0.145
#soil_type p = 0.001
#host_type:soil_type p = 0.678








###### Site type PERMANOVA
#create a subset of ps phyloseq object
du_ps<-subset_samples(ps, location == "Durham_NH")
rye_ps<-subset_samples(ps, location == "Rye_NH")
nh_ps = merge_phyloseq(du_ps, rye_ps)


site_meta <- as(sample_data(nh_ps), "data.frame")
adonis2(distance(nh_ps, method="bray") ~ site_type*year.plot, data = site_meta) 
#site type x year interaction p = 0.001

site19<-subset_samples(nh_ps, year.plot == "2019/Plot1")
site19_meta <- as(sample_data(site19), "data.frame")
adonis2(distance(site19, method="bray") ~ site_type*texture, data = site19_meta) 
#site_type p = 0.001
#texture p = 0.001


site20<-subset_samples(nh_ps, year.plot == "2020/Plot2")

site20_meta <- as(sample_data(site20), "data.frame")
adonis2(distance(site20, method="bray") ~ site_type*texture, data = site20_meta) 
#site_type p = 0.001
#texture p = 0.002




## CAP Plots for Woodman Fungi
wo_psits20<- subset_samples(wo_psits, year.plot == "2020/Plot2")
metaits20 <- as(sample_data(wo_psits20), "data.frame")

wo_psits19<- subset_samples(wo_psits, year.plot == "2019/Plot1")
metaits19 <- as(sample_data(wo_psits19), "data.frame")


ordutype20 = ordinate(wo_psits20, "CAP", "bray", formula = OTU ~ host_type)
cap_hosttype_20<-plot_ordination(wo_psits20, ordutype20, color = "host_type") +
  geom_point(size = 3) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size = 22)) +
  theme(axis.title.y = element_text(color="black", size =22)) +
  theme(axis.text.x = element_text(size =20, color = "black")) +
  theme(axis.text.y = element_text(size =20, color = "black")) +
  theme(legend.text = element_text(size =20)) +
  theme(legend.title = element_text(size =20)) +
  labs(color = "Host type")+
  ggtitle("2020 Plot") +
  theme(plot.title = element_text(size =24, face = "bold"))



ordutype19 = ordinate(wo_psits19, "CAP", "bray", formula = OTU ~ host_type)
cap_hosttype_19<-plot_ordination(wo_psits19, ordutype19, color = "host_type") +
  geom_point(size = 3) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size = 22)) +
  theme(axis.title.y = element_text(color="black", size =22)) +
  theme(axis.text.x = element_text(size =20, color = "black")) +
  theme(axis.text.y = element_text(size =20, color = "black")) +
  theme(legend.text = element_text(size =20)) +
  theme(legend.title = element_text(size =20, face = "bold")) +
  labs(color = "Host type") +
  ggtitle("2019 Plot") +
  theme(plot.title = element_text(size =24, face = "bold")) 

ggarrange(cap_hosttype_19, cap_hosttype_20, nrow =1, ncol =2, 
          common.legend = TRUE, legend = "right")




## Plotting just vermont samples by pH
ordu_vt_16s = ordinate(vt_ps, "PCoA", "bray")
vt_16s_meta <- as(sample_data(ordu_vt_16s), "data.frame")
vt_16s<-plot_ordination(vt_ps, ordu_vt_16s, color = "ph", shape = "texture") +
  geom_point(size = 3) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size = 22)) +
  theme(axis.title.y = element_text(color="black", size =22)) +
  theme(axis.text.x = element_text(size =20, color = "black")) +
  theme(axis.text.y = element_text(size =20, color = "black")) +
  theme(legend.text = element_text(size =20)) +
  theme(legend.title = element_text(size =20)) +
  labs(color = "pH") +
  ggtitle("Bacteria") +
  theme(plot.title = element_text(size =24, face = "bold")) 


ordu_vt_its = ordinate(vt_ps_its, "PCoA", "bray")
vt_its<-plot_ordination(vt_ps_its, ordu_vt_its, color = "om", shape = "texture") +
  geom_point(size = 3) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size = 22)) +
  theme(axis.title.y = element_text(color="black", size =22)) +
  theme(axis.text.x = element_text(size =20, color = "black")) +
  theme(axis.text.y = element_text(size =20, color = "black")) +
  theme(legend.text = element_text(size =20)) +
  theme(legend.title = element_text(size =20)) +
  labs(color = "pH") +
  ggtitle("Fungi") +
  theme(plot.title = element_text(size =24, face = "bold")) 


ggarrange(vt_16s, vt_its, nrow =1, ncol =2)



#### ====================================================================== ####

#### Taking things out of Phyloseq for plotting ####

 
#### ====================================================================== ####


#creating necessary metadata files from subsets
ps_metadata <- data.frame(sample_data(wo_ps)) 
nh_metadata <- data.frame(sample_data(nh_ps)) 
wopsits_metadata <- data.frame(sample_data(wo_psits))
bacsite19_metadata <-data.frame(sample_data(site19))
bacsite20_metadata <-data.frame(sample_data(site20))
funsite19_metadata <-data.frame(sample_data(funsite19))
funsite20_metadata <-data.frame(sample_data(funsite20))

# Extract the taxonomic information - this will be useful for plotting later
ps_tax_info <- data.frame(tax_table(wo_ps))
nh_tax_info <- data.frame(tax_table(nh_ps))
wopsits_tax_info <- data.frame(tax_table(wo_psits))
bacsite19_tax_info <- data.frame(tax_table(site19))
bacsite20_tax_info <- data.frame(tax_table(site20))
funsite19_tax_info <- data.frame(tax_table(funsite19))
funsite20_tax_info <- data.frame(tax_table(funsite20))


# Extract OTU table from global patterns dataset and convert it out of phyloseq format
ps_otu_df <- data.frame(otu_table(wo_ps), check.names = FALSE) 
nh_otu_df <- data.frame(otu_table(nh_ps), check.names = FALSE) 
wopsits_otu_df <- data.frame(otu_table(wo_psits), check.names = FALSE) 
bacsite19_otu_df <- data.frame(otu_table(site19), check.names = FALSE) 
bacsite20_otu_df <- data.frame(otu_table(site20), check.names = FALSE) 
funsite19_otu_df <- data.frame(otu_table(funsite19), check.names = FALSE) 
funsite20_otu_df <- data.frame(otu_table(funsite20), check.names = FALSE)
 

# Rename the OTUs
# I don't like my OTUs to only be named numbers, so I'll append the word "OTU" in front of them
# this makes it easier so R doesn't mistake them for numbers

# Rename them in the otu table
rownames(ps_otu_df) <- paste0("OTU_", rownames(ps_otu_df))
rownames(nh_otu_df) <- paste0("OTU_", rownames(nh_otu_df))
rownames(wopsits_otu_df) <- paste0("OTU_", rownames(wopsits_otu_df))
rownames(bacsite19_otu_df) <- paste0("OTU_", rownames(bacsite19_otu_df))
rownames(bacsite20_otu_df) <- paste0("OTU_", rownames(bacsite20_otu_df))
rownames(funsite19_otu_df) <- paste0("OTU_", rownames(funsite19_otu_df))
rownames(funsite20_otu_df) <- paste0("OTU_", rownames(funsite20_otu_df))


# Rename them in the taxonomy table
rownames(ps_tax_info) <- paste0("OTU_", rownames(ps_tax_info))
rownames(nh_tax_info) <- paste0("OTU_", rownames(nh_tax_info))
rownames(wopsits_tax_info) <- paste0("OTU_", rownames(wopsits_tax_info))
rownames(bacsite19_tax_info) <- paste0("OTU_", rownames(bacsite19_tax_info))
rownames(bacsite20_tax_info) <- paste0("OTU_", rownames(bacsite20_tax_info))
rownames(funsite19_tax_info) <- paste0("OTU_", rownames(funsite19_tax_info))
rownames(funsite20_tax_info) <- paste0("OTU_", rownames(funsite20_tax_info))


# Transform otu table so species are columns and rows are samples - this is the default that 
# vegan and most other non-microbial ecological packages expects
ps_otu_df_t <- t(ps_otu_df) # note: this converts from dataframe into matrix format

nh_otu_df_t <- t(nh_otu_df)

wopsits_otu_df_t <- t(wopsits_otu_df)

bacsite19_otu_df_t <- t(bacsite19_otu_df)

bacsite20_otu_df_t <- t(bacsite20_otu_df)

funsite19_otu_df_t <- t(funsite19_otu_df)

funsite20_otu_df_t <- t(funsite20_otu_df)





# Bray-curtis dissimilarities
ps_bc <- vegdist(x = ps_otu_df_t, method = "bray")
nh_bc <- vegdist(x = nh_otu_df_t, method = "bray")
wopsits_bc <- vegdist(x = wopsits_otu_df_t, method = "bray")
bacsite19_bc <- vegdist(x = bacsite19_otu_df_t, method = "bray")
bacsite20_bc <- vegdist(x = bacsite20_otu_df_t, method = "bray")
funsite19_bc <- vegdist(x = funsite19_otu_df_t, method = "bray")
funsite20_bc <- vegdist(x = funsite20_otu_df_t, method = "bray")


# PCOA oridination, with bray-curtis
ps_bc.pcoa <- pcoa(ps_bc) # this is from the package ape; 
ps_bc.pct_ex <- round((ps_bc.pcoa$values$Relative_eig) * 100, 1) # this line computes the variation on each axis

nh_bc.pcoa <- pcoa(nh_bc) # this is from the package ape; 
nh_bc.pct_ex <- round((nh_bc.pcoa$values$Relative_eig) * 100, 1) # this line computes the variation on each axis

wopsits_bc.pcoa <- pcoa(wopsits_bc) # this is from the package ape; 
wopsits_bc.pct_ex <- round((wopsits_bc.pcoa$values$Relative_eig) * 100, 1) # this line computes the variation on each axis

bacsite19_bc.pcoa <- pcoa(bacsite19_bc) # this is from the package ape; 
bacsite19_bc.pct_ex <- round((bacsite19_bc.pcoa$values$Relative_eig) * 100, 1) # this line computes the variation on each axis

bacsite20_bc.pcoa <- pcoa(bacsite20_bc) # this is from the package ape; 
bacsite20_bc.pct_ex <- round((bacsite20_bc.pcoa$values$Relative_eig) * 100, 1) # this line computes the variation on each axis


funsite19_bc.pcoa <- pcoa(funsite19_bc) # this is from the package ape; 
funsite19_bc.pct_ex <- round((funsite19_bc.pcoa$values$Relative_eig) * 100, 1) # this line computes the variation on each axis

funsite20_bc.pcoa <- pcoa(funsite20_bc) # this is from the package ape; 
funsite20_bc.pct_ex <- round((funsite20_bc.pcoa$values$Relative_eig) * 100, 1) # this line computes the variation on each axis


#### ====================================================================== ####

#### Step 5: Prepare plotting dataframes #####
# here we join our sample metadata with the points that were generated by the ordination calculations above
# For simplicity we'll actually just create 1 dataframe with all of the different sets of ordination points
# but you could also make them separately
#### ====================================================================== ####

ps_ord_plot.df <- ps_metadata %>%
  # First join the pcoa data
  left_join(ps_bc.pcoa$vectors[,1:2] %>% data.frame() %>%
              rownames_to_column(var = "SampleID"), # here we make a SampleID column before joining the dataframes
            by = "SampleID") %>%
  dplyr::rename(PCOA1_BC = Axis.1, PCOA2_BC = Axis.2)

nh_ord_plot.df <- nh_metadata %>%
  # First join the pcoa data
  left_join(nh_bc.pcoa$vectors[,1:2] %>% data.frame() %>%
              rownames_to_column(var = "SampleID"), # here we make a SampleID column before joining the dataframes
            by = "SampleID") %>%
  dplyr::rename(PCOA1_BC = Axis.1, PCOA2_BC = Axis.2)

wopsits_ord_plot.df <- wopsits_metadata %>%
  # First join the pcoa data
  left_join(wopsits_bc.pcoa$vectors[,1:2] %>% data.frame() %>%
              rownames_to_column(var = "SampleID"), # here we make a SampleID column before joining the dataframes
            by = "SampleID") %>%
  dplyr::rename(PCOA1_BC = Axis.1, PCOA2_BC = Axis.2)


bacsite19_ord_plot.df <- bacsite19_metadata %>%
  # First join the pcoa data
  left_join(bacsite19_bc.pcoa$vectors[,1:2] %>% data.frame() %>%
              rownames_to_column(var = "SampleID"), # here we make a SampleID column before joining the dataframes
            by = "SampleID") %>%
  dplyr::rename(PCOA1_BC = Axis.1, PCOA2_BC = Axis.2)

bacsite20_ord_plot.df <- bacsite20_metadata %>%
  # First join the pcoa data
  left_join(bacsite20_bc.pcoa$vectors[,1:2] %>% data.frame() %>%
              rownames_to_column(var = "SampleID"), # here we make a SampleID column before joining the dataframes
            by = "SampleID") %>%
  dplyr::rename(PCOA1_BC = Axis.1, PCOA2_BC = Axis.2)


funsite19_ord_plot.df <- funsite19_metadata %>%
  # First join the pcoa data
  left_join(funsite19_bc.pcoa$vectors[,1:2] %>% data.frame() %>%
              rownames_to_column(var = "SampleID"), # here we make a SampleID column before joining the dataframes
            by = "SampleID") %>%
  dplyr::rename(PCOA1_BC = Axis.1, PCOA2_BC = Axis.2)

funsite20_ord_plot.df <- funsite20_metadata %>%
  # First join the pcoa data
  left_join(funsite20_bc.pcoa$vectors[,1:2] %>% data.frame() %>%
              rownames_to_column(var = "SampleID"), # here we make a SampleID column before joining the dataframes
            by = "SampleID") %>%
  dplyr::rename(PCOA1_BC = Axis.1, PCOA2_BC = Axis.2)



#### ====================================================================== ####

#### Compute Arrows ####
#### ====================================================================== ####
# Learn more about the ins and outs of this here:
# https://www.davidzeleny.net/anadat-r/doku.php/en:suppl_vars_examples

# For PCOA arrows are essentially the correlation between the axis scores and the 
# environmental variable or OTU abundance. For NMDS, it's not quite as simple since
# NMDS's are caluated a little differently

# Note: Arrows assume that there is a linear trend in the variable you are fitting
# across the surface of the plot, this isn't always (or even often) true in unconstrained
# ordination techniques such as PCOA and NMDS, sometimes a better visual representation is 
# to use a more complex visualization of the gradient with ordisurf. 
# See: https://fromthebottomoftheheap.net/2011/06/10/what-is-ordisurf-doing/
# for more information; and https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html
# for an example of how to convert ordisurf output into ggplot

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
  pvals.adj <- p.adjust(fit$vectors$pvals, method = "fdr")
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




ps_bc.pcoa.arrows <- env_fit_wrapper(ordination = ps_bc.pcoa,
                                       ordination_type = "PCOA",
                                       otu_tab = ps_otu_df_t, 
                                       tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but nice
  left_join(ps_tax_info %>% 
              rownames_to_column(var = "OTU_ID"),
            by = c("FactorLabel" = "OTU_ID"))


nh_bc.pcoa.arrows <- env_fit_wrapper(ordination = nh_bc.pcoa,
                                     ordination_type = "PCOA",
                                     otu_tab = nh_otu_df_t, 
                                     tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but nice
  left_join(nh_tax_info %>% 
              rownames_to_column(var = "OTU_ID"),
            by = c("FactorLabel" = "OTU_ID"))


wopsits_bc.pcoa.arrows <- env_fit_wrapper(ordination = wopsits_bc.pcoa,
                                     ordination_type = "PCOA",
                                     otu_tab = wopsits_otu_df_t, 
                                     tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but nice
  left_join(wopsits_tax_info %>% 
              rownames_to_column(var = "OTU_ID"),
            by = c("FactorLabel" = "OTU_ID"))


bacsite19_bc.pcoa.arrows <- env_fit_wrapper(ordination = bacsite19_bc.pcoa,
                                          ordination_type = "PCOA",
                                          otu_tab = bacsite19_otu_df_t, 
                                          tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but nice
  left_join(bacsite19_tax_info %>% 
              rownames_to_column(var = "OTU_ID"),
            by = c("FactorLabel" = "OTU_ID"))


bacsite20_bc.pcoa.arrows <- env_fit_wrapper(ordination = bacsite20_bc.pcoa,
                                            ordination_type = "PCOA",
                                            otu_tab = bacsite20_otu_df_t, 
                                            tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but nice
  left_join(bacsite20_tax_info %>% 
              rownames_to_column(var = "OTU_ID"),
            by = c("FactorLabel" = "OTU_ID"))



funsite19_bc.pcoa.arrows <- env_fit_wrapper(ordination = funsite19_bc.pcoa,
                                            ordination_type = "PCOA",
                                            otu_tab = funsite19_otu_df_t, 
                                            tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but nice
  left_join(funsite19_tax_info %>% 
              rownames_to_column(var = "OTU_ID"),
            by = c("FactorLabel" = "OTU_ID"))


funsite20_bc.pcoa.arrows <- env_fit_wrapper(ordination = funsite20_bc.pcoa,
                                            ordination_type = "PCOA",
                                            otu_tab = funsite20_otu_df_t, 
                                            tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but nice
  left_join(funsite20_tax_info %>% 
              rownames_to_column(var = "OTU_ID"),
            by = c("FactorLabel" = "OTU_ID"))


#### ====================================================================== ####

#### Step 5: Plot data #####
#### ====================================================================== ####


# Woodman samples
ggplot(data = ps_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = year)) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", ps_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", ps_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = ps_bc.pcoa.arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text(data = ps_bc.pcoa.arrows,
            mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel)) +
  # Alternatively
  geom_text_repel(data = ps_bc.pcoa.arrows,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel))




#Samples by site type (nh samples)
ggplot(data = nh_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = site_type, shape = year.plot)) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", nh_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", nh_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = nh_bc.pcoa.arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text(data = nh_bc.pcoa.arrows,
            mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel)) +
  # Alternatively
  geom_text_repel(data = nh_bc.pcoa.arrows,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel))


## ITS Woodman samples
ggplot(data = wopsits_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(shape = year, color = host)) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", wopsits_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", wopsits_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = wopsits_bc.pcoa.arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text(data = wopsits_bc.pcoa.arrows,
            mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel)) +
  # Alternatively
  geom_text_repel(data = wopsits_bc.pcoa.arrows,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel))




#### ====================================================================== ####

# Plot by taxonomic group
#### ====================================================================== ####


# Summarize OTU table by class
ClassTab <- t(ps_otu_df_t) %>% 
  data.frame() %>%
  rownames_to_column(var = "OTU_ID") %>%
  left_join(ps_tax_info %>% rownames_to_column(var = "OTU_ID"),
            by = "OTU_ID") %>%
  mutate(Class = ifelse(is.na(Class), "ClasssNA", Class)) %>% # fix NA classes
  group_by(Class) %>%
  summarize(across(X1102B_16s:X6144Rz, ~sum(.))) %>% # sum class abundances in each sample column
  column_to_rownames(var = "Class") %>%
  t() # retransform data


# Summarize OTU table by class
nhClassTab <- t(nh_otu_df_t) %>% 
  data.frame() %>%
  rownames_to_column(var = "OTU_ID") %>%
  left_join(nh_tax_info %>% rownames_to_column(var = "OTU_ID"),
            by = "OTU_ID") %>%
  mutate(Class = ifelse(is.na(Class), "ClasssNA", Class)) %>% # fix NA classes
  group_by(Class) %>%
  summarize(across(X1102B_16s:RyeA3Rz2, ~sum(.))) %>% # sum class abundances in each sample column
  column_to_rownames(var = "Class") %>%
  t() # retransform data


# Summarize OTU table by class
wopsitsClassTab <- t(wopsits_otu_df_t) %>% 
  data.frame() %>%
  rownames_to_column(var = "OTU_ID") %>%
  left_join(wopsits_tax_info %>% rownames_to_column(var = "OTU_ID"),
            by = "OTU_ID") %>%
  mutate(Class = ifelse(is.na(Class), "ClasssNA", Class)) %>% # fix NA classes
  group_by(Class) %>%
  summarize(across(X1102B_ITS:X691Rz, ~sum(.))) %>% # sum class abundances in each sample column
  column_to_rownames(var = "Class") %>%
  t() # retransform data


# Summarize OTU table by class - Bacteria by site type 
bacsite19ClassTab <- t(bacsite19_otu_df_t) %>% 
  data.frame() %>%
  rownames_to_column(var = "OTU_ID") %>%
  left_join(bacsite19_tax_info %>% rownames_to_column(var = "OTU_ID"),
            by = "OTU_ID") %>%
  mutate(Class = ifelse(is.na(Class), "ClasssNA", Class)) %>% # fix NA classes
  group_by(Class) %>%
  summarize(across(X1101B:RyeA3Rz2, ~sum(.))) %>% # sum class abundances in each sample column
  column_to_rownames(var = "Class") %>%
  t() # retransform data


bacsite20ClassTab <- t(bacsite20_otu_df_t) %>% 
  data.frame() %>%
  rownames_to_column(var = "OTU_ID") %>%
  left_join(bacsite20_tax_info %>% rownames_to_column(var = "OTU_ID"),
            by = "OTU_ID") %>%
  mutate(Class = ifelse(is.na(Class), "ClasssNA", Class)) %>% # fix NA classes
  group_by(Class) %>%
  summarize(across(X1102B_16s:RyeB5B4_16s, ~sum(.))) %>% # sum class abundances in each sample column
  column_to_rownames(var = "Class") %>%
  t() # retransform data



# Summarize OTU table by class - Fungi by site type 
funsite19ClassTab <- t(funsite19_otu_df_t) %>% 
  data.frame() %>%
  rownames_to_column(var = "OTU_ID") %>%
  left_join(funsite19_tax_info %>% rownames_to_column(var = "OTU_ID"),
            by = "OTU_ID") %>%
  mutate(Class = ifelse(is.na(Class), "ClasssNA", Class)) %>% # fix NA classes
  group_by(Class) %>%
  summarize(across(X1101Rz:RyeA3Rz, ~sum(.))) %>% # sum class abundances in each sample column
  column_to_rownames(var = "Class") %>%
  t() # retransform data


funsite20ClassTab <- t(funsite20_otu_df_t) %>% 
  data.frame() %>%
  rownames_to_column(var = "OTU_ID") %>%
  left_join(funsite20_tax_info %>% rownames_to_column(var = "OTU_ID"),
            by = "OTU_ID") %>%
  mutate(Class = ifelse(is.na(Class), "ClasssNA", Class)) %>% # fix NA classes
  group_by(Class) %>%
  summarize(across(X1102B_ITS:RyeB8B_ITS, ~sum(.))) %>% # sum class abundances in each sample column
  column_to_rownames(var = "Class") %>%
  t() # retransform data







# NH samples
nh_bc.pcoa.class.arrows <- env_fit_wrapper(ordination = nh_bc.pcoa,
                                           ordination_type = "PCOA",
                                           otu_tab = nhClassTab, 
                                           tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but  nice
  left_join(nh_tax_info %>% 
              rownames_to_column(var = "OTU_ID") %>%
              select(Phylum:Class) %>% distinct(),
            by = c("FactorLabel" = "Class"))

#Remove the underscores from the beginning of the Factor Labels
nh_bc.pcoa.class.arrows$FactorLabel<-str_sub(nh_bc.pcoa.class.arrows$FactorLabel,6)


## ITS samples - Woodman
# NH samples
wopsits_bc.pcoa.class.arrows <- env_fit_wrapper(ordination = wopsits_bc.pcoa,
                                           ordination_type = "PCOA",
                                           otu_tab = wopsitsClassTab, 
                                           tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but  nice
  left_join(wopsits_tax_info %>% 
              rownames_to_column(var = "OTU_ID") %>%
              select(Phylum:Class) %>% distinct(),
            by = c("FactorLabel" = "Class"))

#Remove the underscores from the beginning of the Factor Labels
wopsits_bc.pcoa.class.arrows$FactorLabel<-str_sub(wopsits_bc.pcoa.class.arrows$FactorLabel,4)
wopsits_bc.pcoa.class.arrows<-subset(wopsits_bc.pcoa.class.arrows, FactorLabel != "ssNA")



# 16s samples by site type - 2019 and 2020
bacsite19_bc.pcoa.class.arrows <- env_fit_wrapper(ordination = bacsite19_bc.pcoa,
                                                ordination_type = "PCOA",
                                                otu_tab = bacsite19ClassTab, 
                                                tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but  nice
  left_join(bacsite19_tax_info %>% 
              rownames_to_column(var = "OTU_ID") %>%
              select(Phylum:Class) %>% distinct(),
            by = c("FactorLabel" = "Class"))

#Remove the underscores from the beginning of the Factor Labels
bacsite19_bc.pcoa.class.arrows$FactorLabel<-str_sub(bacsite19_bc.pcoa.class.arrows$FactorLabel,6)


bacsite20_bc.pcoa.class.arrows <- env_fit_wrapper(ordination = bacsite20_bc.pcoa,
                                                  ordination_type = "PCOA",
                                                  otu_tab = bacsite20ClassTab, 
                                                  tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but  nice
  left_join(bacsite20_tax_info %>% 
              rownames_to_column(var = "OTU_ID") %>%
              select(Phylum:Class) %>% distinct(),
            by = c("FactorLabel" = "Class"))

#Remove the underscores from the beginning of the Factor Labels
bacsite20_bc.pcoa.class.arrows$FactorLabel<-str_sub(bacsite20_bc.pcoa.class.arrows$FactorLabel,6)



# ITS samples by site type - 2019 and 2020
funsite19_bc.pcoa.class.arrows <- env_fit_wrapper(ordination = funsite19_bc.pcoa,
                                                  ordination_type = "PCOA",
                                                  otu_tab = funsite19ClassTab, 
                                                  tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but  nice
  left_join(funsite19_tax_info %>% 
              rownames_to_column(var = "OTU_ID") %>%
              select(Phylum:Class) %>% distinct(),
            by = c("FactorLabel" = "Class"))

#Remove the underscores from the beginning of the Factor Labels
funsite19_bc.pcoa.class.arrows$FactorLabel<-str_sub(funsite19_bc.pcoa.class.arrows$FactorLabel,4)



funsite20_bc.pcoa.class.arrows <- env_fit_wrapper(ordination = funsite20_bc.pcoa,
                                                  ordination_type = "PCOA",
                                                  otu_tab = funsite20ClassTab, 
                                                  tax_cutoff = 10) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but  nice
  left_join(funsite20_tax_info %>% 
              rownames_to_column(var = "OTU_ID") %>%
              select(Phylum:Class) %>% distinct(),
            by = c("FactorLabel" = "Class"))

#Remove the underscores from the beginning of the Factor Labels
funsite20_bc.pcoa.class.arrows$FactorLabel<-str_sub(funsite20_bc.pcoa.class.arrows$FactorLabel,4)
funsite20_bc.pcoa.class.arrows<-subset(funsite20_bc.pcoa.class.arrows, FactorLabel != "ssNA")







#### ====================================================================== ####

#### Step 5: Plot data #####
#### ====================================================================== ####



# Plot Arrows by Class
# PCOA Bray-Curtis
p_bac<-ggplot(data = ps_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = year), size =2) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", ps_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", ps_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = ps_bc.pcoa.class.arrows,
               x = 0, y = 0.06, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text_repel(data = ps_bc.pcoa.class.arrows, max.overlaps = 20,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel, size =45)) +
  theme(axis.title.x = element_text(color="black", size = 22)) +
  theme(axis.title.y = element_text(color="black", size =22)) +
  theme(axis.text.x = element_text(size =20, color = "black")) +
  theme(axis.text.y = element_text(size =20, color = "black")) +
  theme(legend.text = element_text(size =20)) +
  theme(legend.title = element_blank()) +
  scale_fill_discrete(labels = c("2019 Plot", "2020 Plot")) +
  guides(size = "none") +
  ggtitle("Bacteria") +
  theme(plot.title = element_text(size =24, face = "bold"))



## NH Samples
# Plot Arrows by Class
# PCOA Bray-Curtis

# All NH samples 2019 + 2020
nh_bac<-ggplot(data = nh_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = site_type, shape = year.plot), size =2) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", nh_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", nh_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = nh_bc.pcoa.class.arrows,
               x = 0, y = 0.06, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text_repel(data = nh_bc.pcoa.class.arrows, max.overlaps = 20,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel, size =35)) +
  theme(axis.title.x = element_text(color="black", size = 20)) +
  theme(axis.title.y = element_text(color="black", size =20)) +
  theme(axis.text.x = element_text(size =12, color = "black")) +
  theme(axis.text.y = element_text(size =12, color = "black")) +
  theme(legend.text = element_text(size =14)) +
  theme(legend.title = element_blank())# +
  #scale_fill_discrete(labels = c("Plot 1", "Plot 2"))


# 2019 NH Samples - 16s
bacsite19<-ggplot(data = bacsite19_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = site_type, shape = location), size =2) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", bacsite19_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", bacsite19_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = bacsite19_bc.pcoa.class.arrows,
               x = 0, y = 0.06, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text_repel(data = bacsite19_bc.pcoa.class.arrows, max.overlaps = 20,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel, size =35)) +
  theme(axis.title.x = element_text(color="black", size = 20)) +
  theme(axis.title.y = element_text(color="black", size =20)) +
  theme(axis.text.x = element_text(size =12, color = "black")) +
  theme(axis.text.y = element_text(size =12, color = "black")) +
  theme(legend.text = element_text(size =14)) +
  theme(legend.title = element_text(size = 14)) +
  guides(size = "none") +
  ggtitle("Bacteria 2019/Plot1") +
  theme(plot.title = element_text(size =20, face = "bold"))



# 2020 NH Samples -16s
bacsite20<-ggplot(data = bacsite20_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = site_type, shape = location), size =2) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", bacsite20_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", bacsite20_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = bacsite20_bc.pcoa.class.arrows,
               x = 0, y = 0.06, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text_repel(data = bacsite20_bc.pcoa.class.arrows, max.overlaps = 20,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel, size =35)) +
  theme(axis.title.x = element_text(color="black", size = 20)) +
  theme(axis.title.y = element_text(color="black", size =20)) +
  theme(axis.text.x = element_text(size =12, color = "black")) +
  theme(axis.text.y = element_text(size =12, color = "black")) +
  theme(legend.text = element_text(size =14)) +
  theme(legend.title = element_text(size = 14)) +
  guides(size = "none") +
  ggtitle("Bacteria 2020/Plot2") +
  theme(plot.title = element_text(size =20, face = "bold"))



# 2019 NH Samples - ITS
p_funsite19<-ggplot(data = funsite19_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = site_type, shape = location), size =2) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", funsite19_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", funsite19_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = funsite19_bc.pcoa.class.arrows,
               x = 0, y = 0.06, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text_repel(data = funsite19_bc.pcoa.class.arrows, max.overlaps = 20,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel, size =35)) +
  theme(axis.title.x = element_text(color="black", size = 20)) +
  theme(axis.title.y = element_text(color="black", size =20)) +
  theme(axis.text.x = element_text(size =12, color = "black")) +
  theme(axis.text.y = element_text(size =12, color = "black")) +
  theme(legend.text = element_text(size =14)) +
  theme(legend.title = element_text(size = 14)) +
  guides(size = "none") +
  ggtitle("Fungi 2019/Plot1") +
  theme(plot.title = element_text(size =20, face = "bold"))


# 2020 sample by site - ITS
p_funsite20<-ggplot(data = funsite20_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = site_type, shape = location), size =2) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", funsite20_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", funsite20_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = funsite20_bc.pcoa.class.arrows,
               x = 0, y = 0.06, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text_repel(data = funsite20_bc.pcoa.class.arrows, max.overlaps = 20,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel, size =35)) +
  theme(axis.title.x = element_text(color="black", size = 20)) +
  theme(axis.title.y = element_text(color="black", size =20)) +
  theme(axis.text.x = element_text(size =12, color = "black")) +
  theme(axis.text.y = element_text(size =12, color = "black")) +
  theme(legend.text = element_text(size =14)) +
  theme(legend.title = element_text(size = 14)) +
  guides(size = "none") +
  ggtitle("Fungi 2020/Plot2") +
  theme(plot.title = element_text(size =20, face = "bold"))


ggarrange(bacsite19, bacsite20, p_funsite19, p_funsite20,
          ncol = 2, nrow =2, common.legend = TRUE, legend = "right")





## ITS - Woodman
# PCOA Bray-Curtis
p_fun<-ggplot(data = wopsits_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = year.plot), size =2) +
  scale_fill_discrete(labels = c("Plot 1", "Plot 2")) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", wopsits_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", wopsits_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = wopsits_bc.pcoa.class.arrows,
               x = 0, y = 0.06, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text_repel(data = wopsits_bc.pcoa.class.arrows, max.overlaps = 20,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel, size =45)) +
  theme(axis.title.x = element_text(color="black", size = 22)) +
  theme(axis.title.y = element_text(color="black", size =22)) +
  theme(axis.text.x = element_text(size =20, color = "black")) +
  theme(axis.text.y = element_text(size =20, color = "black")) +
  theme(legend.text = element_text(size =20)) +
  theme(legend.title = element_blank()) +
  guides(size = "none") +
  ggtitle("Fungi") +
  theme(plot.title = element_text(size =24, face = "bold"))


ggarrange(p_bac, p_fun, nrow =1 , ncol =2, common.legend = TRUE, legend = "right")

