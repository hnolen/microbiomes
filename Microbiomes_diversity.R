#### Diversity

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
library(car)


#### ====================================================================== ####

#### Load 16S data ####
#### ====================================================================== ####

#reading in feature table, taxonomy table, and metadata
feat_table = read.csv("mergedredo/12500/split_feat_table.csv", sep = ",", row.names =1)
feat_table = as.matrix(feat_table)

taxonomy = read.csv("mergedredo/12500/split_tax.csv", sep = ",", row.names = 1)
taxonomy <- replace(taxonomy, taxonomy == "", NA)
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

physeq_its = phyloseq(OTUits, TAXits, METAits)


##looking at just woodman samples
wo_psits<- subset_samples(physeq_its, location == "Durham_NH")
wo_psits<- subset_samples(wo_psits, site_type == "agricultural")



###16s diversity analyses - looking at diffs in woodman diversity in 2019 and 2020
wo_ps<-subset_samples(ps, site == "Woodman")

#get richness estimates
richness_est<-estimate_richness(wo_ps, split = TRUE)


richness_est_c <- tibble::rownames_to_column(richness_est, "index")
metadata_c <- tibble::rownames_to_column(metadata_full, "index")

#join richness data and metadata
rich_w_meta <- left_join(richness_est_c, metadata_c, by = "index") #now I can do stats on the diversity dataframe!

#first things first: make the first column the rownames
row.names(rich_w_meta) <- rich_w_meta$index
rich_w_meta[1] <- NULL


##want to do anova on sample type (soil type) on each year separately to get p-values
rich_w_meta$year<-as.factor(rich_w_meta$year)
div2019<-subset(rich_w_meta, year == "2019")
mod19<- lm(Shannon~soil_type, data = div2019)
Anova(mod19) # p = 0.21

hmod19<- lm(Shannon~host_type, data = div2019)
Anova(hmod19) #host type p = 0.61

hsmod19<- lm(Shannon~host, data = div2019)
Anova(hsmod19) #host p = 0.8

div2020<-subset(rich_w_meta, year == "2020")
mod20<- lm(Shannon~soil_type, data = div2020)
Anova(mod20) #p = 0.37

hmod20<- lm(Shannon~host_type, data = div2020)
Anova(hmod20) #host p = 0.63

hsmod20<- lm(Shannon~host, data = div2020)
Anova(hsmod20) #host p = 0.73

mod <- lm(Shannon~soil_type*year, data = rich_w_meta)
library(car)
Anova(mod) #year p = 0.001, soil_type p = 0.85


##plot diversity
divp<-ggplot(rich_w_meta, aes(x = year, y = Shannon, fill = year)) + 
  geom_boxplot() +
  theme_bw() +
  ylim(2,8) +
  labs(x = "", y = "Shannon Richness") +
  theme(axis.title.x = element_text(color = "black", size = 22), 
        axis.title.y = element_text(color = "black", size = 22), 
        axis.text = element_text(size = 20, colour = "black")) +
  stat_signif(annotation = "**", textsize = 8, xmin = 1, xmax = 2, y_position = 7.8) +
  scale_x_discrete(labels = c("2019 Plot", "2020 Plot")) +
  theme(legend.position = "none") +
  ggtitle("Bacteria") +
  theme(plot.title = element_text(size =24, face = "bold"))




#### ====================================================================== ####

###ITS diversity analyses - looking at diffs in woodman diversity in 2019 and 2020

#### ====================================================================== ####

richness_its<-estimate_richness(wo_psits, split = TRUE)
write.csv(richness_its, "./QDM011/merged/its/fungi_div_est.csv", row.names = TRUE) #need to get X's out of names
richness_its = read.csv("./QDM011/merged/its/fungi_div_est.csv", sep = ",", row.names =1, check.names = FALSE)



richness_its_c <- tibble::rownames_to_column(richness_its, "index")
metadata_itsc <- tibble::rownames_to_column(metadata_fullits, "index")

rich_w_metaits <- left_join(richness_its_c, metadata_itsc, by = "index") 

#first things first: make the first column the rownames
row.names(rich_w_metaits) <- rich_w_metaits$index
rich_w_metaits[1] <- NULL


##want to do anova on sample type (soil type) on each year separately to get p-values
rich_w_metaits$year.plot<-as.factor(rich_w_metaits$year.plot)
divits2019<-subset(rich_w_metaits, year.plot == "2019/Plot1")
modits19<- lm(Shannon~soil_type, data = divits2019)
Anova(modits19) # p = 0.18

hitsmod19<- lm(Shannon~host_type, data = divits2019)
Anova(hitsmod19) #host type p = 0.61

hsitsmod19<- lm(Shannon~host, data = divits2019)
Anova(hsitsmod19) #host p = 0.81

divits2020<-subset(rich_w_metaits, year.plot == "2020/Plot2")
modits20<- lm(Shannon~soil_type, data = divits2020)
Anova(modits20) #p = 0.32

hitsmod20<- lm(Shannon~host_type, data = divits2020)
Anova(hitsmod20) #host type p = 0.53

hsitsmod20<- lm(Shannon~host, data = divits2020)
Anova(hsitsmod20) #host p = 0.19

modits <- lm(Shannon~soil_type*year.plot, data = rich_w_metaits)
library(car)
Anova(modits) #year p = 0.52, soil_type p = 0.75


##plotting
divpits<-ggplot(rich_w_metaits, aes(x = year.plot, y = Shannon, fill = year.plot)) + 
  geom_boxplot() +
  theme_bw() +
  ylim(2,8) +
  labs(x = "", y = "") +
  theme(axis.title.x = element_text(color = "black", size = 22), 
        axis.title.y = element_text(color = "black", size = 22), 
        axis.text = element_text(size = 20, colour = "black")) +
  scale_x_discrete(labels = c("2019 Plot", "2020 Plot")) +
  theme(legend.position = "none") +
  ggtitle("Fungi") +
  theme(plot.title = element_text(size =24, face = "bold"))




##putting 16s and ITS plots together
ggarrange(divp, divpits, nrow = 1, ncol =2)




#### ====================================================================== ####

#### Location - ITS

#### ====================================================================== ####



#### Location
##diversity - everything but woodman - location experiment

#remove woodman samples from dataset
loc_ps <- subset_samples(ps, site != "Woodman")

#estimate richness
loc_rich<-estimate_richness(loc_ps, split=TRUE)


metadata <- as(sample_data(loc_ps), "data.frame")

loc_div_c <- tibble::rownames_to_column(loc_rich, "index")
metadata_c <- tibble::rownames_to_column(metadata, "index")

loc_div<- left_join(loc_div_c, metadata_c, by = "index")

#first things first: make the first column the rownames
row.names(loc_div) <- loc_div$index
loc_div[1] <- NULL



mod <- lm(Shannon~location+texture, data = loc_div)
Anova(mod) #location p = 0.33, texture p = 0.002
tuk<-LSD.test(mod, "location")

#plotting
divp<-ggplot(loc_div, aes(x = factor(location), y = Shannon)) + 
  geom_boxplot() +
  ylim(0,8) +
  labs(x = "", y = "Shannon Richness") +
  theme(axis.title.x = element_text(color = "black", size = 15), 
        axis.title.y = element_text(color = "black", size = 15), 
        axis.text = element_text(size = 13, colour = "black")) 

tex_nona<- subset(loc_div, texture != "na")
tex_nona<- subset(tex_nona, texture != "clay loam")
tex_nona<- subset(tex_nona, texture != "silty clay")
tex_nona<- subset(tex_nona, texture != "silt loam")
mod <- lm(Shannon~texture, data = tex_nona)
Anova(mod) # location p = 0.0014, 

tuktex<-HSD.test(mod, "texture")
t_tex<-tuktex$groups
t_tex<-tibble::rownames_to_column(t_tex, "trtmnt")
str(t_tex)

divtex<-ggplot(tex_nona, aes(x = factor(texture), y = Shannon)) + 
  geom_boxplot() +
  theme_bw() +
  ylim(0,8) +
  labs(x = "", y = "Shannon Richness") +
  theme(axis.title.y = element_text(color = "black", size = 20), 
        axis.text = element_text(size = 20, colour = "black")) +
  annotate("text", x = 1, y = 6.8, label = "a", size = 8) +
  annotate("text", x = 2, y = 7.3, label = "ab", size = 8) +
  annotate("text", x = 3, y = 6.7, label = "b", size = 8) +
  annotate("text", x = 4, y = 7.2, label = "a", size = 8)




#### ====================================================================== ####
##### Diversity by Site type - agriculture vs wild - 16s

#### ====================================================================== ####
du_ps<-subset_samples(ps, location == "Durham_NH")
rye_ps<-subset_samples(ps, location == "Rye_NH")
nh_ps = merge_phyloseq(du_ps, rye_ps)

site_rich<-estimate_richness(nh_ps, split=TRUE)
write.csv(site_rich, "./mergedredo/12500/site_div.csv", row.names = TRUE)
site_div = read.csv("./mergedredo/12500/site_div.csv", sep = ",", row.names =1, check.names = FALSE)

metadata <- as(sample_data(nh_ps), "data.frame")

site_div_c <- tibble::rownames_to_column(site_div, "index")
metadata_c <- tibble::rownames_to_column(metadata, "index")

site_div<- left_join(site_div_c, metadata_c, by = "index")

#first things first: make the first column the rownames
row.names(site_div) <- site_div$index
site_div[1] <- NULL

#make soil_type an ordered variable (not sure why we need to do this)
site_div$soil_type <- ordered(site_div$soil_type, c("bulk", "rhizosphere"))

mod <- lm(Shannon~site_type*year, data = site_div)
Anova(mod) #site type p < 0.001 year p = 0.019


## subsetting into 2019 and 2020
site19_div<- subset(site_div, year.plot == "2019/Plot1")


site20_div<- subset(site_div, year.plot == "2020/Plot2")

mod19 <- lm(Shannon~soil_type*site_type, data = site19_div)
Anova(mod19) #soil type p = 0.4, site type p < 0.001

mod20 <- lm(Shannon~soil_type*site_type, data = site20_div)
Anova(mod20) #soil type p = 0.8, site type p < 0.001, soil_type:site_type p = 0.05

bsite20_div <- subset(site20_div, soil_type == "bulk")
rzsite20_div <- subset(site20_div, soil_type == "rhizosphere")

bmod20 <- lm(Shannon~site_type, data = site20_div)
Anova(mod20)

div19p<-ggplot(site19_div, aes(x = factor(soil_type), y = Shannon, fill = site_type)) + 
  geom_boxplot() +
  theme_bw() +
  ylim(0,8) +
  labs(x = "", y = "Shannon Richness") +
  theme(axis.title.x = element_text(color = "black", size = 22), 
        axis.title.y = element_text(color = "black", size = 22), 
        axis.text = element_text(size = 20, colour = "black")) +
  ggtitle("2019") +
  theme(plot.title = element_text(size =24, face = "bold")) +
  theme(legend.text = element_text(size =20)) +
  theme(legend.title = element_text(size =20)) +
  labs(fill = "Site type") +
  stat_signif(annotation = "***", textsize = 8, xmin = 1.81, xmax = 2.19, y_position = 7.5) +
  stat_signif(annotation = "***", textsize = 8, xmin = 0.81, xmax = 1.19, y_position = 7.5)


div20p<-ggplot(site20_div, aes(x = factor(soil_type), y = Shannon, fill = site_type)) + 
  geom_boxplot() +
  theme_bw() +
  ylim(0,8) +
  labs(x = "", y = "Shannon Richness") +
  theme(axis.title.x = element_text(color = "black", size = 22), 
        axis.title.y = element_text(color = "black", size = 22), 
        axis.text = element_text(size = 20, colour = "black")) +
  ggtitle("2020") +
  theme(plot.title = element_text(size =24, face = "bold")) +
  theme(legend.text = element_text(size =20)) +
  theme(legend.title = element_text(size =20)) +
  labs(fill = "Site type") +
  stat_signif(annotation = "***", textsize = 8, xmin = 1.81, xmax = 2.19, y_position = 7.8) +
  stat_signif(annotation = "***", textsize = 8, xmin = 0.81, xmax = 1.19, y_position = 7.8)

ggarrange(div19p, div20p, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")


#### ====================================================================== ####

##### Site type - ITS

#### ====================================================================== ####


du_ps<-subset_samples(physeq_its, location == "Durham_NH")
rye_ps<-subset_samples(physeq_its, location == "Rye_NH")
nh_ps = merge_phyloseq(du_ps, rye_ps)

site_rich<-estimate_richness(nh_ps, split=TRUE)
write.csv(site_rich, "./QDM011/merged/its/site_div.csv", row.names = TRUE)
site_div = read.csv("./QDM011/merged/its/site_div.csv", sep = ",", row.names =1, check.names = FALSE)


metadata <- as(sample_data(nh_ps), "data.frame")

site_div_c <- tibble::rownames_to_column(site_div, "index")
metadata_c <- tibble::rownames_to_column(metadata, "index")

site_div<- left_join(site_div_c, metadata_c, by = "index")

#first things first: make the first column the rownames
row.names(site_div) <- site_div$index
site_div[1] <- NULL

#make soil_type an ordered variable (not sure why we need to do this)
site_div$soil_type <- ordered(site_div$soil_type, c("bulk", "rhizosphere"))

mod <- lm(Shannon~site_type*year.plot, data = site_div)
Anova(mod) #site type by year p = 0.03


## subsetting into 2019 and 2020
site19_div<- subset(site_div, year.plot == "2019/Plot1")


site20_div<- subset(site_div, year.plot == "2020/Plot2")

mod19 <- lm(Shannon~soil_type*site_type, data = site19_div)
Anova(mod19) #soil type p = 0.2, site type p =0.15

mod20 <- lm(Shannon~soil_type*site_type, data = site20_div)
Anova(mod20) #soil type p = 0.2, site typep = 0.1

bsite20_div <- subset(site20_div, soil_type == "bulk")
rzsite20_div <- subset(site20_div, soil_type == "rhizosphere")

bmod20 <- lm(Shannon~site_type, data = site20_div)
Anova(mod20)







