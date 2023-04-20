### Microbiome differential abundance work

##Differential abundance refer to this  - http://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html

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



##subset to get Woodman data

wo_ps<- subset_samples(ps, location == "Durham_NH")
wo_ps<- subset_samples(wo_ps, site_type == "agricultural")

#separate into 2019 and 2020 data
wops19<-subset_samples(wo_ps, year == "2019")
wops20<-subset_samples(wo_ps, year == "2020")

###significant effect of soil type (bulk vs rhiz) in 2019 woodman samples - differential abundance
##to create waterfall plots comparing diff. abundance across two groups

##get the taxonomy table from the phyloseq object of interest - I'm looking at genus level taxonomy - make sure any empty string reads as "Unmatched genus"
tax_table(wops19)[tax_table(wops19)[,"Genus"]== "","Genus"] <- "Unmatchedgenus"

#create a new phyloseq object from that genus data
genus_dat = aggregate_taxa(wops19, "Genus")

#perform the differential abundance analyis, formula and group should be the factor you are interested in
# contains: 1) log fold changes; 2) standard errors; 3) test statistics; 
# 4) p-values; 5) adjusted p-values; 6) indicators whether the taxon is differentially abundant (TRUE) or not (FALSE).
out<-ANCOMBC::ancombc(phyloseq = genus_dat, formula = "soil_type", 
                      p_adj_method = "holm", prv_cut = 0.1, #holm is default p-value adjustment, taxa with prevalence less that prv_cut will be excluded, default 0.1
                      lib_cut = 1000, group = "soil_type", #samples with library sizes less than lib cut will be excluded
                      struc_zero = TRUE, neg_lb = TRUE,
                      tol = 1e-5, max_iter = 100, #tol - convergence tolerance - default 1e-5, max_iter-maximum iterations, default 100
                      conserve = TRUE, alpha = 0.05, global = TRUE) #conserve = TRUE uses a conservative variance estimator for test statistic - recommended for small sample sizes
res<-out$res
kable(head(res$diff_abn))
lfc<-res$lfc
head(lfc)

#create log fold change and standard error datasets
df_lfc = data.frame(res$lfc * res$diff_abn, check.names = FALSE) %>% rownames_to_column("taxon_id")
df_se = data.frame(res$se * res$diff_abn, check.names = FALSE) %>% rownames_to_column("taxon_id")
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

#create the dataframe needed to make the plot with log fold change values and standard error joined by taxon
df_fig_soil_type = df_lfc %>% 
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, soil_typerhizosphere) %>%
  dplyr::filter(soil_typerhizosphere != 0) %>%
  dplyr::mutate(direct = ifelse(soil_typerhizosphere > 0, "Positive LFC", "Negative LFC"))
df_fig_soil_type$taxon_id = factor(df_fig_soil_type$taxon_id, levels = df_fig_soil_type$taxon_id)
df_fig_soil_type$direct = factor(df_fig_soil_type$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))

##create figure
p_soil_type = ggplot(data = df_fig_soil_type, aes(x = taxon_id, y = soil_typerhizosphere, fill = direct, color = direct)) +
  geom_bar(stat = "identity", width = 0.7,
           position = position_dodge(width = 0.4)) +
  labs(x = NULL, y = "Log fold change",
       title = "Enriched/Depleted taxa in 2019 rhizosphere vs bulk soil") +
  scale_fill_discrete(name = NULL) + 
  scale_color_discrete(name = NULL) +
  theme_bw() 


##write csv for diff abundance data
write.table(df_fig_soil_type, file = "~/Desktop/mergedredo/12500/soiltypeDA2019_bacteria_genus.txt", sep = "\t", col.names = TRUE, row.names = FALSE)






###significant effect of soil type (bulk vs rhiz) in 2020 woodman samples - differential abundance 
tax_table(wops20)[tax_table(wops20)[,"Genus"]== "","Genus"] <- "Unmatchedgenus"

genus_dat = aggregate_taxa(wops20, "Genus")

out<-ANCOMBC::ancombc(phyloseq = genus_dat, formula = "soil_type", 
                      p_adj_method = "holm", prv_cut = 0.1, 
                      lib_cut = 1000, group = "soil_type",
                      struc_zero = TRUE, neg_lb = TRUE,
                      tol = 1e-5, max_iter = 100,
                      conserve = TRUE, alpha = 0.05, global = TRUE)
res<-out$res
kable(head(res$diff_abn))
lfc<-res$lfc
head(lfc)


df_lfc = data.frame(res$lfc * res$diff_abn, check.names = FALSE) %>% rownames_to_column("taxon_id")
df_se = data.frame(res$se * res$diff_abn, check.names = FALSE) %>% rownames_to_column("taxon_id")
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

df_fig_soil_type = df_lfc %>% 
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, soil_typerhizosphere) %>%
  dplyr::filter(soil_typerhizosphere != 0) %>%
  dplyr::mutate(direct = ifelse(soil_typerhizosphere > 0, "Positive LFC", "Negative LFC"))
df_fig_soil_type$taxon_id = factor(df_fig_soil_type$taxon_id, levels = df_fig_soil_type$taxon_id)
df_fig_soil_type$direct = factor(df_fig_soil_type$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))

p_soil_type = ggplot(data = df_fig_soil_type, aes(x = taxon_id, y = soil_typerhizosphere, fill = direct, color = direct)) +
  geom_bar(stat = "identity", width = 0.7,
           position = position_dodge(width = 0.4)) +
  labs(x = NULL, y = "Log fold change",
       title = "Enriched/Depleted taxa in 2020 rhizosphere vs bulk soil") +
  scale_fill_discrete(name = NULL) + 
  scale_color_discrete(name = NULL) +
  theme_bw() 


##write csv for diff abundance data
write.table(df_fig_soil_type, file = "~/Desktop/mergedredo/12500/soiltypeDA2020_bacteria_genus.txt", sep = "\t", col.names = TRUE, row.names = FALSE)



#### ====================================================================== ####

#### Differences by location (2019 data only)

#### ====================================================================== ####


loc_ps<-subset_samples(ps, site_type != "agricultural")


tax_table(loc_ps)[tax_table(loc_ps)[,"Phylum"]== "","Phylum"] <- "Unmatched_phylum"

library(microbiome)
phylum_dat = aggregate_taxa(loc_ps, "Phylum")

pseq = subset_samples(phylum_dat, location %in% c("Durham_NH", "Rye_NH", "Utah", "Vermont", "Colorado"))
out3<-ANCOMBC::ancombc(phyloseq = pseq, formula = "location", 
                       p_adj_method = "holm", prv_cut = 0.1, 
                       lib_cut = 100, group = "location",
                       struc_zero = TRUE, neg_lb = TRUE,
                       tol = 1e-5, max_iter = 100,
                       conserve = TRUE, alpha = 0.05, global = TRUE)

##got a warning message because Utah sample size is < 5, results are unstable with small sample size
l_res<-out3$res
res_global = out3$res_global #need to get global bc there are more than 2 levels in location

##getting bias-corrected abundances
samp_frac = out3$samp_frac
samp_frac[is.na(samp_frac)] = 0 #replace NA with 0
log_obs_abn = log(abundances(pseq) + 1) # Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn_adj = t(t(log_obs_abn) - samp_frac) # Adjust the log observed abundances


#making a dataframe of the significant taxa and getting the taxonomy into columns
sig_taxa = res_global %>%
  tibble::rownames_to_column("taxon") %>%
  dplyr::filter(diff_abn == TRUE) %>% ##only getting taxa that were differentially abundant
  .$taxon 


#getting a datframe of the metadata and putting sample names as a column
meta_df = meta(phylum_dat) %>%
  tibble::rownames_to_column("sample")


#getting the abundances, adding in metadata, and elongating the dataset
df_sig = as.data.frame(t(log_obs_abn_adj[sig_taxa, ]))%>%
  tibble::rownames_to_column("sample") %>%
  dplyr::left_join(meta(meta_df)) %>%
  dplyr::filter(!is.na(location)) %>%
  tidyr::pivot_longer(cols = -one_of("sample", "SampleID", "location", "input", "filtered", "percentage.of.input.passed.filter", "denoised", "non.chimeric", 
                                     "percentage.of.input.non.chimeric", "soil_type", "location","host_type", "host_sg", "site_type", "texture", "host", "year", 
                                     "ph", "calcium", "magnesium", "potassium", "phosphorus", "lead", "om", "site",
                                     "elevation", "solar_rad", "av_temp", "av_precip", "leaf_wetness_hrs", "rh_90_hrs"), #this is saying to elongate NOT one of these columns - I think you have to do it this double negative way unless you want to use one_of() and list out all of the taxonomy columns
                      names_to = "taxon", values_to = "value")
#creating data frame for heat map
df_heat = df_sig %>%
  dplyr::group_by(location, taxon) %>% #putting location and taxon at the front
  dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) %>% #getting mean of values?
  dplyr::mutate(value = round(value, 2)) %>% #rounding values
  dplyr::arrange(location) 
df_heat$taxon = factor(df_heat$taxon, levels = sig_taxa) #ordering the taxonomy by significant taxa


#Remove the underscores from the beginning of the Factor Labels
df_heat$taxon<-str_sub(df_heat$taxon,6) 


#creating heat map
lo = floor(min(df_heat$value))
up = ceiling(max(df_heat$value))
mid = (lo + up)/2
p_heat = df_heat %>%
  ggplot(aes(x = location, y = taxon, fill = value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(location, taxon, label = value), color = "black", size = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(colour = "black", size = 14)) +
  theme(axis.text.y = element_text(colour = "black", size = 14)) +
  xlab("") +
  ylab("")


###### Differential abundance by texture

tex_ps = subset_samples(loc_ps, texture != "na")

tax_table(tex_ps)[tax_table(tex_ps)[,"Phylum"]== "","Phylum"] <- "Unmatched_phylum"


phylum_dat_tex = aggregate_taxa(tex_ps, "Phylum")
#need to remove silty clay - only 2 instances

pseq_tex = subset_samples(phylum_dat_tex, texture %in% c("clay loam", "loamy sand", "sand", "sandy loam", "loam"))

#perform the differential abundance analyis, formula and group should be the factor you are interested in
# contains: 1) log fold changes; 2) standard errors; 3) test statistics; 
# 4) p-values; 5) adjusted p-values; 6) indicators whether the taxon is differentially abundant (TRUE) or not (FALSE).
out3_tex<-ANCOMBC::ancombc(phyloseq = pseq_tex, formula = "texture", 
                       p_adj_method = "holm", prv_cut = 0.1, #holm is default p-value adjustment, taxa with prevalence less that prv_cut will be excluded, default 0.1
                       lib_cut = 100, group = "texture", #samples with library sizes less than lib cut will be excluded
                       struc_zero = TRUE, neg_lb = TRUE, 
                       tol = 1e-5, max_iter = 100, #tol - convergence tolerance - default 1e-5, max_iter-maximum iterations, default 100
                       conserve = TRUE, alpha = 0.05, global = TRUE) #conserve = TRUE uses a conservative variance estimator for test statistic - recommended for small sample sizes


l_restex<-out3_tex$res
res_globaltex = out3_tex$res_global #need to get global bc there are more than 2 levels in the factor
#res_globaltex has abundance, p val, adjusted p val, and diff_abn

##getting bias-corrected abundances
samp_fractex = out3_tex$samp_frac #getting sample fractions
samp_fractex[is.na(samp_fractex)] = 0 #replace NA with 0
log_obs_abntex = log(abundances(pseq_tex) + 1) # Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn_adjtex = t(t(log_obs_abntex) - samp_fractex) # Adjust the log observed abundances


#making a dataframe of the significant taxa and getting the taxonomy into columns
sig_taxatex = res_globaltex %>%
  tibble::rownames_to_column("taxon") %>%
  dplyr::filter(diff_abn == TRUE) %>% ##only getting taxa that were differentially abundant
  .$taxon 

meta_dftex = meta(pseq_tex) %>%
  tibble::rownames_to_column("sample")

#getting the abundances, adding in metadata, and pivoting sample columns but not metadata colums
df_sigtex = as.data.frame(t(log_obs_abn_adjtex[sig_taxatex, ]))%>%
  tibble::rownames_to_column("sample") %>%
  dplyr::left_join(meta(meta_dftex)) %>%
  dplyr::filter(!is.na(texture)) %>%
  tidyr::pivot_longer(cols = -one_of("sample", "SampleID", "location", "input", "filtered", "percentage.of.input.passed.filter", "denoised", "non.chimeric", 
                                     "percentage.of.input.non.chimeric", "soil_type", "location","host_type", "host_sg", "site_type", "texture", "host", "year", 
                                     "ph", "calcium", "magnesium", "potassium", "phosphorus", "lead", "om", "site",
                                     "elevation", "solar_rad", "av_temp", "av_precip", "leaf_wetness_hrs", "rh_90_hrs"), #this is saying to elongate NOT one of these columns - I think you have to do it this double negative way unless you want to use one_of() and list out all of the taxonomy columns
                      names_to = "taxon", values_to = "value")

#create datframe for heat plot
df_heattex = df_sigtex %>%
  dplyr::group_by(texture, taxon) %>% #putting texture and taxon at the front
  dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) %>% #getting mean of values
  dplyr::mutate(value = round(value, 2)) %>% #rounding values
  dplyr::arrange(texture) 
df_heattex$taxon = factor(df_heattex$taxon, levels = sig_taxatex) #ordering the taxonomy by significant taxa

#Remove the underscores from the beginning of the Factor Labels
df_heattex$taxon<-str_sub(df_heattex$taxon,6) 

#define lo, medium and high levels based on LFC value
lotex = floor(min(df_heattex$value))
uptex = ceiling(max(df_heattex$value))
midtex = (lotex + uptex)/2
p_heat_tex = df_heattex %>%
  ggplot(aes(x = texture, y = taxon, fill = value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = midtex, limit = c(lotex, uptex),
                       name = NULL) +
  geom_text(aes(texture, taxon, label = value), color = "black", size = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(colour = "black", size = 14)) +
  theme(axis.text.y = element_text(colour = "black", size = 14)) +
  xlab("") +
  ylab("")

ggarrange(p_heat, p_heat_tex, nrow = 1, ncol =2, common.legend = TRUE)




#### ====================================================================== ####

#### Differences by site type 

#### ====================================================================== ####

du_ps<-subset_samples(ps, location == "Durham_NH")
rye_ps<-subset_samples(ps, location == "Rye_NH")
nh_ps = merge_phyloseq(du_ps, rye_ps)

nh19_ps<- subset_samples(nh_ps, year.plot == "2019/Plot1")
nh20_ps<- subset_samples(nh_ps, year.plot == "2020/Plot2")


## 2019
tax_table(nh19_ps)[tax_table(nh19_ps)[,"Genus"]== "","Genus"] <- "Unmatchedgenus"

genus_dat19 = aggregate_taxa(nh19_ps, "Genus")

out19<-ANCOMBC::ancombc(phyloseq = genus_dat19, formula = "site_type", 
                      p_adj_method = "holm", prv_cut = 0.1, 
                      lib_cut = 1000, group = "site_type",
                      struc_zero = TRUE, neg_lb = TRUE,
                      tol = 1e-5, max_iter = 100,
                      conserve = TRUE, alpha = 0.05, global = TRUE)
res19<-out19$res
kable(head(res19$diff_abn))
lfc19<-res19$lfc
head(lfc19)


df_lfc19 = data.frame(res19$lfc * res19$diff_abn, check.names = FALSE) %>% rownames_to_column("taxon_id")
df_se19 = data.frame(res19$se * res19$diff_abn, check.names = FALSE) %>% rownames_to_column("taxon_id")
colnames(df_se19)[-1] = paste0(colnames(df_se19)[-1], "SE")

df_fig_site_type19 = df_lfc19 %>% 
  dplyr::left_join(df_se19, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, site_typewild) %>%
  dplyr::filter(site_typewild != 0) %>%
  dplyr::mutate(direct = ifelse(site_typewild > 0, "Positive LFC", "Negative LFC"))
df_fig_site_type19$taxon_id = factor(df_fig_site_type19$taxon_id, levels = df_fig_site_type19$taxon_id)
df_fig_site_type19$direct = factor(df_fig_site_type19$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))

p_site_type19 = ggplot(data = df_fig_site_type19, aes(x = taxon_id, y = site_typewild, fill = direct, color = direct)) +
  geom_bar(stat = "identity", width = 0.7,
           position = position_dodge(width = 0.4)) +
  labs(x = NULL, y = "Log fold change",
       title = "Enriched/Depleted taxa in 2019 wild vs ag soil") +
  scale_fill_discrete(name = NULL) + 
  scale_color_discrete(name = NULL) +
  theme_bw() 


##write csv for diff abundance data
write.table(df_fig_site_type19, file = "~/Desktop/mergedredo/12500/sitetypeDA2019_bacteria_genus.txt", sep = "\t", col.names = TRUE, row.names = FALSE)



## 2020
tax_table(nh20_ps)[tax_table(nh20_ps)[,"Genus"]== "","Genus"] <- "Unmatchedgenus"

genus_dat20 = aggregate_taxa(nh20_ps, "Genus")

out20<-ANCOMBC::ancombc(phyloseq = genus_dat20, formula = "site_type", 
                        p_adj_method = "holm", prv_cut = 0.1, 
                        lib_cut = 1000, group = "site_type",
                        struc_zero = TRUE, neg_lb = TRUE,
                        tol = 1e-5, max_iter = 100,
                        conserve = TRUE, alpha = 0.05, global = TRUE)
res20<-out20$res
kable(head(res20$diff_abn))
lfc20<-res20$lfc
head(lfc20)


df_lfc20 = data.frame(res20$lfc * res20$diff_abn, check.names = FALSE) %>% rownames_to_column("taxon_id")
df_se20 = data.frame(res20$se * res20$diff_abn, check.names = FALSE) %>% rownames_to_column("taxon_id")
colnames(df_se20)[-1] = paste0(colnames(df_se20)[-1], "SE")

df_fig_site_type20 = df_lfc20 %>% 
  dplyr::left_join(df_se20, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, site_typewild) %>%
  dplyr::filter(site_typewild != 0) %>%
  dplyr::mutate(direct = ifelse(site_typewild > 0, "Positive LFC", "Negative LFC"))
df_fig_site_type20$taxon_id = factor(df_fig_site_type20$taxon_id, levels = df_fig_site_type20$taxon_id)
df_fig_site_type20$direct = factor(df_fig_site_type20$direct, 
                                   levels = c("Positive LFC", "Negative LFC"))

p_site_type20 = ggplot(data = df_fig_site_type20, aes(x = taxon_id, y = site_typewild, fill = direct, color = direct)) +
  geom_bar(stat = "identity", width = 0.7,
           position = position_dodge(width = 0.4)) +
  labs(x = NULL, y = "Log fold change",
       title = "Enriched/Depleted taxa in 2019 wild vs ag soil") +
  scale_fill_discrete(name = NULL) + 
  scale_color_discrete(name = NULL) +
  theme_bw() 


##write csv for diff abundance data
write.table(df_fig_site_type20, file = "~/Desktop/mergedredo/12500/sitetypeDA2020_bacteria_genus.txt", sep = "\t", col.names = TRUE, row.names = FALSE)










