######STACKED BAR PLOTS AND BUBBLE PLOTS######
# https://jkzorz.github.io/2019/06/05/stacked-bar-plots.html
# https://jkzorz.github.io/2019/06/05/Bubble-plots.html 

#load packages
library(ggplot2)
library(reshape2)

#first things first: make the first column the rownames
#row.names(level4) <- level4$index
#level4[1] <- NULL

#change data format from "wide" to "long"
lev2_long<-melt(level2, id = c("index", "type", "location", "host", "grow_system", "year", "host_sum", "grow_sum", "loc_sum", "type_sum"))
#make sure the abundance value is a numeric
str(lev2_long$value)
lev2_long$value <-as.numeric(lev2_long$value, na.rm = TRUE)

##this section is to just get sums to populate the columns in the file to calculate relative abundance
#getting sums of different groups for relative abundance calculation
host_df<-aggregate(lev2_long$value, by=list(host=lev2_long$host), FUN=sum)
host_sum <- sum(host_df$x)
sum_all<- sum(lev2_long$value)

album_df<-subset(lev2_long, host == "C. album")
album_df$sum<-sum(album_df$value)

quinoa_df<-subset(lev2_long, host == "C. quinoa")
quinoa_df$sum<-sum(quinoa_df$value)

ficifolium_df<-subset(lev2_long, host == "C. ficifolium")
ficifolium_df$sum<-sum(ficifolium_df$value)

bvmquinoa_df<-subset(lev2_long, host == "C. BVM x C. quinoa")
bvmquinoa_df$sum<-sum(bvmquinoa_df$value)

bvm_df<-subset(lev2_long, host == "C. berlandieri var. macrocalycium")
bvm_df$sum<-sum(bvm_df$value)

chenopod_df<-subset(lev2_long, host == "Chenopodium spp.")
chenopod_df$sum<-sum(chenopod_df$value)

control_df<-subset(lev2_long, host == "Control")
control_df$sum<-sum(control_df$value)

##put all of the sums into the csv "host_sum" column- need to figure out how to do this in R

##now getting sums to group by grow_system in the stacked bar plots
farm_df<-subset(lev2_long, grow_system == "Farm")
farm_df$sum<-sum(farm_df$value)

wild_df<-subset(lev2_long, grow_system == "Wild")
wild_df$sum<-sum(wild_df$value)

growcontrol_df<-subset(lev2_long, grow_system == "Control")
growcontrol_df$sum<-sum(growcontrol_df$value)

##getting sum by location
loccontrol_df<-subset(lev2_long, location == "Control")
loccontrol_df$sum<-sum(loccontrol_df$value)

nh_df<-subset(lev2_long, location == "New Hampshire")
nh_df$sum<-sum(nh_df$value)

co_df<-subset(lev2_long, location == "Colorado")
co_df$sum<-sum(co_df$value)

vt_df<-subset(lev2_long, location == "Vermont")
vt_df$sum<-sum(vt_df$value)

##sample type sums
bulk_df<-subset(lev2_long, type == "Bulk")
bulk_df$sum<-sum(bulk_df$value)

rhiz_df<-subset(lev2_long, type == "Rhizosphere")
rhiz_df$sum<-sum(rhiz_df$value)

ctrl_df<-subset(lev2_long, type == "Control")
ctrl_df$sum<-sum(ctrl_df$value)


#transform abundances to relative abundance
lev2_long$rel_abundance<- 100 * lev2_long$value/lev2_long$host_sum
lev2_long$rel_grow<- 100 * lev2_long$value/lev2_long$grow_sum
lev2_long$rel_loc<- 100 * lev2_long$value/lev2_long$loc_sum
lev2_long$rel_type<- 100 * lev2_long$value/lev2_long$type_sum

#if you want to keep the order of your samples run the following code
#otherwise R will arrange them alphabetically
lev2_long$index<-factor(lev2_long$index, levels=unique(lev2_long$index))

#drop control values
lev2_long_noctrl<-subset(lev2_long, type != "Control")


#stacked bar plot notes:
#two types of bar charts: geom_bar() and geom_col(), the bar makes the height of the bar 
#proportional to the # of cases in each group
#geom_col makes the heights of the bars represent value sin the data
#make the plot!
ggplot(lev2_long_noctrl, aes(x = host, fill = variable, y = rel_abundance)) + 
  geom_bar(stat = "identity") +
  #guides(fill=FALSE) +
  theme(axis.text = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 14, colour = "black"), axis.text.x = element_text(size = 14, angle = 45, hjust = 1)) + 
  scale_y_continuous(limits=c(0,100), expand = c(0,0)) + 
  labs(x = "", y = "Relative Abundance (%)", fill = "ASV") +
  #scale_fill_manual(values = colours)
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey"))  
  #theme(legend.position="bottom", strip.background = element_blank(), strip.placement = "outside", panel.spacing = unit(1, "lines")) + 
  #guides(fill=guide_legend(nrow=5))


###bar plot grouped by grow system
ggplot(lev2_long, aes(x = grow_system, fill = variable, y = rel_grow)) + 
  geom_bar(stat = "identity") +
  #guides(fill=FALSE) +
  theme(axis.text = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 14, colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative Abundance (%)", fill = "ASV") +
  #scale_fill_manual(values = colours)
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey"))  
#theme(legend.position="bottom", strip.background = element_blank(), strip.placement = "outside", panel.spacing = unit(1, "lines")) + 
#guides(fill=guide_legend(nrow=5))



###bar plot grouped by location
ggplot(lev2_long, aes(x = location, fill = variable, y = rel_loc)) + 
  geom_bar(stat = "identity") +
  #guides(fill=FALSE) +
  theme(axis.text = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 14, colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative Abundance (%)", fill = "ASV") +
  #scale_fill_manual(values = colours)
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey"))  
#theme(legend.position="bottom", strip.background = element_blank(), strip.placement = "outside", panel.spacing = unit(1, "lines")) + 
#guides(fill=guide_legend(nrow=5))



###bar plot grouped by location
ggplot(lev2_long, aes(x = type, fill = variable, y = rel_type)) + 
  geom_bar(stat = "identity") +
  #guides(fill=FALSE) +
  theme(axis.text = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 14, colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative Abundance (%)", fill = "ASV") +
  #scale_fill_manual(values = colours)
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey"))  
#theme(legend.position="bottom", strip.background = element_blank(), strip.placement = "outside", panel.spacing = unit(1, "lines")) + 
#guides(fill=guide_legend(nrow=5))



#bubble plot - BY SAMPLE
ggplot(lev2_long, aes(x = index, y = variable)) + 
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) + 
  #next line breaks the abudance scale into groups of 1, 10, 50 75 % abundance
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,17), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
  #guides(fill=FALSE) makes sure there isn't a legend containing each taxon
  guides(fill=FALSE)
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right")  
  #scale_fill_manual(values = colours, guide = FALSE) + 
  #scale_y_discrete(limits = rev(levels(lev4_long$variable))) 


lev2_long$location<-factor(lev2_long$location, levels=unique(lev2_long$location))
  
#bubble plot - LOCATION
ggplot(lev2_long, aes(x = location, y = variable)) + 
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) + 
  #next line breaks the abudance scale into groups of 1, 10, 50 75 % abundance
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,17), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
  #guides(fill=FALSE) makes sure there isn't a legend containing each taxon
  guides(fill=FALSE)
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right")  


lev3_long$host<-factor(lev3_long$host, levels=unique(lev3_long$host))
  
#bubble plot - HOST
ggplot(lev3_long, aes(x = host, y = variable)) + 
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) + 
  #next line breaks the abudance scale into groups of 1, 10, 50 75 % abundance
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,17), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
  #guides(fill=FALSE) makes sure there isn't a legend containing each taxon
  guides(fill=FALSE)
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right")    

  
  
lev3_long$grow_system<-factor(lev3_long$grow_system, levels=unique(lev3_long$grow_system))
  
#bubble plot - GROW SYSTEM
ggplot(lev3_long, aes(x = grow_system, y = variable)) + 
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) + 
  #next line breaks the abudance scale into groups of 1, 10, 50 75 % abundance
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,17), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
  #guides(fill=FALSE) makes sure there isn't a legend containing each taxon
  guides(fill=FALSE)
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right")      
  

  

    
  
