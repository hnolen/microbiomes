###DIVERSITY METRICS AND VISUALIZATION###
#load packages
library("vegan")
library("agricolae")

#### rare_data<-rrarefy(abundance.matrix, 3500)

#drop control values
level5<-subset(level5, type != "Control")

#subset bulk and rhizosphere samples
level5_rz <-subset(level5, type == "Rhizosphere")
level5_b <- subset(level5, type == "Bulk")

#first things first: make the first column the rownames
row.names(level5_rz) <- level5_rz$index
level5_rz[1] <- NULL

row.names(level5_b) <- level5_b$index
level5_b[1] <- NULL
#make grow system an ordered variable (not sure why we need to do this)
level3$grow_system <- ordered(level3$grow_system, c("Farm", "Wild", "Control"))

#Question: does species diversity vary across grow systems?
#abundance matrix - getting columns with taxa abundance
#last 13 columns are metadata info
abundance.matrix.rz <- level5_rz[,1:311]

abundance.matrix.b <- level5_b[,1:311]

# store computed indices in a new data frame called 'indices'
indices_rz <- level5_rz[,c("type", "host","location","grow_system")]
indices_rz$Richness <- rowSums(abundance.matrix.rz>0)
indices_rz$Shannon <- diversity(abundance.matrix.rz) # shannon is default
indices_rz$Rarefied <- c(rarefy(abundance.matrix.rz[1:30,], sample=3500, se = FALSE)) #sample is the subsample size that the researcher chose, based on rarefaction curves

indices_b <- level5_b[,c("type", "host","location","grow_system")]
indices_b$Richness <- rowSums(abundance.matrix.b>0)
indices_b$Shannon <- diversity(abundance.matrix.b) # shannon is default
indices_b$Rarefied <- c(rarefy(abundance.matrix.b[1:57,], sample=3500, se = FALSE))


#probelm with some NTC's gettting in: got the error: requested 'sample' was larger than smallest site maximum (0)

par(mar=c(3,4,3,2), mfrow=c(1,2))
colors = terrain.colors(6)[5:1]
boxplot(Richness~host, data=indices, boxwex=0.5, col=colors, 
        cex.axis=0.5, ylab="Richness")
boxplot(Shannon~host, data=indices, boxwex=0.5, col=colors, 
        cex.axis=0.5, ylab="Shannon diversity")
#boxplot(Rarefied~host, data=indices, boxwex=0.5, col=colors, 
#        cex.axis=0.5, ylab="Rarefied richness")


#can test for differences among habitats statistically using a linear model, with habitat as a predictor of diversity
# fit linear models
mod.richness.rz <- lm(Richness~host, data=indices_rz)
mod.shannon.rz <- lm(Shannon~host, data=indices_rz)
mod.rich.grow.rz<-lm(Richness~grow_system, data=indices_rz)
mod.shann.grow.rz<-lm(Shannon~grow_system, data=indices_rz)
mod.rich.loc.rz<-lm(Richness~location, data=indices_rz)
mod.shann.loc.rz<-lm(Shannon~location, data=indices_rz)
#mod.rich.type.rz<-lm(Richness~type, data=indices_rz)
#mod.shann.type.rz<-lm(Shannon~type, data=indices_rz)
#mod.rarefied <- lm(Rarefied~Habitat, data=indices)

mod.richness.b <- lm(Richness~host, data=indices_b)
mod.shannon.b <- lm(Shannon~host, data=indices_b)
mod.rich.grow.b<-lm(Richness~grow_system, data=indices_b)
mod.shann.grow.b<-lm(Shannon~grow_system, data=indices_b)
mod.rich.loc.b<-lm(Richness~location, data=indices_b)
mod.shann.loc.b<-lm(Shannon~location, data=indices_b)
#mod.rich.type.b<-lm(Richness~type, data=indices_b)
#mod.shann.type.b<-lm(Shannon~type, data=indices_b)
#mod.rarefied <- lm(Rarefied~Habitat, data=indices)

# ANOVA 
anova(mod.richness.rz) # NS
anova(mod.shannon.rz) # NS
anova(mod.rich.grow.rz) # p = 0.005723
anova(mod.shann.grow.rz) # p = 0.001098
anova(mod.rich.loc.rz) # NS
anova(mod.shann.loc.rz) # NS
#anova(mod.rich.type) # p = 3.144e-06
#anova(mod.shann.type) #NS
#anova(mod.rarefied) 

anova(mod.richness.b) # NS
anova(mod.shannon.b) # NS
anova(mod.rich.grow.b) # NS
anova(mod.shann.grow.b) # NS
anova(mod.rich.loc.b) # NS
anova(mod.shann.loc.b) # NS

summary(mod.richness)
summary(mod.shannon)

rich_tuk<-HSD.test(mod.richness, "host")
shann_tuk<-HSD.test(mod.shannon, "host")
rich_loc<-HSD.test(mod.rich.loc, "location")
shann_loc<-HSD.test(mod.shann.loc, "location")
rich_grow<-HSD.test(mod.rich.grow.rz, "grow_system")
shann_grow<-HSD.test(mod.shann.grow.rz, "grow_system")
rich_type<-HSD.test(mod.rich.type, "type")
shann_type<-HSD.test(mod.shann.type, "type")

par(mar=c(3,4,3,2), mfrow=c(1,2))#idk what mar= is doing, but mfrow is making rows and columns for multiple graphs in a panel
colors = terrain.colors(6)[5:1]
boxplot(Richness~grow_system, data=indices, boxwex=0.25, col=colors, 
        cex.axis=1.2, main="Species Richness")
boxplot(Shannon~grow_system, data=indices, boxwex=0.25, col=colors, 
        cex.axis=1.2, main="Shannon diversity")
text(x=5.05,y=1.9,labels="*",cex=2)
#boxplot(Rarefied~host, data=indices, boxwex=0.5, col=colors, 
#cex.axis=0.5, ylab="Rarefied richness")

png(file = "diversity~grow.png", width = 500, height = 480, units = px, bg = "white")
dev.off()


par(mar=c(3,4,3,2), mfrow=c(1,2)) #idk what mar= is doing, but mfrow is making rows and columns for multiple graphs in a panel
colors = terrain.colors(6)[5:1]
boxplot(Richness~location, data=indices, boxwex=0.5, col=colors, 
        cex.axis=0.5, ylab="Richness")
boxplot(Shannon~location, data=indices, boxwex=0.5, col=colors, 
        cex.axis=0.5, ylab="Shannon diversity")
#boxplot(Rarefied~host, data=indices, boxwex=0.5, col=colors, 
#cex.axis=0.5, ylab="Rarefied richness")



par(mar=c(3,4,3,2), mfrow=c(1,2)) #idk what mar= is doing, but mfrow is making rows and columns for multiple graphs in a panel
colors = terrain.colors(6)[5:1]
boxplot(Richness~type, data=indices, boxwex=0.5, col=colors, 
        cex.axis=0.5, ylab="Richness")
boxplot(Shannon~type, data=indices, boxwex=0.5, col=colors, 
        cex.axis=0.5, ylab="Shannon diversity")
#boxplot(Rarefied~host, data=indices, boxwex=0.5, col=colors, 
#cex.axis=0.5, ylab="Rarefied richness")

#trying to make this plot in ggplot:
library(ggplot2)
library(ggsignif)
ggplot(indices, aes(x = factor(type), y = Richness)) + 
        geom_boxplot() +
        labs(x = "Sample type", y = "Richness") +
        theme(axis.title.x = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12), axis.text = element_text(size = 10, colour = "black")) +
        geom_signif(comparisons = list(c("Bulk", "Rhizosphere")), map_signif_level = TRUE)

# rhizosphere farm vs wild:
ggplot(indices_rz, aes(x = factor(grow_system), y = Richness)) + 
        geom_boxplot() +
        labs(x = "Grow system", y = "Richness") +
        theme(axis.title.x = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12), axis.text = element_text(size = 10, colour = "black")) +
        geom_signif(comparisons = list(c("Farm", "Wild")), map_signif_level = TRUE)


library("ggplot2")
library("ggsignif")
ggplot(iris, aes(x=Species, y=Sepal.Length)) +
        geom_boxplot() + 
        geom_signif(comparisons = list(c("versicolor", "virginica")),
                    map_signif_level = TRUE)
ggplot(indices, aes(x=grow_system, y=Shannon), na.rm=TRUE) +
        geom_boxplot() +
        geom_signif(comparisons = list(c("cultivated", "wild")),
                    map_signif_level = TRUE)


##repeating not separating by bulk vs rhizosphere:
row.names(level5) <- level5$index
level5[1] <- NULL
abundance.matrix <- level5[,1:311]

indices <- level5[,c("type", "host","location","grow_system")]
indices$Richness <- rowSums(abundance.matrix>0)
indices$Shannon <- diversity(abundance.matrix) # shannon is default
#indices$Rarefied <- c(rarefy(abundance.matrix[1:30,], sample=3500, se = FALSE)) #sample is the subsample size that the researcher chose, based on rarefaction curves


mod.richness <- lm(Richness~host, data=indices)
mod.shannon <- lm(Shannon~host, data=indices)
mod.rich.grow<-lm(Richness~grow_system, data=indices)
mod.shann.grow<-lm(Shannon~grow_system, data=indices)
mod.rich.loc<-lm(Richness~location, data=indices)
mod.shann.loc<-lm(Shannon~location, data=indices)
mod.rich.type<-lm(Richness~type, data=indices)
mod.shann.type<-lm(Shannon~type, data=indices)
#mod.rarefied <- lm(Rarefied~Habitat, data=indices)

anova(mod.richness) # NS
anova(mod.shannon) # NS
anova(mod.rich.grow) # NS
anova(mod.shann.grow) # NS 
anova(mod.rich.loc) # NS 
anova(mod.shann.loc) # NS
anova(mod.rich.type) # p = 3.144e-06
anova(mod.shann.type) # NS

ggplot(indices, aes(x = factor(type), y = Richness)) + 
        geom_boxplot() +
        labs(x = "Sample type", y = "Richness") +
        theme(axis.title.x = element_text(color = "black", size = 20), axis.title.y = element_text(color = "black", size = 20), axis.text = element_text(size = 20, colour = "black")) +
        geom_signif(comparisons = list(c("Bulk", "Rhizosphere")), map_signif_level = TRUE, textsize=7)
