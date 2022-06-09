#####INDICATOR SPECIES ANALYSIS#####
#https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html

#install indicspecies package
install.packages("indicspecies")
library("indicspecies")

#first things first: make the first column the rownames
row.names(level3) <- level3$index
level3[1] <- NULL

#make ASV dataframe - getting columns - change final number to however many columns you have
abundance.matrix <- level3[,1:52]
#make vector of the treatment variable you are interested in
grow_sys<-level3$grow_system
host<-level3$host
location<-level3$location

#running the indicator species command
#multipatt returns a list of species that are associated with a particular group of samples
#if your group has more than 2 categories, multipatt will also ID species that are
#statistically more abundant in combinations of categories
#parameters: abudance matrix, vector of variable (grouping info), the function that
#multipatt uses to find indicator species "r.g" and the number of permutations used
#in the statistcal test
inv <- multipatt(abundance.matrix, grow_sys, func = "r.g", control = how(nperm=9999))
inv_host<- multipatt(abundance.matrix, host, func = "r.g", control = how(nperm=9999))
inv_loc<- multipatt(abundance.matrix, location, func = "r.g", control = how(nperm=9999))
#view results
summary(inv)
summary(inv_host)
summary(inv_loc)




##
row.names(level4) <- level4$index
level4[1] <- NULL

#make ASV dataframe - getting columns - change final number to however many columns you have
abundance.matrix <- level4[,1:120]
#make vector of the treatment variable you are interested in
grow_sys<-level4$grow_system
host<-level4$host
location<-level4$location

#running the indicator species command
#multipatt returns a list of species that are associated with a particular group of samples
#if your group has more than 2 categories, multipatt will also ID species that are
#statistically more abundant in combinations of categories
#parameters: abudance matrix, vector of variable (grouping info), the function that
#multipatt uses to find indicator species "r.g" and the number of permutations used
#in the statistcal test
inv4 <- multipatt(abundance.matrix, grow_sys, func = "r.g", control = how(nperm=9999))
inv4_host<- multipatt(abundance.matrix, host, func = "r.g", control = how(nperm=9999))
inv4_loc<- multipatt(abundance.matrix, location, func = "r.g", control = how(nperm=9999))
#view results
summary(inv4)
summary(inv4_host)
summary(inv4_loc)




###level 6 - genus level
row.names(level6) <- level6$index
level6[1] <- NULL

#make ASV dataframe - getting columns - change final number to however many columns you have
abundance.matrix <- level6[,1:682]
#make vector of the treatment variable you are interested in
grow_sys<-level6$grow_system
host<-level6$host
location<-level6$location

#running the indicator species command
#multipatt returns a list of species that are associated with a particular group of samples
#if your group has more than 2 categories, multipatt will also ID species that are
#statistically more abundant in combinations of categories
#parameters: abudance matrix, vector of variable (grouping info), the function that
#multipatt uses to find indicator species "r.g" and the number of permutations used
#in the statistcal test
inv6 <- multipatt(abundance.matrix, grow_sys, func = "r.g", control = how(nperm=9999))
inv6_host<- multipatt(abundance.matrix, host, func = "r.g", control = how(nperm=9999))
inv6_loc<- multipatt(abundance.matrix, location, func = "r.g", control = how(nperm=9999))
#view results
summary(inv6)
summary(inv6_host)
summary(inv6_loc)


