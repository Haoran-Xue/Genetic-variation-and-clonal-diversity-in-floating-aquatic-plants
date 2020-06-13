#**************************************************************************************************************************#
#						Plotting a simple dendrogram to see samples' structure								   #
#**************************************************************************************************************************#

install.packages("ggplot2")
install.packages("dendextend")
install.packages("Rmisc")
install.packages("dplyr")
install.packages("ape")
install.packages("phangorn")

library(ggplot2)
library(dendextendRcpp)
library(dendextend)
library(Rmisc)
library(dplyr)
library(ape)
library(phangorn)

# esse arquivo dist ? matriz de distancias completa com as 285 amostras para o mesmo filtro
# reading the file, notice the check.names because of enconding differences from the server
# and mac, further, I did set the the samples as row.names

azurea.distance1.pl30 <- read.table("distance_PVE_azurea_PGAbfp_HC_50miss_indmiss_qualOdp_fqbd_noindels.distance1Table.DpPl30.tbl", 
                                             header = T, check.names=FALSE, row.names = 1)
# just checking file sctructure and else meta.data
str(azurea.distance1.pl30)
head(azurea.distance1.pl30)
range(azurea.distance1.pl30)
dim(azurea.distance1.pl30)

# minDp30maxDp250Pl30Ind60[upper.tri(minDp30maxDp250Pl30Ind60)] <- "NA" 
# now I do convert it into a distance object
D1_azurea.distance1.pl30 <- as.dist(azurea.distance1.pl30)
D1_azurea.distance1.pl30 %>% head
library(cluster)
D1_azurea.distance1.pl30 %>% sizeDiss()
# Calculate the dendrogram
dend1_azurea.dist1 <- as.dendrogram(hclust(D1_azurea.distance1.pl30))
dend1_azurea.dist1 <- set(dend1_azurea.dist1, "labels_cex", 1)

keys <- read.table ("all_keys_mac.txt", sep = "" , header = T, check.names=FALSE)
keysN <- keys
str(keysN)
rownames(keysN)  <- keys[, 1]

keysN_ordered <- keysN[match(labels(dend1_azurea.dist1), row.names(keysN)), ] # notice that now, I am ordering the fucking keys by the dend using match
cbind.data.frame(labels(dend1_azurea.dist1), row.names(keysN_ordered)) # and it is in the same order wohoo...

# Now plotting the dendrogram
# first setting the fram aspects
par(cex=0.4, mar=c(5, 8, 4, 10))

# Then set the colors for each basin and order it accordint to the dendrogram
colorCodesBasin <- c("black", "red", "blue", "green", "grey", "purple") [keysN_ordered$Basin] # [keysN[match(labels(dend1_new), row.names(keysN)), ]$Basin]
labels_colors(dend1_azurea.dist1) <- colorCodesBasin #[order.dendrogram(dend1_new)]

# plot the dendro
png(file="dend1_azurea_dist1_color.png", width=1200, height=2000)
plot(dend1_azurea.dist1, horiz = TRUE, main = "Dendrogram for E. azurea samples", cex = 1)
legend(0.12, 60, legend=c("Araguaia", "Madeira", "Pantanal",  "Parana", "Teles Pires", "Tocantins"), 
       col=c("black", "red", "blue", "green", "grey", "purple"), pch=19, cex=1.5, box.lty=1)
dev.off()
