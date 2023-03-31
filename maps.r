#make Palau shapefile as a ggplot2 object and plot collection coordinates
#note: the mapdata package that is compatible with phytools function phylo.to.map doesn't have a great map of Palau 
#the phylo.to.map palau map is superimposed to this map by lining the points from the phylo.to.map output with the points from the ggplot2 map

#load libraries
library(rgdal)
library(ggplot2)
library(ggthemes)

#create a shape file for Palau that is R readable 
PalauMap<- st_read("~/Desktop/Microarthropod/graphs_R/Palau_Shoreline.shp")
#fortify map for ggplot2
fortPalau_df <- fortify(PalauMap)
#read coordinate file
jitterAdj=read.table("~/Desktop/Arachnid/Locality Data/jitterCordNewNew.txt",header=T,sep="\t")
#plot coordinates on Palau map
ggplot(data=jitterAdj, aes(x=lat, y=long))+ 
geom_path(data = fortPalau_df, aes(x = long, y = lat, group = group),color = 'black', size = .2)+
geom_point(size=1,shape=21,aes(color="black",fill="#368180"))+
theme_void()+theme(legend.position="none")
#save
ggsave("palauCoord_black.png",plot=last_plot(),bg="transparent")


#load tree files and plot with phylo.to.map 

#load libraries
library(ape)
library(castor)
library(treeio)
library(ggtree)
library(phytools)
library(mapdata)

#read support tree (generated in RAxML)
phylo <- read.tree(file = "~/Desktop/mar23_sequenceFiles/allCOI_sulawesi_trimmed.tree.raxml.support")
#collapse tree at 0.05 distance 
collapsed = collapse_tree_at_resolution(phylo, resolution=0.05,rename_collapsed_nodes=TRUE)
#view tree 
ggtree(collapsed$tree) + geom_text(aes(label=node),size=1, hjust=-.3) + geom_tiplab() + xlim(0,.5)
#look at the Palau node of the full tree 
full_palauClade <- tree_subset(phylo, node=343, levels_back=0)
ggtree(full_palauClade) + geom_text(aes(label=node),size=.5, hjust=-.1)+ geom_tiplab(size=1) + xlim(0, .25)
#check which nodes were collapsed - useful for assigning top values to all coordinate points for phylo to map 
collapsed[["new2old_clade"]]
#view the collapsed nodes to a get list of taxa included in each node on the new tree 
#make ggtree object from tree 
p=ggtree(phylo)
#exract clades from ggtree object
viewClade(p + geom_text(aes(label=node), size=1,hjust=-.1)+geom_tiplab(), node=344)
#repeat for clades 1-52
#write newick file for the collapsed tree (at .05 distance)
write.tree(collapsed$tree,file="COI_.05distSubtree_allLocale_guatemala")

#extract subtree from Palau clade - we know it is node 80 from our previous tree above 
palauClade <- tree_subset(collapsed$tree, node=52, levels_back=0)
#generate tree
ggtree(palauClade) + geom_text(aes(label=node), hjust=-.3)+ geom_tiplab() + xlim(0, .5)
#make newick file for the collapsed tree (at .05 distance)
write.tree(palauClade,file="COI_.05distSubtree_guatemala_palauClade")

#phylo to map on these 0.05 distance pruned tree 
#load locality data (All specimens)
coordALl_5dist =read.table("~/Desktop/rDirectory/schizomida/newSulawesi_coords.txt", header = T, sep = '\t')
#matrix
coordALl_5dist=as.matrix(coordALl_5dist)
#extract rows, make into numeric object
latMat05 = coordALl_5dist[,2 ]
latMat05=as.numeric(latMat05)
longMat05 = coordALl_5dist[ ,3]
longMat05=as.numeric(longMat05)
#extract taxa names
taxaNames = coordALl_5dist[ , 1]
#write to new dataframe
coorMat_5dist=cbind(latMat05,longMat05)
colnames(coorMat_5dist) <- c("lat", "long")
rownames(coorMat_5dist) <- taxaNames
#make phylo to map object
obj<-phylo.to.map(collapsed$tree,coorMat_5dist,database="worldHires",xlim=c(90, 165),ylim=c(-40, 18),plot=FALSE,direction="rightwards")
#make color palette
cols<-setNames(sample(met.brewer("Ingres",n=Ntip(collapsed$tree))),
    collapsed$tree$tip.label)
#plot
plot(obj,direction="rightwards",colors=cols,lty="solid",ftype="off",cex.points=c(1,.75))


#plotting the Palau-only clade
#load locality data 
palauClade_coord =read.table("~/Desktop/rDirectory/schizomida/newSulawesi_palauClade_coords.txt", header = T, sep = '\t')
#matrix
palauClade_coord=as.matrix(palauClade_coord)
#extract rows, make into numeric object
latMat_p5 = palauClade_coord[,2 ]
latMat_p5=as.numeric(latMat_p5)
longMat_p5 = palauClade_coord[ ,3]
longMat_p5=as.numeric(longMat_p5)
#extract taxa names
taxaNames_p5 = palauClade_coord[ , 1]
#write to new dataframe
coordMat_p5=cbind(latMat_p5,longMat_p5)
colnames(coordMat_p5) <- c("lat", "long")
rownames(coordMat_p5) <- taxaNames_p5
#make phylo to map object
objP5<-phylo.to.map(palauClade,coordMat_p5,database="worldHires",xlim=c(90, 165),ylim=c(-40, 18),plot=FALSE)
#make color palette
colsP5<-setNames(sample(met.brewer("Navajo",n=Ntip(palauClade))),
    palauClade$tip.label)
#plot
plot(objP5,direction="rightwards",colors=colsP5,lty="solid",ftype="off",cex.points=c(1,.75))

#Palau map  
onlyPalau_coord =read.table("~/Desktop/rDirectory/schizomida/newSulawesi_palauONLY_coords.txt", header = T, sep = '\t')
#matrix
onlyPalau_coord=as.matrix(onlyPalau_coord)
#extract rows, make into numeric object
latMat_palau = onlyPalau_coord[,2 ]
latMat_palau=as.numeric(latMat_palau)
longMat_palau = onlyPalau_coord[ ,3]
longMat_palau=as.numeric(longMat_palau)
#extract taxa names
taxaNames_palau = onlyPalau_coord[ , 1]
#write to new dataframe
coordMat_palau=cbind(latMat_palau,longMat_palau)
colnames(coordMat_palau) <- c("lat", "long")
rownames(coordMat_palau) <- taxaNames_palau
#make phylo to map object
palauObj<-phylo.to.map(palauClade,coordMat_palau,database="worldHires",xlim=c(130, 135),ylim=c(6, 10),plot=FALSE)

#plot
plot(palauObj,direction="rightwards",colors=colsP5,lty="solid",ftype="off",cex.points=c(1,.75))

#change colors for both the pruned tree and the Palau clade subtree based on bioregion, rather than tip label 
#EDBF82 Carolines
#57643C Philippines
#368180 Palau
#D56234 Yap
#10687E Sulawesi 
#A29653 Sahul
#C48231 Guam
#AC3414 Sundaland
#6D6459 Guatemala 

cols[["JCM3039s001_Palau_Ngerkeklau"]] <- "#368180"
cols[["JCM3517s004_Palau_Mecherchar"]] <- "#368180"
cols[["CASENT_9047511_Philippines"]] <- "#57643C"
cols[["JCM3520s003_Palau_Babeldaob_Airai"]] <- "#368180"
cols[["AIRAIs002_Palau_Babeldaob_Airai"]] <- "#368180"
cols[["JCM3051s002_Palau_Ulong"]] <- "#368180"
cols[["JCM3036s002_Palau_Ngesomel"]] <- "#368180"
cols[["809.1_Chuuk"]] <- "#EDBF82"
cols[["818.3_Pohnpei"]] <- "#EDBF82"
cols[["NGS27_Philippines"]] <- "#57643C"
cols[["JCM0262.1_Palau_Merir"]] <- "#368180"
cols[["asv_lco_19"]] <- "#10687E"
cols[["MCZ_IZ-135128_Brunei"]] <- "#AC3414"
cols[["MJWP08.0237_PapuaNewGuinea"]] <- "#A29653"
cols[["NGS26_Philippines"]] <- "#57643C"
cols[["WAM_T120190_Austrailia"]] <- "#A29653"
cols[["WAM_T122128_Austrailia"]] <- "#A29653"
cols[["WAM_T122887_Austrailia"]] <- "#A29653"
cols[["WAM_T122047_Malaysia"]] <- "#AC3414"
cols[["MCZ_IZ_132985_Indonesia"]] <- "#10687E"
cols[["WAM_T122050_Malaysia"]] <- "#AC3414"
cols[["WAM_T122042_Austrailia"]] <- "#A29653"
cols[["WAM_T122056_Malaysia"]] <- "#AC3414"
cols[["MCZ_IZ_89506_Guatemala"]] <- "#6D6459"
cols[["JCM3091s002_Palau_Babeldaob_Ngiwal"]] <- "#368180"
cols[["845.5_Yap"]] <- "#D56234"
cols[["JCM3047s001_Palau_Ngeruktabel"]] <- "#368180"
cols[["JCM3037s001_Palau_Bablomekang"]] <- "#368180"
cols[["JCM3046s001_Palau_unnamedNgeruktabelGroup"]] <- "#368180"
cols[["840.1_Yap"]] <- "#D56234"
cols[["813.1_Chuuk"]] <- "#EDBF82"
cols[["MCZ_IZ-132988_Philippines"]] <- "#57643C"
cols[["NGS49_Philippines"]] <- "#57643C"
cols[["WAM_T122048_Malaysia"]] <- "#AC3414"
cols[["106356_Indonesia"]] <- "#10687E"
cols[["NGS53_Philippines"]] <- "#57643C"
cols[["853.14_Guam"]] <- "#C48231"
cols[["WAM_T120206_Austrailia"]] <- "#A29653"
cols[["WAM_T120140_Austrailia"]] <- "#A29653"
cols[["WAM_T122124_Austrailia"]] <- "#A29653"
cols[["AMNH_LP_3422_Austrailia"]] <- "#A29653"
cols[["WAM_T122132_Austrailia"]] <- "#A29653"
cols[["MCZ_IZ-132983_Indonesia"]] <- "#AC3414"
cols[["WAM_T122051_Malaysia"]] <- "#AC3414"
cols[["WAM_T122054_Malaysia"]] <- "#AC3414"
cols[["WAM_T122052_Malaysia"]] <- "#AC3414"
cols[["WAM_T122049_Malaysia"]] <- "#AC3414"

colsP5[["JCM3039s001_Palau_Ngerkeklau"]] <- "#368180"
colsP5[["JCM3091s002_Palau_Babeldaob_Ngiwal"]] <- "#368180"
colsP5[["JCM3517s004_Palau_Mecherchar"]] <- "#368180"
colsP5[["845.5_Yap"]] <- "#D56234"
colsP5[["CASENT_9047511_Philippines"]] <- "#57643C"
colsP5[["JCM3047s001_Palau_Ngeruktabel"]] <- "#368180"
colsP5[["JCM3520s003_Palau_Babeldaob_Airai"]] <- "#368180"
colsP5[["JCM3037s001_Palau_Bablomekang"]] <- "#368180"
colsP5[["AIRAIs002_Palau_Babeldaob_Airai"]] <- "#368180"
colsP5[["JCM3046s001_Palau_unnamedNgeruktabelGroup"]] <- "#368180"
colsP5[["JCM3051s002_Palau_Ulong"]] <- "#368180"
colsP5[["840.1_Yap"]] <- "#D56234"
colsP5[["JCM3036s002_Palau_Ngesomel"]] <- "#368180"
colsP5[["813.1_Chuuk"]] <- "#EDBF82"
colsP5[["809.1_Chuuk"]] <- "#EDBF82"
colsP5[["MCZ_IZ-132988_Philippines"]] <- "#57643C"
colsP5[["818.3_Pohnpei"]] <- "#EDBF82"
