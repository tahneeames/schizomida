#alignment, tree generation, and OTU identification 
#mostly in R but the RAxML code is in here, too 

#some data wrangling to quickly identify sequences that needed some manouvering to align properly 
#took mafft alignment and removed sequences with gaps 
#now extracting the sequences I removed in geneious from r to re-add them to alignment file 
library(seqinr)
# Read in the full FASTA file with all COI sequences
all_seqs <- read.fasta("~/Desktop/Arachnid/mar23_sequenceFiles/allSequences_coi.fasta")
# Read in the FASTA file with the missing COI sequences
aligned_seqs <- read.fasta("~/Desktop/Arachnid/mar23_sequenceFiles/noGaps_seq.fasta")
# Extract the IDs of the sequences that are not shared 
missing_ids <- setdiff(names(all_seqs), aligned_ids)
# Subset the full FASTA file by the non-missing IDs
missing_seqs <- all_seqs[missing_ids]
# Write the non-missing sequences to a new FASTA file
write.fasta(missing_seqs, file = "missing_seqs.fasta",names = names(missing_seqs))

#same but again with the new combined fasta
more_seqs <- read.fasta("~/Desktop/Arachnid/mar23_sequenceFiles/newFasta_mostSeqs.fasta")
more_ids <- names(more_seqs)
other_ids <- setdiff(names(all_seqs), more_ids)
moreMissing_seqs <- all_seqs[other_ids]
write.fasta(moreMissing_seqs, file = "moreMissing_seqs.fasta",names = names(moreMissing_seqs))

#18s now - but same process 
rdna_seqs <- read.fasta("~/Desktop/Arachnid/Sequence Results/FASTA files Organized/allSequences_18s.fasta")
some_rdna <- read.fasta("~/Desktop/Arachnid/mar23_sequenceFiles/first_18sAlignment.fasta")
rdna_ids <- names(some_rdna)
rdna_missing_ids <- setdiff(names(rdna_seqs), rdna_ids)
missing_rdnaSeqs <- rdna_seqs[rdna_missing_ids]
write.fasta(missing_rdnaSeqs, file = "missing_18sSEQS",names = names(missing_rdnaSeqs))

second_rdna <- read.fasta("~/Desktop/Arachnid/mar23_sequenceFiles/secondPass_18s_alignmentTrim.fasta")
rdna_ids2 <- names(second_rdna)
rdna_missing_ids2 <- setdiff(names(rdna_seqs), rdna_ids2)
missing_rdnaSeqs2 <- rdna_seqs[rdna_missing_ids2]
write.fasta(missing_rdnaSeqs2, file = "missing_18sSEQS_second",names = names(missing_rdnaSeqs2))

third_rdna <- read.fasta("~/Desktop/Arachnid/mar23_sequenceFiles/18s_thirdAlignment.fasta")
rdna_ids3 <- names(third_rdna)
rdna_missing_ids3 <- setdiff(names(rdna_seqs), rdna_ids3)
missing_rdnaSeqs3 <- rdna_seqs[rdna_missing_ids3]
write.fasta(missing_rdnaSeqs3, file = "missing_18sSEQS_third",names = names(missing_rdnaSeqs3))

fourth_rdna <- read.fasta("~/Desktop/Arachnid/mar23_sequenceFiles/18s_fourthTrimmed.fasta")
rdna_ids4 <- names(fourth_rdna)
rdna_missing_ids4 <- setdiff(names(rdna_seqs), rdna_ids4)
missing_rdnaSeqs4 <- rdna_seqs[rdna_missing_ids4]
write.fasta(missing_rdnaSeqs4, file = "missing_18sSEQS_fourth",names = names(missing_rdnaSeqs4))

#now that the alignment is looking good, we can make a tree with the COI and 18s markers 
#RAxML code command line 
raxml-ng --all --msa /Desktop/Arachnid/mar23_sequenceFiles/trimmed_may23/RAxML_coi_18s_gtrRun_concat.txt --model /Desktop/Arachnid/mar23_sequenceFiles/trimmed_may23/RAxML_coi_18s_gtrRun_concat.part.txt --prefix /Desktop/Arachnid/mar23_sequenceFiles/trimmed_may23/coi_18s_gtrRun --seed 292524 --outgroup MCZ_IZ_89506|Guatemala --bs-metric tbe --tree rand{1} --bs-trees 100
#partition file 
GTR+FO+I+G4m+B,       0_part_1_codon1 = 1-681\3
GTR+FO+I+G4m+B,       0_part_1_codon2 = 2-681\3
GTR+FO+I+G4m+B,       0_part_1_codon3 = 3-681\3
K80+G4m+B,      1_part_1 = 682-1613


#in R - after RAxML run, extract distance matrix for ASAP 
require(ape)
require(phangorn)
#read RAxML tree file 
raxml_concatTree <- read.tree("~/Desktop/Arachnid/mar23_sequenceFiles/trimmed_may23/coi_18s_gtrRun.raxml.bestTree.tre")
#create distance matrix                                 
distMatrix <- cophenetic.phylo(raxml_concatTree) 
#save                                      
writeDist(distMatrix, "RAxML_dnadist_for_ABGD_ASAP.txt", format ="phylip")  


#after running the distance file in the ASAP webtool 
#https://bioinfo.mnhn.fr/abi/public/asap/ | phylip dnadist | seq length 681 
#in R - extract OTUs identified in ASAP
library(ape)
library(ggtree)
#load new beast file 
fullTree_ab44=read.nexus("~/Desktop/Arachnid/mar23_sequenceFiles/trimmed_may23/coi_beast_abPriors_summary.txt")
#make a list of tips we want to keep for biogeobears and figures (0.068 distance)
44tips=c("WAM_T120210_Australia","WAM_T122041_Australia","MJWP08.0237_PapuaNewGuinea","853.15_Guam","MCZ_IZ_135129_PapuaNewGuinea","842.1_Yap","asv_lco_348","848.5_Palau_Babeldaob_Melekeok","JCM3021s002_Palau_Ulong","813.1_Chuuk","809.1_Chuuk","819.2_Pohnpei","JCM0127.1_Palau_Kayangel_Kayangel","845.2_Yap","JCM3520s008_Palau_Babeldaob_Airai","JCM3047s001_Palau_Ngeruktabel","WAM_T122045_Malaysia","106356_Indonesia","JCM0262.1_Palau_Merir","asv_lco_19","NGS55_Philippines","MCZ_IZ_135128_Brunei","MCZ_IZ_132985_Indonesia","WAM_T122129_Australia","WAM_T122117_Australia","AMNH_LP_3422_Australia","WAM_T122887_Australia","WAM_T122056_Malaysia","WAM_T122052_Malaysia","WAM_T120186_Australia","WAM_T120148_Australia","NGS27_Philippines","NGS26_Philippines","WAM_T122132_Australia","MCZ_IZ_132988_Philippines","WAM_T122047_Malaysia","MCZ_IZ_132983_Indonesia","CASENT_9047511_Philippines","WAM_T122054_Malaysia","WAM_T122044_Malaysia","WAM_T122051_Malaysia","WAM_T122050_Malaysia","MCZ_IZ_89506_Guatemala")
#make subtree using tips list
bgb_subTree_44=keep.tip(fullTree_ab44, 44tips)
#save 
write.tree(bgb_subTree_44,file="abrams2019_newBeast_44tips")
#same but for the micronesia beast file 
fullTree_ab44_mm=read.nexus("~/Desktop/Arachnid/mar23_sequenceFiles/trimmed_may23/coi_beast_abPriors_mm_summary.txt")
bgb_subTree_44_mm=keep.tip(fullTree_ab44_mm, 44tips)
write.tree(bgb_subTree_44_mm,file="abrams2019_newBeast_44tips_microMono")
