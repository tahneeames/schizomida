#using compute canada software suite to align COI sequence results and build the base phylogenetic tree with all sequences 
#Unix 

#load software tools 
module load CVMFS_CC
module load muscle 
module load trimal 

#read in fasta sequence file 
nano allCOI_sulawesi.fasta

#align in muscie 
muscle -in allCOI_sulawesi.fasta -out allCOI_sulawesi.aln

#trim 
trimal -in allCOI_sulawesi.aln -out allCOI_sulawesi_trimmed.aln -phylip -automated1
#export as fasta file as well; easier to get sequence lengths from here 
trimal -in allCOI_sulawesi.aln -out allCOI_sulawesi_trimmed.fasta -fasta -automated1

#copy file 
cp allCOI_sulawesi_trimmed.fasta allCOI_sulawesi_trimmed_noGap.fasta
#remove - symbol from the file 
sed -i 's/-//g' allCOI_sulawesi_trimmed_noGap.fasta
#get sequence lengths 
awk '/^>/ { if (seqlen) {
              print seqlen
              }
            print

            seqtotal+=seqlen
            seqlen=0
            seq+=1
            next
            }
    {
    seqlen += length($0)
    }     
    END{print seqlen
        print seq" sequences, total length " seqtotal+seqlen
    }' allCOI_sulawesi_trimmed_noGap.fasta > seqLength_noGap.txt
#remove > 
sed -i 's/>//g' seqLength_noGap.txt

#now make a tree! 
#load dependencies 
module load gcc/9.3.0; module load openmpi/4.0.3
#load raxml 
module load raxml-ng/1.0.1
#build
#replace piping symbol with an underscore so linux doesn't get upset 
sed -i 's/|/_/g' allCOI_sulawesi_trimmed.aln
#run raxml 
raxml-ng --all --msa allCOI_sulawesi_trimmed.aln --model GTR --prefix allCOI_sulawesi_trimmed.tree --seed 869276 --outgroup MCZ_IZ_89506_Guatemala --bs-metric tbe --tree rand{1} --bs-trees 100 

#download using sftp 
get allCOI* 
get seqLength_noGap.txt

