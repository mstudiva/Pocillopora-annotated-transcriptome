# Pocillopora spp. Transcriptome Annotation, version December 18, 2024
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (studivanms@gmail.com) for use on FAU's HPC (KoKo)


#------------------------------
# BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster
# terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions -
# please make sure to read them before copy-pasting.


#------------------------------
# setup

# To install Bioperl as a conda environment
conda create -y -n bioperl perl-bioperl

# getting scripts
cd ~/bin
git clone https://github.com/z0on/annotatingTranscriptomes.git
mv annotatingTranscriptomes/* .
rm -rf annotatingTranscriptomes
rm launcher_creator.py

git clone https://github.com/z0on/emapper_to_GOMWU_KOGMWU.git
mv emapper_to_GOMWU_KOGMWU/* .
rm -rf emapper_to_GOMWU_KOGMWU

git clone https://github.com/mstudiva/Pocillopora-annotated-transcriptome.git
mv Pocillopora-annotated-transcriptome/* .
rm -rf Pocillopora-annotated-transcriptome

# creating backup directory
mkdir backup

# creating annotation directory
cd
mkdir annotate
cd annotate


#------------------------------
# getting transcriptomes

# Buitrago-Lopez (2020) https://doi.org/10.1093/gbe/evaa184
# cds and protein translations from genome at https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036669915.1/

# Renaming gene identifiers for ease
sed -i 's/lcl|NC_/Pocillopora/' Pocillopora.fasta
sed -i 's/LOC/Pocillopora/' Pocillopora.fasta

# transcriptome statistics
conda activate bioperl
echo "seq_stats.pl Pocillopora.fasta > seqstats_Pocillopora.txt" > seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats.slurm

Pocillopora.fasta
-------------------------
37406 sequences.
1859 average length.
63966 maximum length.
204 minimum length.
N50 = 2481
69.5 Mb altogether (69529273 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
# Extracting contig, isogroup, and protein IDs into lookup tables

grep ">" Pocillopora.fasta | awk '{
    header = substr($1, 2);             # Remove the leading ">" from the first field
    match($0, /gene=([^]]+)/, gene);    # Extract the gene name
    print header "\t" gene[1];          # Print the full header and gene name
}' > Pocillopora_seq2iso.tab

grep ">" Pocillopora.fasta | awk -F'[][]' '{
    for (i=1; i<=NF; i++) {
        if ($i ~ /gene=/) { gsub("gene=", "", $i); gene=$i }
        if ($i ~ /protein_id=/) { gsub("protein_id=", "", $i); protein=$i }
    }
    if (gene && protein) { print protein, gene }
}' > Pocillopora_seq2pro.tab


#------------------------------
# GO annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

# first rename protein headers as gene ID based on lookup table
awk 'BEGIN {
    while (getline < "Pocillopora_seq2pro.tab") {
        map[$1] = $2
    }
}
/^>/ {
    protein_id = substr($0, 2, index($0, " ") - 2)
    if (protein_id in map) {
        sub(protein_id, map[protein_id])
    }
}
{ print }' protein.faa > Pocillopora_pro.fasta

# selecting the longest contig per isogroup
fasta2SBH.pl Pocillopora_pro.fasta >Pocillopora_pro_out.fasta

# scp your *_pro_out.fas file to laptop, submit it to
http://eggnog-mapper.embl.de
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/\*_pro_out.fasta .

# copy link to job ID status and output file, paste it below instead of current link:
# http://eggnog-mapper.embl.de/job_status?jobname=MM_nw1y1jd7

# once it is done, download results to HPC:
wget http://eggnog-mapper.embl.de/MM_nw1y1jd7/out.emapper.annotations

# GO:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$10 }' out.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Pocillopora_iso2go.tab
# gene names:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$8 }' out.emapper.annotations | grep -Ev "\tNA" >Pocillopora_iso2geneName.tab


#------------------------------
# KOG annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

cp ~/bin/kog_classes.txt .

#  KOG classes (single-letter):
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$7 }' out.emapper.annotations | grep -Ev "\tNA" >Pocillopora_iso2kogClass1.tab
# converting single-letter KOG classes to text understood by KOGMWU package (must have kog_classes.txt file in the same dir):
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Pocillopora_iso2kogClass1.tab > Pocillopora_iso2kogClass.tab


#------------------------------
# KEGG annotations:

# first rename fasta headers as gene ID (isogroup) rather than contig ID
awk 'BEGIN {
    while (getline < "Pocillopora_seq2iso.tab") {
        map[$1] = $2
    }
}
/^>/ {
    gene_id = substr($0, 2, index($0, " ") - 2)
    if (gene_id in map) {
        sub(gene_id, map[gene_id])
    }
}
{ print }' Pocillopora.fasta > Pocillopora_iso.fasta

# selecting the longest contig per isogroup:
fasta2SBH.pl Pocillopora_iso.fasta >Pocillopora_iso_out.fasta

# Sanity check: How many unique isogroups do we have?
grep ">" Pocillopora.fasta | sort | uniq | wc -l            # 37406
grep ">" Pocillopora_pro.fasta | sort | uniq | wc -l        # 30050
grep ">" Pocillopora_pro_out.fasta | sort | uniq | wc -l    # 25866
grep ">" Pocillopora_iso.fasta | sort | uniq | wc -l        # 37406
grep ">" Pocillopora_iso_out.fasta | sort | uniq | wc -l    # 25866
# The two _out files should roughly match

# scp *_iso_out.fasta to your laptop
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/\*_iso_out.fasta .
# use web browser to submit 4kegg.fasta file to KEGG's KAAS server (http://www.genome.jp/kegg/kaas/)
# select SBH method, upload nucleotide query
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1734547392&key=541MLaVJ

# Once it is done, download to HPC - it is named query.ko by default
wget https://www.genome.jp/tools/kaas/files/dl/1734547392/query.ko

# selecting only the lines with non-missing annotation:
cat query.ko | awk '{if ($2!="") print }' > Pocillopora_iso2kegg.tab

# the KEGG mapping result can be explored for completeness of transcriptome in terms of genes found,
# use 'html' output link from KAAS result page, see how many proteins you have for conserved complexes and pathways, such as ribosome, spliceosome, proteasome etc


#------------------------------
# file transfer

# copy all files to local machine
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/\* .
