# Raw data
```
/home/ben/projects/rrg-ben/ben/2022_Tyrone/RADseq
```

# Demultiplex with sabre

```
#!/bin/sh
#SBATCH --job-name=sabre_plate1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=2gb
#SBATCH --output=sabre_plate1.%J.out
#SBATCH --error=sabre_plate1.%J.err
#SBATCH --account=def-ben


module load nixpkgs/16.09  intel/2016.4
module load module load sabre/1.00

sabre pe -f NS.2027.004.D701---B504.Hayes082308GBS01_R1.fastq.gz -r NS.2027.004.D701---B504.Hayes082308GBS01_R2.fastq.gz -b Tyrone_sabre_barcodes.txt -u
 no_bc_match_R1.fq -w no_bc_match_R2.fq
 ```
 
 # Cutadapt (removes adapter seqs)
 ```
 #!/bin/sh
#SBATCH --job-name=cutadapt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=2gb
#SBATCH --output=cutadapt.%J.out
#SBATCH --error=cutadapt.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2020_cutadapt.sh ../raw_data/plate1


module load python/3.7

# this allows for 20% mismatch of adapters, and trims from the 5' and 3'
# ends based on quality scores. the -B flag tells cutadapt to trim
# adaptors anywhere in the read.  The adapter sequences are the same for
# both directions because the end of the primer that would be
# incorporated into the seq in the 3' end of each read is identical.

#~/.local/bin/cutadapt -a "AGATCGGAAGAGC;max_error_rate=0.2" -A "AGATCGGAAGAGC;max_error_rate=0.2" -B -q 15,10 -o
 out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq


v=1
#  Always use for-loop, prefix glob, check if exists file.
for file in $1/*.fq ; do         # Use ./* ... NEVER bare *
  if [ -e "$file" ] ; then   # Check whether file exists.
  	if [[ $v -eq 1 ]]
	then # if/then branch
	    ~/.local/bin/cutadapt -b "AGATCGGAAGAGC" -B "AGATCGGAAGAGC" -e 0.2 -q 15,10 -o ${file::-5}cut.R1.fq -
p ${file::-5}cut.R2.fq ${file::-5}R1.fq ${file::-5}R2.fq
		  v=0
	else # else branch
  		v=1
	fi
  fi
done 
```
# Trimmomatic
```
#!/bin/sh
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=8gb
#SBATCH --output=trimmomatic.%J.out
#SBATCH --error=trimmomatic.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2020_trimmomatic.sh ../raw_data/plate1


module load StdEnv/2020
module load trimmomatic/0.39

#v=1
#  Always use for-loop, prefix glob, check if exists file.
for file in $1/*_cut.R1.fq ; do         # Use ./* ... NEVER bare *
  if [ -e "$file" ] ; then   # Check whether file exists.
  	#if [[ $v -eq 1 ]]
	#then # if/then branch
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE ${file::-10}_cut.R1.fq ${file::-10}_cut.R2.
fq ${file::-10}_cuttrim.R1.fq.gz ${file::-10}_cuttrim.R1_single.fq.gz ${file::-10}_cuttrim.R2.fq.gz ${fi
le::-10}_cuttrim.R2_single.fq.gz SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:3
	#	  v=0
	#else # else branch
  	#	v=1
	#fi
  fi
done 
```

# Cuttrim fq files
```
/home/ben/projects/rrg-ben/ben/2022_Tyrone/RADseq/raw_data/cuttrim
```

# Ref genome
```
/home/ben/projects/rrg-ben/ben/2022_Tyrone/XL_v10_concatscaf/XL_v10.1_concatenatedscaffolds.fa
```
# Align to ref with concatenated scaffolds
```
#!/bin/sh
#SBATCH --job-name=bwa_align
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --mem=32gb
#SBATCH --output=bwa_align.%J.out
#SBATCH --error=bwa_align.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2020_align_paired_fq_to_ref.sh pathandname_of_ref path_to_paired_fq_filez
# sbatch 2020_align_paired_fq_to_ref.sh /home/ben/projects/rrg-ben/ben/2018_Austin_XB_genome/Austin_geno
me/Xbo.v1.fa.gz pathtofqfilez
# or for XL genome use:
# sbatch 2020_align_paired_fq_to_ref.sh /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_
genome.fa.gz pathtofqfilez

module load bwa/0.7.17
module load samtools/1.10


for file in ${2}/*.R1.fq.gz ; do         # Use ./* ... NEVER bare *    
    if [ -e "$file" ] ; then   # Check whether file exists.
	#${file::-9}
	echo bwa mem ${1} ${file::-9}.R1.fq.gz ${file::-9}.R2.fq.gz -t 16 | samtools view -Shu - | samto
ols sort - -o ${file::-9}_sorted.bam
	bwa mem ${1} ${file::-9}.R1.fq.gz ${file::-9}.R2.fq.gz -t 16 | samtools view -Shu - | samtools s
ort - -o ${file::-9}_sorted.bam
	samtools index ${file::-9}_sorted.bam
  fi
done
```

# Add readgroups and index bam
```
#!/bin/sh
#SBATCH --job-name=readgroups
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=8gb
#SBATCH --output=readgroups.%J.out
#SBATCH --error=readgroups.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2021_picard_add_read_groups.sh /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_la
ev/raw_data/cutaddapted_by_species_across_three_plates/clivii/ 

module load picard/2.23.3

for file in ${1}*_sorted.bam
do
    java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${file} O=${file}_rg.bam RGID=4 RGLB=$(b
asename $file) RGPL=ILLUMINA RGPU=$(basename $file) RGSM=$(basename $file)
done

module load StdEnv/2020 samtools/1.12


for file in ${1}*_sorted.bam_rg.bam
do
    samtools index ${file}
done
```

# Call genotypes on each bam file
```
#!/bin/sh
#SBATCH --job-name=HaplotypeCaller
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
#SBATCH --mem=30gb
#SBATCH --output=HaplotypeCaller.%J.out
#SBATCH --error=HaplotypeCaller.%J.err
#SBATCH --account=def-ben


# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute the GATK command "RealignerTargetCreator" on these files. 

# execute like this:
# sbatch 2021_HaplotypeCaller.sh /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa /home/ben/projec
ts/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/clivii/ 

module load nixpkgs/16.09 gatk/4.1.0.0

for file in ${2}*_sorted.bam_rg.bam
do
    gatk --java-options -Xmx24G HaplotypeCaller  -I ${file} -R ${1} -O ${file}.g.vcf -ERC GVCF
done
```

# Use Genomics DB to combine vcf files
```
#!/bin/sh
#SBATCH --job-name=GenomicsDBImport
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=60:00:00
#SBATCH --mem=24gb
#SBATCH --output=GenomicsDBImport.%J.out
#SBATCH --error=GenomicsDBImport.%J.err
#SBATCH --account=def-ben


# This script will read in the *.g.vcf file names in a directory, and 
# make and execute the GATK command "GenotypeGVCFs" on these files. 

# execute like this:
# sbatch 2021_GenomicsDBImport.sh /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa /home/ben/proje
cts/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/clivii/vcf/ chr1L tem
p_path_and_dir db_path_and_dir_prefix

module load nixpkgs/16.09 gatk/4.1.0.0

commandline="gatk --java-options -Xmx20G GenomicsDBImport -R ${1}"
for file in ${2}*g.vcf
do
    commandline+=" -V ${file}"
done

commandline+=" -L ${3} --tmp-dir=${4} --batch-size 50 --genomicsdb-workspace-path ${5}_${3}"

${commandline}
```

# GenotypeGVCFs
```
#!/bin/sh
#SBATCH --job-name=GenotypeGVCFs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=23gb
#SBATCH --output=GenotypeGVCFs.%J.out
#SBATCH --error=GenotypeGVCFs.%J.err
#SBATCH --account=def-ben


# This script will execute the GATK command "GenotypeGVCFs" on
# folders generated by GenomicDBImport

# execute like this:
# sbatch 2021_GenotypeGVCFs_DB.sh ref DB_path_and_prefix chr1L

module load nixpkgs/16.09 gatk/4.1.0.0

commandline="gatk --java-options -Xmx18G GenotypeGVCFs -R ${1} -V gendb://${2}_${3} -O ${2}_${3}_out.vcf"

${commandline}
```

# Test for phenotypic associations with Plink

## Concatenate vcf files
```
module load StdEnv/2020  gcc/9.3.0 bcftools/1.13

bcftools concat DB_Chr1L_Chr1L_out.vcf DB_Chr1S_Chr1S_out.vcf DB_Chr2L_Chr2L_out.vcf DB_Chr2S_Chr2S_out.vcf DB_Chr3L_Chr3L_out.vcf DB_Chr3S_Chr3S_out.vcf DB_Chr4L_Chr4L_out.vcf DB_Chr4S_Chr4S_out.vcf DB_Chr5L_Chr5L_out.vcf DB_Chr5S_Chr5S_out.vcf DB_Chr6L_Chr6L_out.vcf DB_Chr6S_Chr6S_out.vcf DB_Chr7L_Chr7L_out.vcf DB_Chr7S_Chr7S_out.vcf DB_Chr8L_Chr8L_out.vcf DB_Chr8S_Chr8S_out.vcf DB_Chr9_10L_Chr9_10L_out.vcf DB_Chr9_10S_Chr9_10S_out.vcf DB_Scaffolds_Scaffolds_out.vcf -O z -o TyroneRADseq_unfiltered_allChrs.vcf.gz
```
## run plink
```
module load nixpkgs/16.09 plink/1.9b_5.2-x86_64

plink --vcf TyroneRADseq_unfiltered_allChrs.vcf.gz --recode --const-fid 0 --chr-set 37 no-y no-xy no-mt --allow-extra-chr --out TyroneRADseq_unfiltered_allChrs.vcf.gz_myplink

plink --file TyroneRADseq_unfiltered_allChrs.vcf.gz_myplink --pheno sex_phenotype --assoc --allow-no-sex --allow-extra-chr

mv plink.assoc TyroneRADseq_unfiltered_allChrs.vcf.gz_myplink.assoc
```

# Plot and interpret results
```R
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Tyrone/plink")
library(ggplot2)
dat <-read.table("TyroneRADseq_unfiltered_allChrs.vcf.gz_myplink.assoc", header=TRUE)

#newdat <- dat[!is.na(dat$P),]
pdf("./unfiltered_mappedtoXL_plink_assoc_long.pdf",w=8, h=24.0, version="1.4", bg="transparent")
p<-ggplot(dat, aes(x=BP, y=-log10(P))) + 
    # add points
    geom_point(size=2, alpha = 0.7 ) +
    # color the stuff the way I want
    facet_wrap(~CHR, ncol = 1) +
    # get rid of gray background
    theme_bw()
p
dev.off()



# get top ten max -logP values
# make a -logP column
dat$minuslogP <- -log10(dat$P)
# reverse sort the df by -logP
head(dat[order(-dat$minuslogP),],n=20)

# get a specific region
subset<-dat[(dat$CHR == 'Chr5L'),]

# get a specific region (for allofraseri)
#subset<-dat[(dat$CHR == 'chr7L')&(dat$BP >= 14368095)&(dat$BP <= 33688731),]

# plot it
p<-ggplot(subset, aes(x=BP, y=-log10(P))) + 
    # add points
    geom_point(size=2, alpha = 0.7 ) +
    # get rid of gray background
    theme_bw()
p

pdf("./Chr5L_closeup.pdf",w=8, h=4.0, version="1.4", bg="transparent")
    p
dev.off()    



################################################
#
#
#  Deterine which RADseq SNPs are in exons
#
#
################################################


#BiocManager::install("genomation")
#library(genomation)
library("GenomicRanges")
library(genomation)
#setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/plink_somemappedtoXLv9.2_othersmappedtoAustinXB")
# read in the data from plink output
#dat <-read.table("pygmaeus_filtered_removed_allchrs.vcf.gz_plink_noNAs.assoc", header=TRUE)
# make a -logP column
dat$minuslogP <- -log10(dat$P)
# reverse sort the df by -logP
head(dat[order(-dat$minuslogP),],n=20)


gff <- gffToGRanges("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/plink_somemappedtoXLv9.2_othersmappedtoAustinXB/XENLA_9.2_Xenbase.gtf", filter = NULL, zero.based = FALSE, ensembl = FALSE)

# Get rid of rows with NA for P value
dat_NoNAs <- dat[!is.na(dat$minuslogP),]

# get the RADseq SNPs that have a strong bias
subset<-dat_NoNAs[(dat_NoNAs$CHR == 'Chr5L')&(dat_NoNAs$minuslogP >= 6),];subset
dim(subset)
# subset gtf and keep only the ranges exons in chr8L that are exons
Chr5L_exons <- gff %>% subset(seqnames == 'Chr5L') %>% subset(type == 'exon')
# make a GenomicRanges object from the RADseq SNPs
SL_sites <- makeGRangesFromDataFrame(subset,
                                     keep.extra.columns=FALSE,
                                     ignore.strand=FALSE,
                                     seqinfo=NULL,
                                     seqnames.field=c("seqnames", "seqname",
                                                      "chromosome", "chrom",
                                                      "chr", "chromosome_name",
                                                      "seqid"),
                                     start.field="BP",
                                     end.field=c("BP"),
                                     strand.field="strand",
                                     starts.in.df.are.0based=FALSE)

o = findOverlaps(SL_sites,Chr5L_exons);o
SL_sites = split(SL_sites[queryHits(o)], 1:length(o)) # You can't mendoapply on a GRanges object
Chr5L_exons = split(chr8L_exons[subjectHits(o)], 1:length(o))
foo = function(x, y) {
    rv = x
    start(rv) = max(start(x), start(y))
    end(rv) = min(end(x), end(y))
    return(rv)
}
unlist(mendoapply(foo, SL_sites, y=chr8L_exons))





















################################################
#
#
#  Deterine which RADseq SNPs are in exons
#
#
################################################


# reverse sort the df by -logP
head(dat[order(-dat$minuslogP),],n=20)


library("GenomicRanges")
library(genomation)
#makeGRangesFromGFF(
#    "XENLA_9.2_Xenbase.gtf",
#    level = c("genes", "transcripts"),
#    ignoreVersion = TRUE,
#    synonyms = FALSE
#)
#BiocManager::install("genomation")
#library(genomation)

gff <- gffToGRanges("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/plink_somemappedtoXLv9.2_othersmappedtoAustinXB/XENLA_9.2_Xenbase.gtf", filter = NULL, zero.based = FALSE, ensembl = FALSE)


# get the RADseq SNPs that have a strong bias
subset<-dat[(dat$CHR == 'Chr5L')&(dat$minuslogP >= 6),];subset
    
dim(subset)
# subset gtf and keep only the ranges exons in chr8L that are exons
chr8L_exons <- gff %>% subset(seqnames == 'chr8L') %>% subset(type == 'exon')
# make a GenomicRanges object from the RADseq SNPs
SL_sites <- makeGRangesFromDataFrame(subset,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="BP",
                         end.field=c("BP"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)

o = findOverlaps(SL_sites,chr8L_exons);o
SL_sites = split(SL_sites[queryHits(o)], 1:length(o)) # You can't mendoapply on a GRanges object
chr8L_exons = split(chr8L_exons[subjectHits(o)], 1:length(o))
foo = function(x, y) {
    rv = x
    start(rv) = max(start(x), start(y))
    end(rv) = min(end(x), end(y))
    return(rv)
}
unlist(mendoapply(foo, SL_sites, y=chr8L_exons))






#grhits=intersect(SL_sites,chr8L_exons);grhits

# test df
chr <- rep("chr8L",2);chr
begin <- c(1,49869)
end <- c(1,49869)
test_df<- data.frame(chr, begin, end)
SL_sites <- makeGRangesFromDataFrame(test_df,
                                     keep.extra.columns=FALSE,
                                     ignore.strand=FALSE,
                                     seqinfo=NULL,
                                     seqnames.field=c("seqnames", "seqname",
                                                      "chromosome", "chrom",
                                                      "chr", "chromosome_name",
                                                      "seqid"),
                                     start.field="begin",
                                     end.field=c("end"),
                                     strand.field="strand",
                                     starts.in.df.are.0based=FALSE)

#grhits=intersect(SL_sites,chr8L_exons);grhits
o = findOverlaps(SL_sites,chr8L_exons);o
SL_sites = split(SL_sites[queryHits(o)], 1:length(o)) # You can't mendoapply on a GRanges object
chr8L_exons = split(chr8L_exons[subjectHits(o)], 1:length(o))
foo = function(x, y) {
         rv = x
         start(rv) = max(start(x), start(y))
         end(rv) = min(end(x), end(y))
         return(rv)
     }
unlist(mendoapply(foo, SL_sites, y=chr8L_exons))

```

# Check coverage of high association positions
```
[ben@gra-login3 combined_vcfs]$ zcat TyroneRADseq_unfiltered_allChrs.vcf.gz | grep 'Chr5L' | grep '77572949' | grep -o '\./\.' | wc -l
53
[ben@gra-login3 combined_vcfs]$ zcat TyroneRADseq_unfiltered_allChrs.vcf.gz | grep 'Chr5L' | grep '98706723' | grep -o '\./\.' | wc -l
115 (this is not correct)
[ben@gra-login3 combined_vcfs]$ zcat TyroneRADseq_unfiltered_allChrs.vcf.gz | grep 'Chr5L' | grep '76155133' | grep -o '\./\.' | wc -l
23
[ben@gra-login3 combined_vcfs]$ zcat TyroneRADseq_unfiltered_allChrs.vcf.gz | grep 'Chr5L' | grep '76155153' | grep -o '\./\.' | wc -l
23
```


