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
```

