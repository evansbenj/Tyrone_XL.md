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


