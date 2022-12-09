# Tyrone_XL.md

This is a repo for the XL project with Tyrone, Adonis, Cara, and Dan.

After trimming with trimmomatic, I mapped WGS data from 5 individuals to the XL v10 genome using bwa and samtools.  I called genotypes using GATK.  All methods were the same as detailed in 2020_GBS.

# genome file
```
/home/ben/projects/rrg-ben/ben/2022_Tyrone/XL_v10_concatscaf/XL_v10.1_concatenatedscaffolds.fa
```

# Individual polymorphism

I'm using ANGSD to calculate individual polymorphism.  I began by making a bamfilelist file; each line is the path to a bam file:
```
3396/3396_DNA_S122_L004_trim_sorted.bam_rg.bam
3399/3399_DNA_S123_L004_trim_sorted.bam_rg.bam
3732/3732_DNA_S124_L004_trim_sorted.bam_rg.bam
3733/3733_DNA_S125_L004_trim_sorted.bam_rg.bam
LG_3_1/LG_3_1_DNA_S126_L004_trim_sorted.bam_rg.bam
```
The first 4 individuals are two male and female ZZ individuals that produce ~50% sex reversed offspring. Then I'm using angsd to calculate polymorphism stats as detailed here (http://popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests#Full_command_list_for_below_examples):

```
module load angsd/0.936
./angsd -bam bamfilelist -doSaf 1 -anc ../../XL_v10_concatscaf/XL_v10.1_concatenatedscaffolds.fa -GL 1 -out out
realSFS out.saf.idx -P 24 -fold 1 > out.sfs
thetaStat do_stat out.thetas.idx
thetaStat do_stat out.thetas.idx -win 50000 -step 10000  -outnames theta.thetasWindow.gz
```
