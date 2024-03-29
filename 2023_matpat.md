# MatPat analysis

The combined files for each family are now here:
```
/home/ben/projects/rrg-ben/ben/2022_Tyrone/RADseq/RADseq_and_WGS_combined/family1
/home/ben/projects/rrg-ben/ben/2022_Tyrone/RADseq/RADseq_and_WGS_combined/family2
```

I'd like to explore the association signals of maternal and paternal sites separately. To do this I need to combine the WGS data (from the parents) with the RADseq data (from the offspring).

After running haplotype caller, I used GenomicsDB to merge the files, and then GenotypeGVCFs to genotype them. Then I used SelectVariants and VariantFiltration to filter them

Then I used vcftools to filter sites based on missingness like this:

```
vcftools --gzvcf DB__Chr9_10L_out.vcf_filtered.vcf.gz_filtered_removed.vcf.gz --max-missing 0.9 --recode --out Tyrone_XL_Chr9_10L_maxmissing_0.9
```

This gave me much smaller files that have only high quality sites that are in the parents and most offspring.

Then I separated out each family like this:
```
bcftools view -S bam.filelist_family2 Tyrone_XL_Chr1L_maxmissing_0.9.recode.vcf > Tyrone_XL_Chr1L_maxmissing_0.9_family2.vcf
```
the first family has these individuals (B10_TAGCGGA_sample_153_M_cuttrim_sorted.bam is a typo and should be B10_TAGCGGA_sample_163_M_cuttrim_sorted.bam):
```
3732_DNA_S124_L004_trim_sorted.bam
3733_DNA_S125_L004_trim_sorted.bam
A1_CTCG_sample_1_F_cuttrim_sorted.bam
A2_AGCG_sample_22_M_cuttrim_sorted.bam
A3_TTCTG_sample_33_F_cuttrim_sorted.bam
A4_ATTGA_sample_46_M_cuttrim_sorted.bam
A5_TCGTT_sample_57_F_cuttrim_sorted.bam
A8_CTTGCTT_sample_129_M_cuttrim_sorted.bam
B1_TGCA_sample_2_F_cuttrim_sorted.bam
B2_GATG_sample_23_M_cuttrim_sorted.bam
B3_AGCCG_sample_36_MQ_cuttrim_sorted.bam
B4_CATCT_sample_48_M_cuttrim_sorted.bam
B5_GGTTGT_sample_58_F_cuttrim_sorted.bam
C10_TCGAAGA_sample_173_M_cuttrim_sorted.bam
C1_ACTA_sample_4_M_cuttrim_sorted.bam
C2_TCAG_sample_24_M_cuttrim_sorted.bam
C3_GTATT_sample_37_F_cuttrim_sorted.bam
C4_CCTAG_sample_49_M_cuttrim_sorted.bam
C5_CCACGT_sample_60_F_cuttrim_sorted.bam
D1_CAGA_sample_5_M_cuttrim_sorted.bam
D2_TGCGA_sample_25_F_cuttrim_sorted.bam
D3_CTGTA_sample_38_F_cuttrim_sorted.bam
D4_GAGGA_sample_50_F_cuttrim_sorted.bam
D5_TTCAGA_sample_61_F_cuttrim_sorted.bam
E1_AACT_sample_6_M_cuttrim_sorted.bam
E2_CGCTT_sample_26_F_cuttrim_sorted.bam
E3_ACCGT_sample_40_F_cuttrim_sorted.bam
E4_GGAAG_sample_52_F_cuttrim_sorted.bam
E5_TAGGAA_sample_62_MQ_cuttrim_sorted.bam
E7_GTTGAA_sample_102_M_cuttrim_sorted.bam
F1_GCGT_sample_8_M_cuttrim_sorted.bam
F2_TCACG_sample_28_F_cuttrim_sorted.bam
F3_GCTTA_sample_41_F_cuttrim_sorted.bam
F4_GTCAA_sample_53_F_cuttrim_sorted.bam
F5_GCTCTA_sample_64_F_cuttrim_sorted.bam
F7_TAACGA_sample_105_M_cuttrim_sorted.bam
G1_CGAT_sample_9_F_cuttrim_sorted.bam
G2_CTAGG_sample_29_M_cuttrim_sorted.bam
G3_GGTGT_sample_42_F_cuttrim_sorted.bam
G4_TAATA_sample_54_F_cuttrim_sorted.bam
G5_CCACAA_sample_65_M_cuttrim_sorted.bam
G7_TGGCTA_sample_122_M_cuttrim_sorted.bam
H10_TAGGCCAT_sample_193_M_cuttrim_sorted.bam
H1_GTAA_sample_10_M_cuttrim_sorted.bam
H2_ACAAA_sample_31_M_cuttrim_sorted.bam
H3_AGGAT_sample_43_F_cuttrim_sorted.bam
H4_TACAT_sample_56_F_cuttrim_sorted.bam
H5_CTTCCA_sample_66_F_cuttrim_sorted.bam
H7_TATTTTT_sample_125_MQ_cuttrim_sorted.bam
B10_TAGCGGA_sample_153_M_cuttrim_sorted.bam
```
the second family has these individuals:
```
3396_DNA_S122_L004_trim_sorted.bam
3399_DNA_S123_L004_trim_sorted.bam
A6_GAGATA_sample_84_F_cuttrim_sorted.bam
B6_ATGCCT_sample_87_M_cuttrim_sorted.bam
C6_AGTGGA_sample_89_F_cuttrim_sorted.bam
D6_ACCTAA_sample_91_F_cuttrim_sorted.bam
E6_ATATGT_sample_93_M_cuttrim_sorted.bam
F6_ATCGTA_sample_94_F_cuttrim_sorted.bam
G6_CATCGT_sample_95_F_cuttrim_sorted.bam
H6_CGCGGT_sample_96_M_cuttrim_sorted.bam
A7_CTATTA_sample_97_F_cuttrim_sorted.bam
B7_GCCAGT_sample_99_F_cuttrim_sorted.bam
C7_GGAAGA_sample_100_F_cuttrim_sorted.bam
D7_GTACTT_sample_101_F_cuttrim_sorted.bam
B8_ATGAAAG_sample_131_F_cuttrim_sorted.bam
C8_AAAAGTT_sample_132_F_cuttrim_sorted.bam
D8_GAATTCA_sample_133_F_cuttrim_sorted.bam
E8_GAACTTG_sample_134_M_cuttrim_sorted.bam
F8_GGACCTA_sample_142_M_cuttrim_sorted.bam
G8_GTCGATT_sample_143_F_cuttrim_sorted.bam
H8_AACGCCT_sample_144_F_cuttrim_sorted.bam
A9_AATATGG_sample_145_M_cuttrim_sorted.bam
B9_ACGTGTT_sample_146_F_cuttrim_sorted.bam
C9_ATTAATT_sample_147_F_cuttrim_sorted.bam
D9_ATTGGAT_sample_148_F_cuttrim_sorted.bam
E9_CATAAGT_sample_149_F_cuttrim_sorted.bam
F9_CGCTGAT_sample_150_F_cuttrim_sorted.bam
G9_CGGTAGA_sample_151_F_cuttrim_sorted.bam
H9_CTACGGA_sample_152_F_cuttrim_sorted.bam
A10_GCGGAAT_sample_153_M_cuttrim_sorted.bam
D10_TCTGTGA_sample_178_F_cuttrim_sorted.bam
E10_TGCTGGA_sample_179_F_cuttrim_sorted.bam
F10_ACGACTAG_sample_180_M_cuttrim_sorted.bam
G10_TAGCATGG_sample_181_F_cuttrim_sorted.bam
A11_TGCAAGGA_sample_194_F_cuttrim_sorted.bam
B11_TGGTACGT_sample_195_F_cuttrim_sorted.bam
C11_TCTCAGTG_sample_196_F_cuttrim_sorted.bam
D11_CGCGATAT_sample_197_F_cuttrim_sorted.bam
E11_CGCCTTAT_sample_199_M_cuttrim_sorted.bam
F11_AACCGAGA_sample_202_M_cuttrim_sorted.bam
G11_ACAGGGA_sample_208_F_cuttrim_sorted.bam
H11_ACGTGGTA_sample_209_F_cuttrim_sorted.bam
A12_CCATGGGT_sample_210_F_cuttrim_sorted.bam
B12_CGCGGAGA_sample_214_M_cuttrim_sorted.bam
C12_CGTGTGGT_sample_217_F_cuttrim_sorted.bam
D12_GCTGTGGA_sample_218_F_cuttrim_sorted.bam
E12_GGATTGGT_sample_221_M_cuttrim_sorted.bam
F12_GTGAGGGT_sample_222_M_cuttrim_sorted.bam
G12_TATCGGGA_sample_223_M_cuttrim_sorted.bam
H12_TTCCTGGA_sample_230_M_cuttrim_sorted.bam
```

Then I used this perl script to make lists of sites with parent-specific heterozygosity:
```
#!/usr/bin/env perl
use strict;
use warnings;

# This program reads in a vcf file with genotypic information from
# a family  and identifies positions that
# are heterozygous in only the mother or only the father.
# module load StdEnv/2023 perl/5.36.1
# execute like this:
# ./Gets_matpat_positions_from_vcf_file.pl vcf mat pat matout patout

# where mat is the sample position of the mother 
# where pat is the sample position of the father
# and matout and patout are the output files (with chromsoome name in them) 

my $vcf = $ARGV[0];
my $mat = $ARGV[1];
my $pat = $ARGV[2];
my $outputfile = $ARGV[3];
my $outputfile2 = $ARGV[4];
my @columns;
my @mat;
my @pat;
my @mat1;
my @pat1;


unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile   $!\n\n";
	exit;
}
print "Creating output file: $outputfile\n";

unless (open(OUTFILE2, ">$outputfile2"))  {
	print "I can\'t write to $outputfile2  $!\n\n";
	exit;
}
print "Creating output file: $outputfile2\n";



#### Prepare the input file  with values for svl

#unless (open DATAINPUT, $vcf) {
#	print "Can not find the vcf file, jackass.\n";
#	exit;
#}


if ($vcf =~ /.gz$/) {
	#open DATAINPUT, '<:gzip', $vcf or die "Could not read from $vcf: $!";
	open(DATAINPUT, "gunzip -c $vcf |") || die "can’t open pipe to $vcf";
}
else {
	open DATAINPUT, $vcf or die "Could not read from $vcf: $!";
}


while ( my $line = <DATAINPUT>) {
	@columns=split("	",$line);
		#print $line,"\n";
		if($columns[0] =~ m/^#/){
			if($columns[0] eq '#CHROM'){
				print $columns[$mat+8],"\t",$columns[$pat+8],"\n";
			}	
		}
		else{
			@mat = split(":",$columns[$mat+8]);
			@pat = split(":",$columns[$pat+8]);
			# select positions that are not missing in either parent
			print $mat[0]," ",$pat[0],"\n";
			if(($mat[0] ne './.')&&($pat[0] ne './.')&&($mat[0] ne '.|.')&&($pat[0] ne '.|.')){
				@mat1 = split(/[\|\/]/,$mat[0]);
				@pat1 = split(/[\|\/]/,$pat[0]);
				if(($mat1[0] ne $mat1[1])&&($pat1[0] eq $pat1[1])){ # this is a mat site
					print OUTFILE $columns[0],"\t",$columns[1],"\n";
					print "mat\t",$columns[0],"\t",$columns[1],"\t",$mat1[0],"\t",$mat1[1],"\t",$pat1[0],"\t",$pat1[1],"\n";
				}
				elsif(($mat1[0] eq $mat1[1])&&($pat1[0] ne $pat1[1])){ # this is a pat site
					print OUTFILE2 $columns[0],"\t",$columns[1],"\n";
					print "pat\t", $columns[0],"\t",$columns[1],"\t",$mat1[0],"\t",$mat1[1],"\t",$pat1[0],"\t",$pat1[1],"\n";
				}
			}
		} # end else
} # end while
close DATAINPUT;
close OUTFILE;
close OUTFILE2;
```
