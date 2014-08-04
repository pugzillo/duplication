#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper; 
use Cwd;
use threads;

###parse_MCL_outfile_lineageSpecificDups.pl info:
#Perl script to find lineage specific duplicates in MCL clustering output, align them with muscle, create species guided phylogenies with TreeBest, and parse each tree to find orthologs and paralogs for each gene family. Also, for each gene family, perform paml to look at evolutionary rates  
#
##Programs which must be loaded:
#Perl
#Python
#EMBOSS
#PHYML
#paml
#
#In the directory, there must be:
#cds files for the species that MCL was performed on (multifasta file should be linearized)
#the files name should be in the format of speciesname_version_cds.fa

################################ VARIABLES ################################ 
###Subroutines
sub codeml;
sub UsageandExit; 

###Input variables
if ($#ARGV != 12 ) {;}
my $MCLout = $ARGV[0]; 			##mcl output from blast all v all and then single linkage clustering with MCL
my $genePair = $ARGV[1];		##Prefix for gene family
my $numProc = $ARGV[2];			##Number of processors
my $familydir = $ARGV[3];		##Name of directory for all gene family files
my $SC_file = $ARGV[4]; 		###File for single copy genes
my $speciesTree = $ARGV[5];		###File with Species tree for TreeBeST
my $species1 = $ARGV[6];		###Species 1, must match with the species tree
my $species2 = $ARGV[7];		###Species 2, must match with the species tree
my $species3 = $ARGV[8];		###Species 3, must match with the species tree
my $sp1 = $ARGV[9]; 			###Sequence prefix for species 1(Make sure sequence ids are less that 10 char; will be altered by phyml)
my $sp2 = $ARGV[10]; 			###Sequence prefix for species 2
my $sp3 = $ARGV[11]; 			###Sequence prefix for species 3
my $foreground = $ARGV[12]; 	###PAML foreground branch; Input the sequence prefix for the species used as the foreground,branch site test
$foreground =~ s/[0-9]//g;

###Variables
my $cwd=getcwd();
my $count = 1; 

###Data Structures
my (%family, %family2, %sequence, @Familylist, @Familylist2, %single_copy, @singleCopy)=();

###Directories
mkdir $familydir;

my $paml_dir = $cwd . "/paml/"; 
mkdir $paml_dir;

my $paml_nulldir = $cwd . "/paml/null"; 
mkdir $paml_nulldir;

my $paml_freeratio_dir = $cwd . "/paml/freeRatio"; 
mkdir $paml_freeratio_dir;

my $paml_BSnull_dir = $cwd . "/paml/BSnull"; 
mkdir $paml_BSnull_dir;

my $paml_BSpos_dir = $cwd . "/paml/BSpos"; 
mkdir $paml_BSpos_dir;

###Populate Hashes with Sequence Data for each gene set; make sure that there is no sequence wrapping!!!!
my @cds_files=glob("*.fa"); #pulls all fasta files into an array (READ IN CDS!!!!!!!)

######## READING IN SEQUENCE FILES AND MCL CLUSTERING OUTPUT FILES ##############
for my $fasta(@cds_files){	####Reads through each fasta file and creates hash with gene id and the sequence
	open (IN, $fasta) || die "$fasta file is not found!\n";
	print "Reading in $fasta\n";
	$/ = ">"; 
	my $junk = <IN>; #Discared the ">" at the beginning of the file

	while (my $record = <IN>){ #Read in multifasta file
		chomp $record; 			##Remove the ">" from the end of $record 
		my ($seq_id, @seqLines) = split("\n", $record);		###put sequence ID and sequence into variables
		my $seq = join("", @seqLines);  	
		my $seqLength = length($seq); 		###length of sequence
		#If there are multiple sequences for a gene, checks if it's the longest and uses that for the sequence hash
		if (exists $sequence{$seq_id}){
			if ($seqLength > length($sequence{$seq_id})){
				$sequence{$seq_id} = $seq;
			}else{
				next; 
			}
		}else{
			$sequence{$seq_id} = $seq; ###Hash with the id and sequence for each gene in the genome file
		}
	}
	close(IN);
}

###Reading in MCL output file
open (INPUT, $MCLout) || die "$MCLout is not found!\n";
print "Reading in $MCLout\n";
$/ = "\n";
foreach (<INPUT>){ #Read in the MCL output lines (1 line= 1gene family) 
	chomp $_;	
	my @geneFam = split("\t", $_);  ###array containing all the genes in one family	
	foreach(@geneFam){				###Remove extra spaces in the file; can fuck up parsing later on
		$_=~s/\s//g;
	}
	my $number = scalar(@geneFam); 	###Number of genes in that gene family
	
	if ($number >= 2){				###For gene families size = 2 or more
		my $FamID = $genePair."_".$count; ###gene family ID
		for my $geneID(@geneFam){
			if (exists $sequence{$geneID}){		###Is this gene from clsutering in the sequence hash?
				my $geneID2; 
				if ($geneID =~ $sp1){
					$geneID2 = $geneID."_".$species1; 
				}
				if ($geneID =~ $sp2){
					$geneID2 = $geneID."_".$species2; 
				}
				if ($geneID =~ $sp3){
					$geneID2 = $geneID."_".$species3; 
				}
				my $line = $geneID."\t".$sequence{$geneID}; ####Hash with gene family ID as the key and the geneIDs and sequences as the values delimited by a new line
				my $line2 = $geneID2."\t".$sequence{$geneID}; ####Hash with gene family ID as the key and the geneIDs (with species) and sequences as the values delimited by a new line
				push @{$family{$FamID}}, $line; 
				push @{$family2{$FamID}}, $line2;				
			}
		}
		$count++; 
	}else{
		for my $singleGene(@geneFam){ ###An array for single copy genes (gene family n==1) along with sequence
			if (exists $sequence{$singleGene}){
				my $singleLine = ">".$singleGene."\n".$sequence{$singleGene}; 
				push (@singleCopy, $singleLine); 
			}
		}
	}
}close(INPUT);
 

####Print out a multi-fasta file for each gene Family 
for my $fam (keys %family){
	open (OUT, ">$familydir"."/"."$fam.fa") || die "$fam outfile cannot be opened!\n";
	push @Familylist,$fam; #Array will the list of gene families
	foreach(@{$family{$fam}}){
		my ($id, $sequence) = split("\t", $_);
		print OUT ">$id\n$sequence\n";	
	}
}close(OUT);

####Print out a multi-fasta file for each gene Family 
for my $fam2 (keys %family2){
	open (OUT, ">$familydir"."/"."$fam2.sp.fa") || die "$fam2 outfile cannot be opened!\n";
	push (@Familylist2,$fam2); #Array will the list of gene family files
	foreach(@{$family2{$fam2}}){
		my ($id, $sequence) = split("\t", $_);
		print OUT ">$id\n$sequence\n";	
	}
}close(OUT);

print Dumper\@Familylist2; 

####Print out file with single copy genes
open (OUTPUT, ">$SC_file") || die $!;
for my $singleCopyGene (@singleCopy){
	print OUTPUT $singleCopyGene ."\n";
	}
close (OUTPUT);

####Change to orthofile directory
chdir ($familydir); 

################### MSA, PHYLOGENY, and Homology FOR EACH GENE FAMILY ###################
my $switch1 = 0; 

for my $family (@Familylist2){
	if($switch1 ==1){
		while(1){ ### loop that runs infinitely until one of the running threads completes.  waits N seconds between each loop
			my @running = threads->list(threads::running);
			my @joinable = threads->list(threads::joinable);
			$_->join() for(@joinable);
			last if(@joinable > 0 or @running < $numProc);  ###exit the stationary loop if any threads are joined (end) in the last pass
			sleep 3;
		}
	}
	async{
		unless(-e $family."_homology.txt"){
		
			##Translates nuc to AA
			print "doing alignment for $family\n";
			system("transeq -sequence ".$family.".sp.fa -outseq ".$family.".sp.pep");

			##Removes _1; formatting
			system("sed 's/_1\$//g' ".$family.".sp.pep > ".$family.".sp.faa");
			# system("rm ".$family.".pep");
			
			##MUSCLE: protein msa
			system("muscle3.8.31 -in ".$family.".sp.faa -out ".$family.".sp.afa");
			# # system("rm ".$family.".faa");
			
			##PAL2NAL: convert protein msa to nucleotide 
			system("perl /nv/hp10/lchau6/data/bin/pal2nal.v14/pal2nal.pl ".$family.".sp.afa ".$family.".sp.fa -output fasta > ".$family.".sp.codon.aln");
			
			##TreeBest: build gene trees guided by a Species Tree for detecting homology
			print "Creating species guided gene tree for $family\n"; 
			system("/nv/hp10/lchau6/data/bin/treebest/treebest best -f $speciesTree ".$family.".sp.codon.aln > ".$family.".nhx");
			
			##Python script that uses the ete2 module to find orthologs and paralogs in gene trees
			print "Parsing Genetree $family\n"; 
			system("python /nv/hp10/lchau6/data/bin/geneTreeParse.py ".$family.".nhx ".$family."_homology.txt");
							
		}
	}; 
	$switch1=1; 
}

###Joins running threads 
while(1){
	my @running = threads->list(threads::running);
	my @joinable = threads->list(threads::joinable);
	$_->join() for(@joinable); 
	last unless (@running > 0 or @joinable > 0);
	sleep(5); 
}

###################PAML###################
# my $switch2 = 0; 

# for my $family (@Familylist){
	# if($switch2 ==1){
		# while(1){ ### loop that runs infinitely until one of the running threads completes.  waits N seconds between each loop
			# my @running = threads->list(threads::running);
			# my @joinable = threads->list(threads::joinable);
			# $_->join() for(@joinable);
			# last if(@joinable > 0 or @running < $numProc);  ###exit the stationary loop if any threads are joined (end) in the last pass
			# sleep 3;
		# }
	# }
	# async{
		# unless(-e $family.".codon.aln.gb" && -e $family."_nbs.codon.aln.gb.phylip_phyml_tree.txt"){
			# ##Translates nuc to AA
			# print "doing alignment for $family\n";
			# system("transeq -sequence ".$family.".fa -outseq ".$family.".pep");

			# ##Removes _1; formatting
			# system("sed 's/_1\$//g' ".$family.".pep > ".$family.".faa");
			# # system("rm ".$family.".pep");
			
			# ##MUSCLE: protein msa
			# system("muscle3.8.31 -in ".$family.".faa -out ".$family.".afa");
			# # # system("rm ".$family.".faa");
			
			# ##PAL2NAL: convert protein msa to nucleotide 
			# system("perl /nv/hp10/lchau6/data/bin/pal2nal.v14/pal2nal.pl ".$family.".afa ".$family.".fa -output fasta > ".$family.".codon.aln");		
			
			# ##Use Gblocks to trim poorly aligned segments of the msa
			# system("Gblocks ".$family.".codon.aln -t=c -b5=n -e=.gb -s=y -p=y -n=n -k=y");
			# system("rm ".$family.".codon.aln.gbMask"); 
			# system("rm ".$family.".codon.aln.gb.htm"); 

			# ##Use readseq to convert Gblocks output to phylip format (phylip interleaved)
			# system("java -cp /nv/hp10/lchau6/data/bin/readseq.jar run -inform=8 -format=12 ".$family.".codon.aln.gb");

			# ##Use Phyml to produce maximum likelihood phylogeny trees; HKY85(transitions/transversions, unequal base frequencies)
			# system("phyml --input ".$family.".codon.aln.gb.phylip --bootstrap 100 --datatype nt --model HKY85 --no_memory_check");
			# system("rm ".$family.".codon.aln.gb.phylip_phyml_boot_trees.txt"); 
			# system("rm ".$family.".codon.aln.gb.phylip_phyml_boot_stats.txt");
			# system("rm ".$family.".codon.aln.gb.phylip_phyml_stats.txt");	
			
			# ##Remove Bootstrap values from the phylogeny (fucks up paml by using bootstrap labels as sp labels)
			# system("java -jar /nv/hp10/lchau6/data/bin/PareTree1.0.2.jar -nbs -f ".$family.".codon.aln.gb.phylip_phyml_tree.txt"); 
		# }
		
		# ##Branch Test to get branch lengths from paml Model=0; NS=0; fized kappa=0; kappa=2; fixed omega=0, omega=1
		# my $seqfile = $family.".codon.aln.gb";
		# my $treefile = $family."_nbs.codon.aln.gb.phylip_phyml_tree.txt"; ##Tree with bootstrap values removed 
		# my $null_outctrl = $family."_nullbranch.ctl";
		# my $null_outpaml = $family."_nullbranch.mlc";
		# my $null_model = 0;
		# my $null_NSsites = 0;
		# my $null_fixedKappa=0; 
		# my $null_kappa=2;
		# my $null_fixedomega =0; 
		# my $null_omega =1;
		# &codeml($seqfile ,$treefile, $null_outpaml,$null_outctrl, $null_model, $null_NSsites, $null_fixedKappa, $null_kappa, $null_fixedomega, $null_omega);

		# system("grep -P '\\([A-Za-z].+\\);' ".$null_outpaml." > ".$family."_paml_tree.txt"); ###Pull paml tree from paml output file

		# my $paml_tree = $family."_paml_tree.txt";

		# system("mv ".$family."_nullbranch.mlc ".$paml_nulldir.""); ##Move null paml result to null directory
		# system("mv ".$family."_nullbranch.ctl ".$paml_nulldir.""); ##ctrl file

		# ##Branch Test (Free Ratio; Using tree generated by Paml) Model=1; NS=0; fized kappa=0; kappa=2; fixed omega=0, omega=1
		# my $freeRatio_outctrl = $family."_freeRatio.ctl";
		# my $freeRatio_outpaml = $family."_freeRatio.mlc";
		# my $freeRatio_model = 1;
		# my $freeRatio_NS = 0;
		# my $freeRatio_fixedK=0; 
		# my $freeRatio_k=2;
		# my $freeRatio_fixedw =0; 
		# my $freeRatio_w =1;
		# &codeml($seqfile ,$paml_tree, $freeRatio_outpaml, $freeRatio_outctrl, $freeRatio_model, $freeRatio_NS, $freeRatio_fixedK, $freeRatio_k, $freeRatio_fixedw, $freeRatio_w );

		# system("mv ".$family."_freeRatio.mlc ".$paml_freeratio_dir.""); ##Move free ratio paml result to null directory
		# #system("rm ".$family."_freeRatio.ctl"); ##rm free ratio ctrl file

		# ##Branch-Site Test--->NULL MODEL (model a1: model = 2, NSsites = 2, fix_omega = 1, omega = 1)PAML branch site tree 
		# system("sed 's/".$foreground."[0-9]\{5\}/\#1 &/' ".$family."_paml_tree.txt > ".$family."_bspaml_tree.txt"); #add #1 to the foreground branch
		# my $BSpamltree = $family."_bspaml_tree.txt"; 

		# my $BSnull_outctrl = $family."_BSnull.ctl";
		# my $BSnull_outpaml = $family."_BSnull.mlc";
		# my $BSnull_model = 2;
		# my $BSnull_NS = 2;
		# my $BSnull_fixedK=0; 
		# my $BSnull_k=2;
		# my $BSnull_fixedw =1; 
		# my $BSnull_w =1;

		# &codeml($seqfile ,$BSpamltree, $BSnull_outpaml, $BSnull_outctrl, $BSnull_model, $BSnull_NS, $BSnull_fixedK, $BSnull_k, $BSnull_fixedw, $BSnull_w);

		# system("mv ".$family."_BSnull.mlc ".$paml_BSnull_dir.""); ##Move free ratio paml result to null directory
		# #system("rm ".$family."_BSnull.ctl"); ##rm free ratio ctrl file

		# ##Branch-Site Test--->POSITIVE SELECTION MODEL (model a: model = 2, NSsites = 2, fix_omega = 0)

		# my $BSpos_outctrl = $family."_BSpos.ctl";
		# my $BSpos_outpaml = $family."_BSpos.mlc";
		# my $BSpos_model = 2;
		# my $BSpos_NS = 2;
		# my $BSpos_fixedK=0; 
		# my $BSpos_k=2;
		# my $BSpos_fixedw =0; 
		# my $BSpos_w =1;

		# &codeml($seqfile ,$BSpamltree, $BSpos_outpaml, $BSpos_outctrl, $BSpos_model, $BSpos_NS, $BSpos_fixedK, $BSpos_k, $BSpos_fixedw, $BSpos_w);

		# system("mv ".$family."_BSpos.mlc ".$paml_BSpos_dir.""); ##Move free ratio paml result to null directory
		# system("rm ".$family."_BSpos.ctl"); ##rm free ratio ctrl file
# }	

# $switch2=1; }

# ####Joins running threads 
# while(1){
	# my @running2 = threads->list(threads::running);
	# my @joinable2 = threads->list(threads::joinable);
	# $_->join() for(@joinable2); 
	# last unless (@running2 > 0 or @joinable2 > 0);
	# sleep(5); 
# }

############################ Classify orthologs and paralogs #######################
my (%ortholog_one2one, %ortholog_one2many, %ortholog_many2many, %withinSpecies_paralog, %BtwnSpecies_paralog)=(); 

my @homologyFiles=glob("*_homology.txt"); #pulls all fasta files into an array 

for my $file (@homologyFiles){
	open (INPUT, $file) || die "$file file is not found!\n";
	my $geneTreeID = ($file=~ s/homology.txt//g); 			###geneTree id
	while (<INPUT>){ 	##Read in homology files
		chomp $_;
		if ($_ =~ /^ORTHOLOGY RELATIONSHIP\:/){	###Take in orthologs and classifys them into orthology type 
			$_ =~ s/ORTHOLOGY RELATIONSHIP\://;
			my ($set1, $set2) = split("<====>", $_); 	##seperates the two groups of orthologs
			my @OrthoArray1 = split(",", $set1); 
			my @OrthoArray2 = split(",", $set2);
			if( scalar @OrthoArray1 == 1 && scalar @OrthoArray2 == 1 ){			###Find one to one orthologs
				my $var = $geneTreeID."_".$OrthoArray1[0]; 
				$ortholog_one2one{$var} = $OrthoArray2[0]; 	
			}elsif(( scalar @OrthoArray1 == 1 || scalar @OrthoArray2 == 1) && scalar @OrthoArray1 != scalar @OrthoArray2 ){	###Find one to many orthologs
				unshift(@OrthoArray1, $geneTreeID); 
				my $string1 = join(',', @OrthoArray1); 
				my $string2 = join(',', @OrthoArray2);
				$ortholog_one2many{$string1} = $string2; 
			}else{										###FInds many to many orthologs
				unshift(@OrthoArray1, $geneTreeID);
				my $string3 = join(',', @OrthoArray1);
				my $string4 = join(',', @OrthoArray2);			
				$ortholog_many2many{$string3} = $string4; 
			}
		}elsif($_ =~/^PARALOGY RELATIONSHIP:(.*)/){			###Take in paralogs and classifys them into paralog type
			$_ =~ s/PARALOGY RELATIONSHIP://;
			my ($Pset1, $Pset2) = split("<====>", $_); 	##seperates the two groups of paralogs
			my @ParaArray1 = split(",", $Pset1); 
			my @ParaArray2 = split(",", $Pset2);
			my (%paraHash1, %paraHash2)=(); 
			my $switchP =0; 		###count used to determine if paralog species are the same are not in the parawise comparison
			foreach (@ParaArray1){		###seperates geneID from speciesID, create a AoH with array = paralog set; hash element,geneID=species
				my ($geneID1, $species1) = split(/_([^_]+)$/, $_);
				$paraHash1{$geneID1} = $species1; 
			}
			
			foreach (@ParaArray2){		###seperates geneID from speciesID, create a AoH with array = paralog set; hash element,geneID=species
				my ($geneID2, $species2) = split(/_([^_]+)$/, $_);
				$paraHash2{$geneID2} = $species2; 
			}
			
			for my $paraID_1(keys %paraHash1){			####Compare the species of the two paralog sets
				for my $paraID_2(keys %paraHash2){
					if ($paraHash1{$paraID_1} ne $paraHash2{$paraID_2}){			###If species aren't the same between the two sets,add 1 to switch
						$switchP +=1; 
					}
				}
			}
			if ($switchP > 0){	###Classifies the paralog set as between species (switch >0) and witihin species(switchP =0)
				$BtwnSpecies_paralog{$Pset1}=$Pset2;
			}else{
				$withinSpecies_paralog{$Pset1}=$Pset2;
			}
			
		}
	}close(INPUT);
}

###Print out result files
open (OUTPUT1, ">ortholog_one2one.txt") || die $!;
while(my($key1, $val1) = each %ortholog_one2one){
	print OUTPUT1 "$key1 ;; $val1\n"; 
}close(OUTPUT1); 

open (OUTPUT2, ">ortholog_one2many.txt") || die $!;
while(my($key2, $val2) = each %ortholog_one2many){
	print OUTPUT2 "$key2 ;; $val2\n"; 
}close(OUTPUT2); 

open (OUTPUT3, ">ortholog_many2many.txt") || die $!;
while(my($key3, $val3) = each %ortholog_many2many){
	print OUTPUT3 "$key3 ;; $val3\n"; 
}close(OUTPUT3); 

open (OUTPUT4, ">paralog_withinSpecies.txt") || die $!;
while(my($key4, $val4) = each %withinSpecies_paralog){
	print OUTPUT4 "$key4 ;; $val4\n"; 
}close(OUTPUT4); 

open (OUTPUT5, ">paralog_BtwnSpecies.txt") || die $!;
while(my($key5, $val5) = each %BtwnSpecies_paralog){
	print OUTPUT5 "$key5 ;; $val5\n"; 
}close(OUTPUT5); 

######Subroutines

sub codeml{
my ($seqfile ,$treefile, $outpaml, $outctrl, $model, $NSsites, $fixedKappa, $kappa, $fixedomega, $omega)=@_;

open(CTRL, ">$outctrl");
	print CTRL 
	 "seqfile = $seqfile
     treefile = $treefile
      outfile = $outpaml

        noisy = 0
      verbose = 0
      runmode = 0

      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0   * 0:no clock, 1:global clock; 2:local clock; 3:TipDate
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
   aaRatefile = /../wag.dat

        model = $model
                   * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = $NSsites  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
    fix_kappa = $fixedKappa  * 1: kappa fixed, 0: kappa to be estimated
        kappa = $kappa  * initial or fixed kappa
    fix_omega = $fixedomega  * 1: omega or omega_1 fixed, 0: estimate 
        omega = $omega * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * different alphas for genes
        ncatG = 10  * # of categories in dG of NSsites models

        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

   Small_Diff = .5e-6
*    cleandata = 1
*       method = 1   * 0: simultaneous; 1: one branch at a time";

	close (CTRL);

system("codeml ".$outctrl.""); 
}

sub UsageandExit{
	print "Usage:\n";
	print "./parse_MCL_outfile_lineageSpecificDups.pl <mclOutput> <genefamily Prefix> <# processors> <directory for gene families> <single copy gene File> <species tree> <Species 1> <Species 2> <Species 3> <sp1 seqID prefix> <sp2 seqID prefix> <sp3 seqID prefix> <paml foreground species>\n";
    exit;
	}
	