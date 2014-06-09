#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper; 
use Cwd;
use threads;

###Perl script to find lineage specific duplicates in orthoMCL output, align them with PRANK, create species guided phylogenies with TreeBest, and parse each tree to find orthologs and paralogs for each gene family.  
##make sure both perl and python are loaded!!!

################################ VARIABLES ###################################################
###Input variables
if ($#ARGV != 11 ) {;}
my $MCLout = $ARGV[0]; 			##mcl output from blast all v all and then single linkage clustering with MCL
my $genePair = $ARGV[1];		##Prefix for gene pair
my $numProc = $ARGV[2];			##Number of processors
my $orthodir = $ARGV[3];		##Name of directory for all gene family files
my $SC_file = $ARGV[4]; 		###File for single copy genes
my $speciesTree = $ARGV[5];		###File with Species tree for TreeBeST
my $species1 = $ARGV[6];
my $species2 = $ARGV[7];
my $species3 = $ARGV[8];
my $sp1 = $ARGV[9]; 
my $sp2 = $ARGV[10]; 
my $sp3 = $ARGV[11]; 

###Variables
my $seq_id;
my $seq; 
my $cwd=getcwd();
my $count = 1; 

###Data Structures
my (%family, %sequence, @Familylist, %single_copy, @singleCopy, @homologyFiles)=();

mkdir $orthodir;

#Populate Hashes with Sequence Data for each gene set; make sure that there is no sequence wrapping!!!!
my @cds_files=glob("*.fa"); #pulls all fasta files into an array 

################# READ IN SEQUENCE FILES AND MCL CLUSTERING OUTPUT FILES ##############
for my $fasta(@cds_files){			####Reads through each fasta file and creates hash with gene id and the sequence
	open (IN, $fasta) || die "$fasta file is not found!\n";
	print "Reading in $fasta\n";
	while (<IN>){ #Read in multifasta file
		chomp $_;
		my @fileShit = split("_", $fasta);
		my $species = $fileShit[0];			###extract species name from fasta file name 
		if ($_ =~ /^>(.*)/){
			$seq_id = $_; #Put the sequence ID here!
			$seq_id =~ s/\W//g;
			$seq_id = $seq_id."_".$species; 
			$seq = $1; 
		} else{
			$seq =$_; #Sequence goes here!
		}
		$sequence{$seq_id} = $seq; ###Hash with the id and sequence for each gene in the genome file
	}close(IN);
}

open (INPUT, $MCLout) || die "$MCLout is not found!\n";
print "Reading in $MCLout\n";
while (<INPUT>){ #Read in the MCL output lines with only two genes or more in the gene family
	chomp $_;
	my @geneFam = split("\t", $_); 
	foreach my $member(@geneFam){
		if ($member =~ $sp1){
			$member = $member."_".$species1; 
		}
		if ($member =~ $sp2){
			$member = $member."_".$species2; 
		}
		if ($member =~ $sp3){
			$member = $member."_".$species3; 
		}
	}
	if (scalar @geneFam >= 2){
		my $FamID = $genePair."_".$count; 
		for my $geneID(@geneFam){
			if (exists $sequence{$geneID}){
				my $line = $geneID."\t".$sequence{$geneID}; ####Hash with gene family ID as the key and the geneIDs and sequences as the values delimited by a new line
				push @{$family{$FamID}}, $line;  
			}
		}
		$count++; 
	}else{
		for my $singleGene(@geneFam){ ###An array for single copy genes (gene family n==1) along with sequence
			if (exists $sequence{$singleGene}){
				my $singleLine = $singleGene."\t".$sequence{$singleGene};
				push @singleCopy, $singleLine; 
			}
		}
	}
}close(INPUT);

#####Print out a multi-fasta file for each gene Family
for my $fam (keys %family){
	open (OUT, ">", $orthodir."/".$fam.".fa") || die "$fam outfile cannot be opened!\n";
	push @Familylist,$fam; #Array will the list of gene families
	foreach(@{$family{$fam}}){
		my ($id, $sequence) = split("\t", $_);
		print OUT ">$id\n$sequence\n";	
	}
}close(OUT);

#####Print out file with single copy genes
open (OUTPUT, ">$SC_file") || die $!;
for my $singleCopyGene (@singleCopy){
	print OUT $singleCopyGene ."\n";
	}
close (OUTPUT);

###Change to orthofile directory
chdir ($orthodir); 

##################################### MSA AND PHYLOGENY FOR EACH GENE FAMILY ##############################################
my $switch = 0; 

for my $family (@Familylist){
	if($switch ==1){
		while(1){ ### loop that runs infinitely until one of the running threads completes.  waits N seconds between each loop
			my @running = threads->list(threads::running);
			my @joinable = threads->list(threads::joinable);
			$_->join() for(@joinable);
			last if(@joinable > 0 or @running < $numProc);  ###exit the stationary loop if any threads are joined (end) in the last pass
			sleep 3;
		}
	}
	async{
		unless(-e $family."_codon_interleaved.best.fas.gb.phylip"){
			##Translates nuc to AA
			print "doing alignment for $family\n";
			system("transeq -sequence ".$family.".fa -outseq ".$family.".pep");

			##Removes _1; formatting
			system("sed 's/_1//g' ".$family.".pep > ".$family.".faa");

			##Removes in sequence stop codons
			system("perl /nv/hp10/lchau6/data/apis_molevo/testing/thread_testing/testStop.pl --pep ".$family.".faa --nuc ".$family.".fa --output ".$family.".nostop.fa");
			
			##PRANK: codon msa
			system("prank -d=".$family.".nostop.fa -o=".$family."_codon -codon -F");
			system("prank -convert -d=".$family."_codon.best.fas -f=paml -o=".$family." -keep");
			
			##TreeBest: build gene trees guided by a Species Tree
			print "Creating species guided gene tree for $family\n"; 
			system("treebest best -f $speciesTree ".$family."_codon.best.fas > ".$family.".nhx");
			
			##Python script that uses the ete2 module to find orthologs and paralogs in gene trees
			print "Parsing Genetree $family\n"; 
			system("python geneTreeParse.py ".$family.".nhx ".$family."_homology.txt");
			
			##An array with all the homology files
			push(@homologyFiles, $family."_homology.txt");
		}

	}; 
	$switch=1; 
}
###Joins running threads 
while(1){
	my @running = threads->list(threads::running);
	my @joinable = threads->list(threads::joinable);
	$_->join() for(@joinable); 
	last unless (@running > 0 or @joinable > 0);
	sleep(5); 
}

############################# Classify orthologs and paralogs #######################
my (%ortholog_one2one, %ortholog_one2many, %ortholog_many2many, %withinSpecies_paralog, %BtwnSpecies_paralog)=(); 

for my $file (@homologyFiles){
	open (INPUT, $file) || die "$file file is not found!\n";
	while (<INPUT>){ 	##Read in homology files
		chomp $_;
		if ($_ =~ /^ORTHOLOGY RELATIONSHIP\:/){	###Take in orthologs and classifys them into orthology type 
			$_ =~ s/ORTHOLOGY RELATIONSHIP\://;
			my ($set1, $set2) = split("<====>", $_); 	##seperates the two groups of orthologs
			my @OrthoArray1 = split(",", $set1); 
			my @OrthoArray2 = split(",", $set2);
			if( scalar @OrthoArray1 == 1 && scalar @OrthoArray2 == 1 ){			###Find one to one orthologs
				$ortholog_one2one{$OrthoArray1[0]} = $OrthoArray2[0]; 	
			}elsif(( scalar @OrthoArray1 == 1 || scalar @OrthoArray2 == 1) && scalar @OrthoArray1 != scalar @OrthoArray2 ){	###Find one to many orthologs
				my $string1 = join(',', @OrthoArray1);
				my $string2 = join(',', @OrthoArray2);
				$ortholog_one2many{$string1} = $string2; 
			}else{										###FInds many to many orthologs
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
			my $switch =0; 				###count used to determine if paralog species are the same are not in the parawise comparison
			foreach (@ParaArray1){		###seperates geneID from speciesID, create a AoH with array = paralog set; hash element,geneID=species
				my ($geneID1, $species1) = split("_", $_);
				$paraHash1{$geneID1} = $species1; 
			}
			
			foreach (@ParaArray2){		###seperates geneID from speciesID, create a AoH with array = paralog set; hash element,geneID=species
				my ($geneID2, $species2) = split("_", $_);
				$paraHash2{$geneID2} = $species2; 
			}
			
			for my $paraID_1(keys %paraHash1){			####Compare the species of the two paralog sets
				for my $paraID_2(keys %paraHash2){
					if ($paraHash1{$paraID_1} ne $paraHash2{$paraID_2}){			###If species aren't the same between the two sets,add 1 to switch
						$switch +=1; 
					}
				}
			}
			if ($switch > 0){			###Classifies the paralog set as between species (switch >0) and witihin species(switch =0)
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


