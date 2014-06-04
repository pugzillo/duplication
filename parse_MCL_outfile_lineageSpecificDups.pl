#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper; 
use Cwd;
use threads;

###Perl script to find lineage specific duplicates in orthoMCL output, align them with muscle, and then use paml yn00 to calculate dS

################################ VARIABLES ###################################################
###Input variables
if ($#ARGV != 4 ) {;}
my $MCLout = $ARGV[0]; 			##mcl output from blast all v all and then single linkage clustering with MCL
my $genePair = $ARGV[1];		##Prefix for gene pair
my $numProc = $ARGV[2];			##Number of processors
my $orthodir = $ARGV[3];		##Name of directory for all gene family files
my $SC_file = $ARGV[4]; 		###File for single copy genes
my $speciesTree = $ARGV[5];		###File with Species tree for TreeBeST

###Variables
my $seq_id;
my $seq; 
my $cwd=getcwd();
my $count = 1; 

###Data Structures
my (%family, %sequence, @Familylist, %single_copy, @singleCopy, )=();

mkdir $orthodir;

#Populate Hashes with Sequence Data for each gene set; make sure that there is no sequence wrapping!!!!
my @cds_files=glob("*.fa"); #pulls all fasta files into an array 

################# READ IN SEQUENCE FILES AND MCL CLUSTERING OUTPUT FILES ###################
for my $fasta(@cds_files){			####Reads through each fasta file and creates hash with gene id and the sequence
	open (IN, $fasta) || die "$fasta file is not found!\n";
	print "Reading in $fasta\n";
	while (<IN>){ #Read in multifasta file
		chomp $_;
		if ($_ =~ /^>(.*)/){
			$seq_id = $_; #Put the sequence ID here!
			$seq_id =~ s/\W//g;
			$seq = $1; 
		} else{
			$seq =$_; #Sequence goes here!
		}
		$sequence{$seq_id} = $seq; ###Hash with the id and sequence for each gene in the genome file
	}close(IN);
}

open (IN, $MCLout) || die "$MCLout is not found!\n";
print "Reading in $MCLout\n";
while (<IN>){ #Read in the MCL output lines with only two genes or more in the gene family
	chomp $_;
	my @geneFam = split("\t", $_); 
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
}close(IN);

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
			print "doing alignment and paml";
			system("transeq -sequence ".$family.".fa -outseq ".$family.".pep");

			##Removes _1; formatting
			system("sed 's/_1//g' ".$family.".pep > ".$family.".faa");

			##Removes in sequence stop codons
			system("perl /nv/hp10/lchau6/data/apis_molevo/testing/thread_testing/testStop.pl --pep ".$family.".faa --nuc ".$family.".fa --output ".$family.".nostop.fa");
			
			##PRANK: codon msa
			system("prank -d=".$family.".nostop.fa -o=".$family."_codon -codon -F");
			system("prank -convert -d=".$family."_codon.best.fas -f=paml -o=".$family." -keep");
			
			##TreeBest: build gene trees guided by a Species Tree
			system("treebest best -f $speciesTree ".$family."_codon.best.fas > ".$family.".nhx");
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




