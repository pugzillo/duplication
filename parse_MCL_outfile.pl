#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper; 


my (%family, %sequence, @Familylist)=();

#Populate Hashes with Sequence Data for each gene set; make sure that there is no sequence wrapping!!!!
my @cds_files=glob("*.fa"); #pulls all fasta files into an array

for my $files(@cds_files){ #Loops through each sequence file
	print "$files sequence file is being read\n";
	open (IN, $files) || die "$files sequence file is not found!\n";
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

#Populate Array with MCL clustering output files
my @clusterFiles=glob("*.mci.I20"); 

for my $file(@clusterFiles){ #Loops through each MCL output file
	print "$file MCL output is being read\n";
	open (IN, $file) || die "$file is not found!\n";
	while (<IN>){ #Read in the MCL output
		chomp $_;
		my @geneFam = split("\t", $_); 
		my $FamID = ; 		##create random id for each gene family
		for my $geneID(@geneFam){
			if (exists $sequence{$geneID}){
				$family{$FamID} = $geneID."\t".$sequence{$geneID}; ####Hash with gene family ID as the key and the geneIDs and sequences as the values delimited by a new line
		}
	}close(IN);
}
}

#Print out a multi-fasta file for each gene Family
for my $fam (keys %family){
	open (OUT, ">", $fam.".fa") || die "$ortho outfile cannot be opened!\n";
	push @Familylist,$fam. ".fa"; #Array will the list of ortholog files
	for my $seq_id (keys %{$family{$fam}}){
		my ($id, $sequence) = split("\t", $seq_ID);
		print OUT ">$id\n$sequence\n";	
	}
}close(OUT);

