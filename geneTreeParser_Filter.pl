#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper; 

my @homologyFiles=glob("*_homology.txt"); #pulls all fasta files into an array 

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
