#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use Data::Dumper; 

###GeneFamily_BlastMCL_pipeline.pl User notes:
# 
#1. Load in fasta file with all the sequences you want to used to find gene femilies (protein or cds)
#2. Run reciprocal on every gene against every other (both self and non-species) 
#3. Generate clusters using MCL (based off of blast eval), which has little insight about the phylogeny of species
#
#Multi-fasta files must be linearized!
#sed -e 's/\(^>.*$\)/#\1#/' file.fa | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > linearFile.fa

###########Input Options#################
if ($#ARGV != 0 ) {;}
my $file = $ARGV[0];  ####Fasta file with sequences you want to use to find gene families 
my $numProc =$ARGV[1]; #### number of processors to run with

my $cwd=getcwd();			 #Current directory
my (%options, )= ();
my ($seqtype, $prg) = (); 
my (@blastOUT, @ClusterFiles) = ();

getopts("pn", \%options);	 #enter if protein (-p) or nucleotide (-n)

##takes in user option of protein or nucleotide blast
 if (defined $options{p}){  
 	$seqtype = "prot"; 
 	$prg = "blastp"; 
 }else{
 	$seqtype = "nucl"; 
 	$prg= "blastn"; 
 }

#################Directories##################

my $blastResults = $cwd . "/BLASTfiles"; 
mkdir $blastResults;

#######################################BLAST#####################################################

###Reading in Fasta file and making a blast dataset from it
print "$file geneset is being read\n";
if (! -e $file.".nhr" && ! -e $file.".nin" && ! -e $file.".nsq"){
	print "Creating $file BLAST database\n";
	system ("makeblastdb -in $file -dbtype $seqtype"); 
}else{
	print "$file BLAST database exists\n";
}

#############Blast#############
my $blastout = $file."blastOUTPUT.txt"; ####File with blast output

###Checks for the type of BLAST and performs it (eval cutoff og 0.00001)
if ($prg eq "blastn"){
	system("blastn -outfmt=6 -evalue=0.00001 -task=blastn -query ".$file." -db ".$file." -out ".$blastoutput.""); 
}elsif($prg eq "blastp"){
	system("blastp -outfmt=6 -evalue=0.00001 -task=blastp -query ".$file." -db ".$file." -out ".$blastoutput."");
}

system("mv ".$blastoutput." ".$blastResults.""); ###Move blast result files to the blastoutput directory

chdir ($blastResults); 
	
#####Create MCL input with just query and subject and e-value
system("cut -f 1,2,11 ".$blastoutput." > ".$blastoutput."MCLinput.txt"); 


########################MCL clustering#########################################
#####Create network and dictionary file; --abc-neg-log10 tranforms the numerical values in the input (the BLAST E-values) by taking the logarithm in base 10 and subsequently negating the sign; the transformed values are capped so that any E-value below 1e-200 is set to a maximum allowed edge weight of 200
system("mcxload --stream-mirror -abc ".$blastoutput."MCLinput.txt --stream-neg-log10 -stream-tf 'ceil(200)' -o ".$blastoutput."_blastMCL.mci -write-tab ".$blastresult."_blastMCL.tab");

####Clustering; -I is the inflation constant; range from 1.2-5, with 5 having the most fine granularity in clusteringand 1.2 having the largest (check Mammal_ensembl_GeneFams.xlsx for parameter tests)
system("mcl ".$blastoutput."_blastMCL.mci -I 2.1 -tf 'gq(0.7)' -use-tab ".$blastoutput."_blastMCL.tab");

######################################################################################################

sub UsageandExit{
	print "Usage:\n";
	print "./GeneFamily_BlastMCL_pipeline.pl -n/-p <fasta.fa> <numProcessors>\n";
    exit;
	}
	