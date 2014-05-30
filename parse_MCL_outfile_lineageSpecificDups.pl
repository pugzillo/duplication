#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper; 
use Cwd;
use threads;

###Perl script to find lineage specific duplicates in orthoMCL output, align them with muscle, and then use paml yn00 to calculate dS

my (%family, %sequence, @Familylist)=();
my $count = 1; 

if ($#ARGV != 4 ) {;}

my $MCLout = $ARGV[0]; 		##mcl output
my $genePair = $ARGV[1];	##Prefix for gene pair
my $numProc = $ARGV[2];
my $orthodir = $ARGV[3];
my $SC_file = $ARGV[4]; 		###File for single copy genes

my $seq_id;
my $seq; 
my $cwd=getcwd();
my (%sequence, %family, %single_copy, @singleCopy, )=();

mkdir $orthodir;

#Populate Hashes with Sequence Data for each gene set; make sure that there is no sequence wrapping!!!!
my @cds_files=glob("*.fa"); #pulls all fasta files into an array 

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
				push @singleCopy, $line; 
			}
		}
	}
}close(IN);

# #####Print out a multi-fasta file for each gene Family
# for my $fam (keys %family){
	# open (OUT, ">", $orthodir."/".$fam.".fa") || die "$fam outfile cannot be opened!\n";
	# push @Familylist,$fam; #Array will the list of gene families
	# foreach(@{$family{$fam}}){
		# my ($id, $sequence) = split("\t", $_);
		# print OUT ">$id\n$sequence\n";	
	# }
# }close(OUT);

#####Print out file with single copy genes
open (OUTPUT, ">$SC_file") || die $!;
for my $singleCopyGene (@singleCopy){
	print OUT $singleCopyGene ."\n";
	}
close (OUTPUT);

# ##Change to orthofile directory
# chdir ($orthodir); 

# #############################################################################################
# my $switch=0; 

# for my $family (@Familylist){
	# if($switch ==1){
		# while(1){ ### loop that runs infinitely until one of the running threads completes.  waits N seconds between each loop
			# my @running = threads->list(threads::running);
			# my @joinable = threads->list(threads::joinable);
			# $_->join() for(@joinable);
			# last if(@joinable > 0 or @running < $numProc);  ###exit the stationary loop if any threads are joined (end) in the last pass
			# sleep 3;
		# }
	# }
	# async{
		# unless(-e $family."_codon_interleaved.best.fas.gb.phylip"){
			# ##Translates nuc to AA
			# print "doing alignment and paml";
			# system("transeq -sequence ".$family.".fa -outseq ".$family.".pep");

			# ##Removes _1; formatting
			# system("sed 's/_1//g' ".$family.".pep > ".$family.".faa");

			# ##Removes in sequence stop codons
			# system("perl /nv/hp10/lchau6/data/apis_molevo/testing/thread_testing/testStop.pl --pep ".$family.".faa --nuc ".$family.".fa --output ".$family.".nostop.fa");
			
			# ##PRANK: codon msa
			# system("prank -d=".$family.".nostop.fa -o=".$family."_codon -codon -F");
			# system("prank -convert -d=".$family."_codon.best.fas -f=paml -o=".$family." -keep");

		# }

		# ####PAML: yn00
		
		# ##Branch Test to get branch lengths from paml Model=0; NS=0; fized kappa=0; kappa=2; fixed omega=0, omega=1
		# my $seqfile = $family.".phy"; 
		# my $outctrl = $family."_yn00.ctl";
		# my $outpaml = $family."_yn00.result";
		# &yn00($seqfile, $outctrl, $outpaml);

	# }; 
	# $switch=1; 
# }
# ###Joins running threads 
# while(1){
	# my @running = threads->list(threads::running);
	# my @joinable = threads->list(threads::joinable);
	# $_->join() for(@joinable); 
	# last unless (@running > 0 or @joinable > 0);
	# sleep(5); 
# }


sub yn00{
my ($seqfile, $outctrl, $outpaml)=@_;

open(CTRL, ">$outctrl");
	print CTRL 
	 " seqfile = $seqfile * sequence data file name 
	   outfile = $outpaml * main result file 
       verbose = 0 * 1: detailed output (list sequences), 0: concise output 
 
       icode = 0 * 0:universal code; 1:mammalian mt; 2-10:see below 
       weighting = 0 * weighting pathways between codons (0/1)? 
       commonf3x4 = 0 * use one set of codon freqs for all pairs (0/1)?
	   ndata = 1";

	close (CTRL);

system("yn00 ".$outctrl.""); 
}

