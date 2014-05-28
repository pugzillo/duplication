#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use threads;
use Data::Dumper; 

# ###GeneOrthology_pipeline_threaded.pl: 
# #1. Load gene sets from the genomes in which you want to find homologs in
# #2. Run blastall on every gene against every other (both self and non-species) in a genome wise manner
# #3. Generate clusters using MCL, which has little insight about the phylogeny of species


# ####Multi-fasta files must be linearized
# ####sed -e 's/\(^>.*$\)/#\1#/' aflor_testseq.fa | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > aflor_testseq_linear.fa

###Input Options####
if ($#ARGV != 0 ) {;}

my $numProc =$ARGV[0]; #### number of processors to run with

###Variables####
sub maxval;
sub minval;

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

########################################Directories#############################################################

my $blastResults = $cwd . "/BLASTfiles"; 
mkdir $blastResults;

#######################################BLAST_ALL################################################################

#Populate Hashes with Sequence Data for each gene set; make sure that there is no sequence wrapping!!!!
my @cds_files=glob("*.fa"); #pulls all fasta files into an array 

#make blast db for each gene set
for my $file(@cds_files){ 
	print "$file geneset is being read\n";
	if (! -e $file.".nhr" && ! -e $file.".nin" && ! -e $file.".nsq"){
		print "Creating $file BLAST database\n";
		system ("makeblastdb -in $file -dbtype $seqtype"); 
	}else{
		print "$file BLAST database exists\n";
	}
}

###blastall: blast each genome against each other genome
for my $queryFile(@cds_files){
	for my $DBfile(@cds_files){
		my ($Species, $version, $crap) = ($queryFile =~ /(.*)_(.*)_(.*.fa)/ );
		my $querySpecies = $Species."_".$version; 
		print $querySpecies."\n";
		my ($SubSpecies, $dbversion, $dbcrap) = ($DBfile =~ /(.*)_(.*)_(.*.fa)/ );
		my $dbSpecies = $SubSpecies."_".$dbversion; 
		print "Blast $querySpecies against $dbSpecies\n";
		my $blastoutput = $querySpecies."_blast_".$dbSpecies.".txt"; 
		system("".$prg." -outfmt=6 -task=".$prg." -strand plus -query ".$queryFile." -db ".$DBfile." -out ".$blastoutput.""); 
		push (@blastOUT, $blastoutput); 		###An array which includes all of the blastoutput files
		system("mv ".$blastoutput." ".$blastResults.""); ###Move blast result files to the blastoutput directory
}
}

##########Filtering all of the blastall results (critera: $identity >= 50 &&  $percentOverlap >= 0.6 && $eval <= 9e-3)

chdir ($blastResults); 

my $switch1=0; 

for my $blastresult(@blastOUT){
	print "Filtering $blastresult blast results.\n";
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
	##Opens the blast output file and reads into a hash; filters out to keep only top hit for each query-subject match
		open (IN, "<$blastresult") || die $!; 
		my (%id, %criteria, %filtered) = ();
		
		while(<IN>){ 	
			chomp $_;
			my $line = $_; 
			my ($query, $subject, $id_percent, $algn_length, $mismatch_no, $gap_no, $querystart, $queryend, $subjectstart, $subjectend, $evalue, $bit) = split("\t", $line);
			if (exists $id{$query}{$subject} && $id{$query}{$subject}{"eval"} < $evalue){	#for each query-subject pair, it finds the line blast match that has the most significant evalue	
				$id{$query}{$subject}{"line"} = $line;
				$id{$query}{$subject}{"eval"} = $evalue;
			}else{
				$id{$query}{$subject}{"line"} = $line;
				$id{$query}{$subject}{"eval"} = $evalue;
			}
		}close IN;
		
		# #Goes through blast hash; makes it easier to parse later on by simplifying the hash to just the query and the entire line as the value
		for my $value (keys %id){    
			for my $hit (keys %{$id{$value}}){
				my $val = $id{$value}{$hit}{"line"};
				if (exists $filtered{$value}){
					push( @{$filtered{$value}}, $val);
				}else{
					$filtered{$value} = [$val];
				}
			}
		}
		
		# #Loops through the blast hash and finds those blast pairs that are significant
		for my $gene (keys %filtered){ 
			foreach(@{$filtered{$gene}}){
				my ($query_id, $subject_id, $identity, $length, $mismatch, $gap, $q_start, $q_end, $s_start, $s_end, $eval, $bitscore) = split("\t", $_);
				my $longestGene; 
				my $query_length =$q_end-$q_start; 	
				my $subject_length =$s_end-$s_start;
				my $overlap =  minval($s_end, $q_end) - maxval($s_start, $q_start);
				if ($query_length >= $subject_length){	###Finds whether the query or subject is the longest gene
					$longestGene = $query_length;
				}else{
					$longestGene = $subject_length; 
				}
				my $percentOverlap = $overlap / $longestGene; 			###Percent Overlap 
				if ($identity >= 50 &&  $percentOverlap >= 0.6 && $eval <= 9e-3){ ###Filters based on criteria set by Assis and Bachtrog 2013
					push (@{ $criteria{$gene} }, $_); 	####Hash with significant matches
				}
			}
		
		my $filtered_outfile = $blastresult."_filtered.txt"; 
		
		##Print out filtered blast results
		open (OUT, ">$filtered_outfile") || die $!;
		for my $filtered_gene (keys %criteria){
			for my $filtered_val (@{$criteria{$filtered_gene}}){
				print OUT $filtered_val ."\n";
			}
		}close OUT;
		
		#####Create MCL input with just query and subject and e-value
		system("cut -f 1,2,11 ".$filtered_outfile." > ".$blastresult."MCLinput.txt"); 


		########################MCL clustering#########################################
		#####Create network and dictionary file; --abc-neg-log10 tranforms the numerical values in the input (the BLAST E-values) by taking the logarithm in base 10 and subsequently negating the sign; the transformed values are capped so that any E-value below 1e-200 is set to a maximum allowed edge weight of 200
		system("mcxload -abc ".$blastresult."MCLinput.txt --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o ".$blastresult."_blastMCL.mci -write-tab ".$blastresult."_blastMCL.tab");

		####Clustering; -I is the inflation constant; range from 1.2-5, with 5 having the most fine granularity in clusteringand 1.2 having the largest
		system("mcl ".$blastresult."_blastMCL.mci -I 2 -use-tab ".$blastresult."_blastMCL.tab");
	
	 }
		push(@ClusterFiles, "out.".$blastresult."_blastMCL.mci.I20");
}}
		$switch1=1; 
	
# ########Joins running threads 
while(1){
	my @running = threads->list(threads::running);
	my @joinable = threads->list(threads::joinable);
	$_->join() for(@joinable); 
	last unless (@running > 0 or @joinable > 0);
	sleep(5); 
}
######################################################################################################

sub UsageandExit{
	print "Usage:\n";
	print "./genefamily_pipeline_threaded.pl -n/-p <numProcessors> \n";
    exit;
	}
	
sub maxval {
	my ($max_one, $max_two)= @_;
	if ($max_one >= $max_two){
		return $max_one;
	}else{
		return $max_two;
	}
}
	
sub minval {
	my ($min_one, $min_two)= @_;
	if ($min_one <= $min_two){
		return $min_one;
	}else{
		return $min_two;
	}
}