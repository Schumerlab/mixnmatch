#perl! -w

my $infile=shift(@ARGV); chomp $infile;
open IN, $infile or die "cannot open configuration file\n";

my $genome1=""; my $genome2=""; my $mixture_proportion=""; my $rec_rate_Morgans=""; my $num_indivs=""; my $gens_since_admix=""; my $rate_shared_poly=""; my $per_bp_indel=""; my $number_reads=""; my $gens_drift_parental=""; my $sequence_error=""; my $read_type=""; my $read_length=""; my $poly_par1=""; my $poly_par2=""; my $chr=""; my $freq_diff=""; my $job_submit_cmd=""; my $job_params=""; my $num_indiv_per_job=""; my @commands_array=();

while(my $line = <IN>){

    chomp $line;
    my @elements=split(/=/,$line);

    if($line =~ /genome1=/g){
        $genome1=$elements[1]; chomp $genome1;
        print "parent genome 1 is $genome1\n";
    }#define genome1
    if($line =~ /genome2=/g){
        $genome2=$elements[1]; chomp $genome2;
        print "parent genome 2 is $genome2\n";
    }#define genome2
    if($line =~ /mixture_prop_par1=/g){
	$mixture_proportion=$elements[1]; chomp $mixture_proportion;
    }#define mixture proportion
    if($line =~ /rec_rate_Morgans_kb=/g){
	$rec_rate_Morgans=$elements[1]; chomp $rec_rate_Morgans;
	print "recombination rate is $rec_rate_Morgans\n";
    }#define recombination rate
    if($line =~ /num_indivs/g){
	$num_indivs=$elements[1]; chomp $num_indivs;
	print "simulating data from $num_indivs hybrid individuals\n";
    }#number of individuals to simulate
    if($line =~ /gens_since_admixture=/g){
	$gens_since_admix=$elements[1]; chomp $gens_since_admix;
	print "simulating $gens_since_admix generations since admixture\n";
    }#generations since initial admixture
    if($line =~ /chr_to_simulate=/g){
	$chr=$elements[1]; chomp $chr;
	print "simulating data from chromosome $chr\n";
    }#chromosome to simulate
    if($line =~ /rate_shared_poly_at_aims/g){
	$rate_shared_poly=$elements[1]; chomp $rate_shared_poly;
	print "shared polymorphism rate at AIMs is $rate_shared_poly\n";
    }#shared poly rate
    if($line =~/poly_perbp_par1/g){
	$poly_par1=$elements[1]; chomp $poly_par1;
    }#par1 poly rate
    if($line =~/poly_perbp_par2/g){
	$poly_par2=$elements[1]; chomp $poly_par2;
    }#par2 poly rate 
    if($line =~ /read_type=/g){
        $read_type=$elements[1]; chomp $read_type;
        if(($read_type ne 'SE') && ($read_type ne 'PE')){
            die "read type must be SE or PE\n";
        }
        print "read type is $read_type\n";
    }#define read type
    if($line =~ /per_bp_indels=/g){
	$per_bp_indel=$elements[1]; chomp $per_bp_indel; 
    }#define indel rate
    if($line =~ /number_reads=/g){
	$number_reads=$elements[1]; chomp $number_reads;
	print "simulating $number_reads reads per individual\n";
    }#define number of reads
    if($line =~ /read_length=/g){
	$read_length=$elements[1]; chomp $read_length;
    }#define read length
    if($line =~ /job_submit_cmd=/g){
	$job_submit_cmd=$elements[1]; chomp $job_submit_cmd;
	print "job submit command is $job_submit_cmd\n";
    }#job submission command
    if($line =~ /job_header=/g){
	@commands_array=split(/#/,$line);
     	#print "job submission parameters are $job_params\n";
    }#job header
    if($line =~ /num_indiv_per_job=/g){
	$num_indiv_per_job=$elements[1]; chomp $num_indiv_per_job;
	print "splitting into $num_indiv_per_job individuals per job\n";
    }#num indiv per job 
    if($line =~ /gens_drift_parental=/g){
	$gens_drift_parental=$elements[1]; chomp $gens_drift_parental;
    }#generations drift
    if($line =~ /sequencing_error=/g){
	$sequence_error=$elements[1]; chomp $sequence_error;
    }#error rate
    if($line =~ /aim_freq_cutoff=/g){
	$freq_diff=$elements[1]; chomp $freq_diff;
    }#aim frequency difference imposed in parents

}#for all lines in the infile


my $snp_freqs="shared_polymorphism_distribution";
my $start=1;
my $stop=1;
my $error=1-$freq_diff;

#catalog AIMs, extract focal chromosome
my $chr1="$chr"."_select_par1.fa";
my $chr2="$chr"."_select_par2.fa";
system("perl getScaffold_samtools.pl $genome1 $chr > $chr1");
system("perl getScaffold_samtools.pl $genome2 $chr > $chr2");

my $chr_length=qx(cat $chr1 | tail -n +2 | perl -p -e 's/\n//g' | wc -c | perl -p -e 's/ +/\t/g' | cut -f 1); chomp $chr_length;

print "number of basepairs in $chr is $chr_length\n";

my $aims="simulation_ancestry_informative_sites_"."$genome1"."_"."$genome2"."_"."$chr";
system("perl identify_AIMs_sims.pl $chr1 $chr2 > $aims");

#generate frequency file based on rate of shared polymorphism                                                                                     
#print "Rscript poly_dist.R $rate_shared_poly $poly_par1 $poly_par2 $chr_length $num_indivs $aims $error\n";
system("Rscript poly_dist.R $rate_shared_poly $poly_par1 $poly_par2 $chr_length $num_indivs $aims $error");

#remove previous fastahack index files
if(-e '*.fai'){
    system("rm *.fai");
}#remove

open LIST, ">simulated_hybrids_readlist_"."gen"."$gens_since_admix"."_prop_par1_"."$mixture_proportion";
my $reads_folder="simulated_hybrids_reads_"."gen"."$gens_since_admix"."_prop_par1_"."$mixture_proportion";

if(-e $reads_folder){
    print "WARNING: removing previous results in $reads_folder\n";
    system("rm -r $reads_folder");
}#remove previous results
system("mkdir $reads_folder");

#define files and counters
my $current_outfile="split_file_list_1"; my $current_slurm="slurm_batch1.sh"; my $counter=0; my $track=0;

while($counter<$num_indivs){
    $counter=$counter+1; $track=$track+1;
    #print "$counter\t$track\n";
    
    if(($track>=$num_indiv_per_job) or ($counter eq $num_indivs) or ($counter eq 1)){

	if($track eq $num_indiv_per_job){
	
	    #print commands for current job
	    print SLURM "perl generate_genomes_and_reads.pl $genome1 $genome2 $chr $chr1 $chr2 $chr_length $poly_par1 $poly_par2 $aims $error $reads_folder $mixture_proportion $gens_since_admix $rec_rate_Morgans $snp_freqs $current_outfile $number_reads $read_length $sequence_error $per_bp_indel\n";
	    
	    #submit current job
	    system("$job_submit_cmd $current_slurm");
	}#submit job
	$current_outfile="split_file_list_"."$counter";
	$current_slurm="slurm_batch"."$counter".".sh"; 
	open OUT, ">$current_outfile";
	open SLURM, ">$current_slurm";
	for my $b (1..scalar(@commands_array)-1){
	print SLURM "#"."$commands_array[$b]"."\n";
	}#print out all slurm header elements
	$track=0;
	}#open individual files

    #print to individuals file for the shell
     print "submitting individual $counter\n";
     print OUT "$counter\n";

	my $r1="$reads_folder"."/"."indiv"."$counter"."_read1.fq"; my $r2="$reads_folder"."/"."indiv"."$counter"."_read2.fq";
 
    print LIST "$r1".".gz"."\t$r2".".gz"."\n";
 
}#for all jobs
