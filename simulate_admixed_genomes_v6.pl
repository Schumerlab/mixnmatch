#perl! -w

my $infile=shift(@ARGV); chomp $infile;
open IN, $infile or die "cannot open configuration file\n";

my $genome1=""; my $genome2=""; my $mixture_proportion=""; my $rec_rate_Morgans=""; my $num_indivs=""; my $gens_since_admix=""; my $rate_shared_poly=""; my $per_bp_indel=""; my $number_reads=""; my $parental_drift=""; my $sequence_error=""; my $read_type=""; my $read_length=""; my $poly_par1=""; my $poly_par2=""; my $chr=""; my $freq_diff=""; my $job_submit_cmd=""; my $job_params=""; my $num_indiv_per_job=""; my @commands_array=(); my $macs=0; my $selam_param=""; my $macs_params=""; my $seq_params=""; my $total_par1=0; my $total_par2=0; my $program_path=""; my $use_map=0; my $recombination_map=""; my $par1_aim=0; my $par2_aim=0; my $use_ancestral=0; my $par1_drift=0; my $par2_drift=0;

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
    if($line =~ /use_ancestral=/g){
	$use_ancestral=$elements[1]; chomp $use_ancestral;
	print "using genome1 as ancestral sequence\n";
    }#use an ancestral sequence with seq-gen?
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
    if($line =~ /parental_drift=/g){
	$parental_drift=$elements[1]; chomp $parental_drift;
    }#add parental drift
    if($line =~ /macs_par1_aims_pop/g){
	my $par1_drift_tmp=$elements[1]; chomp $par1_drift_tmp;
	if(length($par1_drift_tmp)>0){$par1_drift=$par1_drift_tmp}
    }#par1 drift population
    if($line =~ /macs_par2_aims_pop/g){
        my $par2_drift_tmp=$elements[1]; chomp $par2_drift_tmp;
        if(length($par2_drift_tmp)>0){$par2_drift=$par2_drift_tmp}
    }#par2 drift population  
    if($line =~ /sequencing_error=/g){
	$sequence_error=$elements[1]; chomp $sequence_error;
    }#error rate
    if($line =~ /aim_freq_cutoff=/g){
	$freq_diff=$elements[1]; chomp $freq_diff;
    }#aim frequency difference imposed in parents
    if($line =~ /use_macs=/g){
	$macs=$elements[1]; chomp $macs;
    }#use macs or don't
    if($line =~ /use_map=/g){
	$use_map=$elements[1]; chomp $use_map;
    }#use recombination map
    if($line =~ /SELAM_param_file=/g){
	$selam_param=$elements[1]; chomp $selam_param;
	if(length($selam_param) > 0){
	    print "using $selam_param as input to SELAM\n";
	}#print which file is being used
    }#selam param file provided?
    if($line =~ /macs_params=/g){
	$macs_params=$elements[1]; chomp $macs_params;
    }#macs param file
    if($line =~ /seq_params=/g){
	$seq_params=$elements[1]; chomp $seq_params;
    }#sequence parameters
    if($line =~/par1_for_aims/g){
	$par1_aim=$elements[1]; chomp $par1_aim;
    }#how many par1 to use for aims
    if($line =~/par2_for_aims/g){
	$par2_aim=$elements[1];chomp $par2_aim;
    }#how many par2 to use for aims
    if($line =~ /program_path=/g){
	$program_path=$elements[1]; chomp $program_path;
	if(length($program_path)<1){$program_path="./";}
    }#set path for simulator
}#for all lines in the infile

if(($use_macs eq 1) && (length($macs_params) eq 0)){
    die "must provided macs parameters if use_macs is 1\n";
   }#die if use macs is 1 but no parameters are provided

if(($parental_drift eq 1) && ($par1_drift + $par2_drift) eq 0){
    die "must provide population ids for AIMs definition if drift between hybridizing and AIM pops specified\n";
}#parental drift specified but no pop ids specified

#check for SELAM parameters:

if(length($selam_param) > 0){
    system("cp $selam_param selam_demography.txt");
    print "selam param is $selam_param\n";
}#copy to demography file name
else{
    $selam_param="selam_demography.txt";
    open SLP, ">selam_demography.txt";
    print SLP "Pop1\tPop2\tSex\t0\t1\n";
    my $par1mix=$mixture_proportion;
    my $par2mix=1-$par1mix;
    print SLP "0\t0\tA\t5000\t5000\n";
    print SLP "0\ta0\tA\t$par1mix\t0\n";
    print SLP "0\ta1\tA\t$par2mix\t0\n";
}#no parameter file provided, use default 

#set output file for SELAM based on cfg input
if(-e 'admix_simulation_demography_output_results.txt'){system("rm admix_simulation_demography_output_results.txt");}
#remove previous results since SELAM appends

open SELAMOUT,">selam_simulation_output_parameters.txt";
$num_per_sex=$num_indivs;
print SELAMOUT "$gens_since_admix\t0\t$num_per_sex\t$num_per_sex\tadmix_simulation_demography_output_results.txt\n";

#set parameters used later
my $snp_freqs="shared_polymorphism_distribution";
my $error=1-$freq_diff;

#run macs here if needed, convert to sequences for downstream
if($macs eq 1){

    my $total=2*$num_indivs;
    
    $macs_params=~ s/ +/ /g;

    my @splitmacsparams=split(/ /,$macs_params);
    my $trackid=2; my $lengthid=0; my $thetaid=1; my $simchr=""; my $simtheta=0; my $recid=1; my $total_parental_haps=0; my $pop_tracker=0; my $focal_prev=""; my @track_pop=();
    for my $p (0..scalar(@splitmacsparams)-1){
	my $focal=$splitmacsparams[$p]; chomp $focal;
      
	$trackid++; $lengthid++; $thetaid++; $recid++;
  
	#print "$trackid\t$focal\n";

	if($lengthid eq 2){$simchr=$focal; print "simulated chromosome length is $focal\n"}
	if($focal eq '-t'){$thetaid= -1;}
	if(($focal eq '-I') & ($parental_drift eq 0)){$trackid= -2; $pop_tracker=1;}
	if(($focal eq '-I') & ($parental_drift eq 1) & ($par1_drift >0) & ($par2_drift eq 0)){$trackid= -3; $pop_tracker=1;} #covers cases with par1 drift parameter
	#this is not quite right as it won't deal properly with sampling individuals from population 1
	if(($focal eq '-I') & ($parental_drift eq 1) & ($par1_drift eq 0) & ($par2_drift > 0)){$trackid= -4; $pop_tracker=1;} #covers cases with par2 drift parameter 
	if(($focal eq '-I') & ($parental_drift eq 1) & ($par1_drift>0) & ($par2_drift>0)){$trackid= -4; $pop_tracker=1;} #covers cases with both	
	if($focal eq '-R'){$recid= -1;}
	
	if($thetaid eq 0){
	    $simtheta=$focal; print "simulated theta is $simtheta\n";
	}#sim theta set
	if($recid eq 0){
	    $recombination_map=$focal; chomp $recombination_map; print "Recombination map used in macs and ancestry tract simulations is $recombination_map\n";
	}#store recombination map file for later use

	if(($focal ne '-I') & ($focal_prev ne '-I') & ($pop_tracker eq 1) & ($trackid le 1)){
	    $total_parental_haps=$total_parental_haps+$focal;
	    print "total current parental haplotype count is $total_parental_haps\n";
	    push(@track_pop,$focal);
	}#count total parental haplotypes

	$focal_prev=$focal;

    }#store parameters and par1 and par2 sample sizes

    #sort through stored par1 and par2 sample sizes to correctly set total numbers for par1 and par
    for my $popinfo (0..scalar(@track_pop)-1){
	$focalinfo=$track_pop[$popinfo];

	print "$popinfo\t$focalinfo\t$total_par1\t$total_par2\n";

	if(($popinfo eq 0) & ($par1_drift eq 0)){$total_par1=$focalinfo; print "simulating $total_par1 haplotypes from parent1\n";}
	if(($popinfo eq 1) & ($par2_drift eq 0)){$total_par2=$focalinfo; print "simulating $total_par2 haplotypes from parent2\n";}

	if(($popinfo eq 2) & ($par1_drift > 0)){$total_par1=$focalinfo; print "simulating $total_par1 haplotypes from parent1\n";}
	if(($popinfo eq 3) & ($par2_drift > 0)){$total_par2=$focalinfo; print "simulating $total_par2 haplotypes from parent2\n";}
    }#cycle through pop info and assign total haplotype numbers based on whether there is drift for that parent

    if($par1_aim>$total_par1){die "cannot run program with more sampled par1 haplotypes than simulated haplotypes ($par1_aim vs $total_par1)\nPlease check your configuration files\n";}
    if($par2_aim>$total_par2){die "cannot run program with more sampled par2 haplotypes than simulated haplotypes ($par2_aim vs $total_par2)\nPlease check your configuration files\n";}

    #generate sequence from trees
    system("macs $macs_params -T | msformatter | grep \'\\[\' > macs_simulation_results_trees.txt");

    print "macs $macs_params -T | msformatter | grep \'\\[\' > macs_simulation_results_trees.txt\n";

    my $partitions=qx(wc -l macs_simulation_results_trees.txt | perl -p -e 's/ +/\t/g' | cut -f 1); chomp $partitions;
    
    if($partitions eq 0){
	die "macs run was unsuccessful, check your macs parameters command and try again\n";
    }# do not procede if the macs command did not work

    if($seq_params !~ /a/g){
	$seq_params="-a 0.597 -f 0.305566 0.194762 0.194557 0.305115";
    }#set base composition etc

    if($use_ancestral ne 1){
    system("seq-gen -mHKY -l $simchr -s $simtheta -p $partitions -q $seq_params <macs_simulation_results_trees.txt > macs_simulation_results_trees.phy");
    }else{
	my $ancestral_raw=qx(cat $genome1 | tail -n +2); chomp $ancestral_raw;
	$ancestral_raw=~ s/^\s+//; $ancestral_raw=~ s/[\n\r]//g; #remove spaces, new lines
	$ancestral_raw=~ s/[NRMWSYK]/join "", map { ("A","T","C","G")[rand(4)] } 1 .. 1/eg;
 
	my $count = $ancestral_raw =~ s/(.)/$1/sg;
	open ANC, ">ancestral_seq_for_seq-gen";
	print ANC "1\t$count\nanc\    $ancestral_raw\n1\n";
        my $trees=qx(cat macs_simulation_results_trees.txt); chomp $trees;
	print ANC "$trees\n";
	system("seq-gen -mHKY -l $count -s $simtheta -p $partitions -q $seq_params -k 1 <ancestral_seq_for_seq-gen > macs_simulation_results_trees.phy");
    }#use ancestral sequence or don't

    print "seq-gen -mHKY -l $simchr -s $simtheta -p $partitions -q $seq_params <macs_simulation_results_trees.txt > macs_simulation_results_trees.phy\n";
   
    my $seq_params=qx(wc -l macs_simulation_results_trees.phy | perl -p -e 's/ +/\t/g' | cut -f 1); chomp $seq_params;

    if($seq_params eq 0){
	die "seq-gen sequence generation failed, check seq-gen or memory parameters\n";
    }#do not proceed in seq-gen didn't work

#####################################
####identify aims passing threshold:
#####################################

    my $select1=$total_par1;
    my $select2=$total_par1-1;
    my $max_par2=$select2+$total_par2;
    my $select_floor=-1;
    my $par1_index=$total_parental_haps-$total_par1-$total_par2;
    my $par2_index=$total_parental_haps-$total_par2;

    #modify here to deal with drift scenario
    my $par1_start=0; my $par2_start=$total_par1-1;
    if($parental_drift >0){
	if(($par1_drift>0) & ($par2_drift>0)){
	my $index1=$track_pop[0]+$track_pop[1];
	my $index2=$track_pop[0]+$track_pop[1]+$track_pop[2];
	$select_floor=$index1;
	$select1=$index2-1;
	$select2=$index2;
	$max_par2=$select2+$total_par2-1;
	print "selecting parent 1 haplotypes from haplotype $select_floor to $select1 and parent 2 haplotypes from haplotype $select2 to $max_par2\n";

	$par1_index=$select_floor;
	$par2_index=$select2;

	}elsif($par1_drift>0){
	    $select2=$total_par1;
	    $max_par2=$select2+$total_par2-1; #parent 2 is not modified
	    
	    my $index1=$track_pop[0]+$track_pop[1];
	    my $index2=$track_pop[0]+$track_pop[1]+$track_pop[2];
	    $select_floor=$index1;
	    $select1=$index2-1; #parent 1 is modified
	    print "selecting parent 1 haplotypes from haplotype $select_floor to $select1 and parent 2 haplotypes from haplotype $select2 to $max_par2\n";

	    $par1_index=$select_floor;
	    $par2_index=$select2;

	}elsif($par2_drift>0){
	    $select1=$total_par1-1;
	    $select_floor=0; #parent 1 is unmodified
	    
	    my $index2=$track_pop[0]+$track_pop[1]+$track_pop[2];
	    $select2=$index2;
	    $max_par2=$select2+$total_par2-1; #parent 2 is modified
	    print "selecting parent 1 haplotypes from haplotype $select_floor to $select1 and parent 2 haplotypes from haplotype $select2 to $max_par2\n";

	    $par1_index=$select_floor;
	    $par2_index=$select2;

	}#assign appropriate coordinates based on drift populations
    }#re-do select to accomodate this

    #Define & print reference genomes
    print "Parent 1 reference genome comes from individual $par1_index\n";
    my $par1_ref=qx(grep -w $par1_index macs_simulation_results_trees.phy); chomp $par1_ref;
    print "Parent 2 reference genome comes from individual $par2_index\n";
    my $par2_ref=qx(grep -w $par2_index macs_simulation_results_trees.phy); chomp $par2_ref;

    $par1_ref=~ s/ +/\t/g;
    $par2_ref=~ s/ +/\t/g;

    my @par1elements=split(/\t/,$par1_ref);
    my @par2elements=split(/\t/,$par2_ref);

    open SIM1, ">macs_simulated_parent1.fa";
    open SIM2, ">macs_simulated_parent2.fa";

    if(-e 'macs_simulated_parent1.fa.fai'){
        system("rm macs_simulated_parent1.fa.fai macs_simulated_parent2.fa.fai");
    }#remove previous files                                                                                                                                                   
    print SIM1 ">chr1\n";
    print SIM2 ">chr1\n";

    print SIM1 "$par1elements[1]\n";
    print SIM2 "$par2elements[1]\n";

    $genome1="macs_simulated_parent1.fa";
    $genome2="macs_simulated_parent2.fa";

    $chr="chr1";

    #subset parent 1 and parent2 haplotypes
    system("cat macs_simulation_results_trees.phy | tail -n +2 | awk \'\$1 <= $select1\' | awk \'\$1 >= $select_floor\' | shuf -n $par1_aim > macs_simulation_results_trees.par1.phy");
    system("cat macs_simulation_results_trees.phy | awk \'\$1 >= $select2\' | awk \'\$1 <= $max_par2\' | shuf -n $par2_aim > macs_simulation_results_trees.par2.phy");
    
    print "par1 select command is:","awk \'\$1 <= $select1\' macs_simulation_results_trees.phy | awk \'\$1 >= $select_floor\' | shuf -n $par1_aim > macs_simulation_results_trees.par1.phy","\n";
    print "par2 select command is:","awk \'\$1 >= $select2\' macs_simulation_results_trees.phy | awk \'\$1 <= $max_par2\' | tail -n +2 | shuf -n $par2_aim > macs_simulation_results_trees.par2.phy","\n";

    #convert to fasta format for downstream
    system("perl $program_path/Phylip2Fasta.pl macs_simulation_results_trees.par1.phy macs_simulation_results_trees.par1.fa");
    system("perl $program_path/Phylip2Fasta.pl macs_simulation_results_trees.par2.phy macs_simulation_results_trees.par2.fa");

    #identify variant sites
    system("python $program_path/Variable_sites_extractor_mod.py macs_simulation_results_trees.par1.fa -c -v -o macs_simulation_results_trees.par1.variable");
    system("python $program_path/Variable_sites_extractor_mod.py macs_simulation_results_trees.par2.fa -c -v -o macs_simulation_results_trees.par2.variable");

    system("cat macs_simulation_results_trees.par1.variable | perl -p -e 's/,/\n/g' | tail -n +2 > macs_simulation_results_trees.par1.variable.coordinates");
    system("cat macs_simulation_results_trees.par2.variable | perl -p -e 's/,/\n/g' | tail -n +2 > macs_simulation_results_trees.par2.variable.coordinates");
    system("cat macs_simulation_results_trees.par1.variable.coordinates macs_simulation_results_trees.par2.variable.coordinates | uniq | sort -nk 1 > macs_simulation_results_trees.bothpar.variable.coordinates");
    system("perl $program_path/identify_AIMs_sims.pl macs_simulated_parent1.fa macs_simulated_parent2.fa | cut -f 2 > macs_simulation_results_trees.aims");

    open SHAREDPOLY, ">macs_simulation_results_trees.bothpar.variable.coordinates.cor";
    open CURR, "macs_simulation_results_trees.bothpar.variable.coordinates" or die "cannot open macs_simulation_results_trees.bothpar.variable.coordinates file\n";

    while(my $currentpoly=<CURR>){
	chomp $currentpoly;
	my $corrected=$currentpoly+1;
	print SHAREDPOLY "$corrected\n";
    }#correct coords

    #FAS scriptome snippet
    system("perl $program_path/fas_merge_by_column_snippet.pl macs_simulation_results_trees.bothpar.variable.coordinates.cor macs_simulation_results_trees.aims");

    system("perl $program_path/seq-gen_sequences_to_aim_definition.pl macs_simulation_results_trees.aims.sharedpoly macs_simulation_results_trees.par1.fa macs_simulation_results_trees.par2.fa $freq_diff");

    print "perl $program_path/seq-gen_sequences_to_aim_definition.pl macs_simulation_results_trees.aims.sharedpoly  macs_simulation_results_trees.par1.fa macs_simulation_results_trees.par2.fa $freq_diff\n";

    #mask the parental genomes
    system("seqtk mutfa macs_simulated_parent1.fa macs_simulation_results_trees.mask.insnp > macs_simulation.par1.masked.fa");
    system("seqtk mutfa macs_simulated_parent2.fa macs_simulation_results_trees.mask.insnp > macs_simulation.par2.masked.fa");

    #cleanup here
    system("mv macs_simulation.par1.masked.fa macs_simulated_parent1.fa");
    system("mv macs_simulation.par2.masked.fa macs_simulated_parent2.fa");

   system("rm macs_simulation_results_trees.par1.variable.coordinates macs_simulation_results_trees.par2.variable.coordinates macs_simulation_results_trees.mask.insnp macs_simulation_results_trees.bothpar.variable.coordinates macs_simulation_results_trees.bothpar.variable.coordinates.cor macs_simulation_results_trees.aims");

}#macs command & sequence generation

#catalog AIMs, extract focal chromosome                                                                                            
my $chr1="$chr"."_select_par1.fa";
my $chr2="$chr"."_select_par2.fa";
system("perl $program_path/getScaffold_samtools.pl $genome1 $chr > $chr1");
system("perl $program_path/getScaffold_samtools.pl $genome2 $chr > $chr2");

my $chr_length=qx(cat $chr1 | tail -n +2 | perl -p -e 's/\n//g' | wc -c | perl -p -e 's/ +/\t/g' | cut -f 1); chomp $chr_length;

print "number of basepairs in $chr is $chr_length\n";

#define genomes and AIMs
my $aims="simulation_ancestry_informative_sites_"."$genome1"."_"."$genome2"."_"."$chr";
system("perl $program_path/identify_AIMs_sims.pl $chr1 $chr2 > $aims");

if($macs eq 0){
    system("perl $program_path/identify_AIMs_sims.pl $chr1 $chr2 > $aims");
    #generate frequency file based on rate of shared polymorphism if not using macs
    system("Rscript $program_path/poly_dist.R $rate_shared_poly $poly_par1 $poly_par2 $chr_length $num_indivs $aims $error");
}#determine aims differently depending on whether you are using macs or not
if($macs eq 1){
    #new script here 
    system("perl $program_path/overlap_AIMs_and_counts_sim.pl simulated_parental_counts_for_AncestryHMM_sharedpoly $aims $rec_rate_Morgans $total_par1 $total_par2 $read_length");
}#add back in fixed AIMs for macs simulation

#run SELAM to generate ancestry tracts for later use
my $length_morgans= ($rec_rate_Morgans)*($chr_length/1000);
print "length of the chromosome in morgans is $length_morgans\n";
system("SELAM -d selam_demography.txt -o selam_simulation_output_parameters.txt -c 2 $length_morgans 0");

print "SELAM command is: SELAM -d selam_demography.txt -o selam_simulation_output_parameters.txt -c 2 $length_morgans 0\n";

open LIST, ">simulated_hybrids_readlist_"."gen"."$gens_since_admix"."_prop_par1_"."$mixture_proportion";
my $reads_folder="simulated_hybrids_reads_"."gen"."$gens_since_admix"."_prop_par1_"."$mixture_proportion";

if(-e $reads_folder){
    print "WARNING: removing previous results in $reads_folder\n";
    system("rm -r $reads_folder");
}#remove previous results
system("mkdir $reads_folder");

#define files and counters
my $current_outfile="split_file_list_1"; my $current_slurm="slurm_batch1.sh"; my $counter=0; my $track=0; my $string="";

if(length($recombination_map) eq 0){
    $recombination_map=0;
}#if map is absent, set to zero as script placeholder

while($counter<$num_indivs){
    $counter=$counter+1; $track=$track+1;
    #print "$counter\t$track\n";
    
    if(($track>=$num_indiv_per_job) or ($counter eq $num_indivs) or ($counter eq 1)){

	if($track eq $num_indiv_per_job){
	
	    #print commands for current job
	    if($macs eq 1){
	    print SLURM "perl $program_path/generate_genomes_and_reads_v2.pl $genome1 $genome2 $chr $chr1 $chr2 $chr_length $poly_par1 $poly_par2 $aims $error $reads_folder $mixture_proportion $gens_since_admix $rec_rate_Morgans $snp_freqs $current_outfile $number_reads $read_length $sequence_error $per_bp_indel $macs $total_par1 $total_par2 $use_map $rec_rate_Morgans $recombination_map $program_path\n";

	    } elsif($macs eq 0){
	    print SLURM "perl $program_path/generate_genomes_and_reads_v2.pl $genome1 $genome2 $chr $chr1 $chr2 $chr_length $poly_par1 $poly_par2 $aims $error $reads_folder $mixture_proportion $gens_since_admix $rec_rate_Morgans $snp_freqs $current_outfile $number_reads $read_length $sequence_error $per_bp_indel $macs 1 1 $use_map $rec_rate_Morgans $recombination_map $program_path\n";
	    }#submit current job
	    my $jobid=qx($job_submit_cmd $current_slurm); chomp $jobid;
	    my @jobarray=split(/ /,$jobid);
	    #print "$jobarray[-1]\n";
	    if(length($string) eq 0){
		$string="$jobarray[-1]";
	    } else{
	    $string="$string".","."$jobarray[-1]";
	    }#final cleanup dependencies
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

open CLEANUP, ">cleanup.sh";
for my $k (1..scalar(@commands_array)-1){
    print CLEANUP "#"."$commands_array[$k]"."\n";
}#print out all slurm header elements

print CLEANUP "rm admix_simulation_demography_output_results*"."\n";
print CLEANUP "rm slurm_batch*"."\n";
print CLEANUP "rm split_file_list_*"."\n";
print CLEANUP "rm macs_simulation_results_trees*"."\n";
print CLEANUP "rm par*_coordinates_fastahack"."\n";
print CLEANUP "rm *fai"."\n";
print CLEANUP "rm indiv*_log"."\n";

if($use_map eq 1){
    print CLEANUP "perl $program_path/prior_recombination_map_ancestryHMM.pl simulated_parental_counts_for_AncestryHMM $recombination_map $rec_rate_Morgans $chr_length $program_path\n";
}#correct map if needed

#print "cleaning up after: $string\n";

#!system("sbatch --dependency=afterok:$string cleanup.sh");
