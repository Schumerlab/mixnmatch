#perl! -w

if(@ARGV<20){
print "error in command line options for genome generation"; exit;
}#print usage, exit
         
my $genome1=shift(@ARGV); chomp $genome1;
my $genome2=shift(@ARGV); chomp $genome2;
my $chr=shift(@ARGV); chomp $chr;
my $chr1=shift(@ARGV); chomp $chr1;
my $chr2=shift(@ARGV); chomp $chr2;
my $chr_length=shift(@ARGV); chomp $chr_length;
my $poly_par1=shift(@ARGV); chomp $poly_par1;
my $poly_par2=shift(@ARGV); chomp $poly_par2;
my $aims=shift(@ARGV); chomp $aims;
my $error=shift(@ARGV); chomp $error;
my $reads_folder=shift(@ARGV); chomp $reads_folder;
my $mixture_proportion=shift(@ARGV); chomp $mixture_proportion;
my $gens_since_admix=shift(@ARGV); chomp $gens_since_admix;
my $rec_rate_Morgans=shift(@ARGV); chomp $rec_rate_Morgans;
my $snp_freqs=shift(@ARGV); chomp $snp_freqs;
my $id_list=shift(@ARGV); chomp $id_list;
my $number_reads=shift(@ARGV); chomp $number_reads;
my $read_length=shift(@ARGV); chomp $read_length;
my $sequence_error=shift(@ARGV); $sequence_error;
my $per_bp_indel=shift(@ARGV); chomp $per_bp_indel;
my $use_macs=shift(@ARGV); chomp $use_macs;
my $total_par1=shift(@ARGV); chomp $total_par1;
my $total_par2=shift(@ARGV); chomp $total_par2;
my $use_map=shift(@ARGV); chomp $use_map;
my $rec_rate_Morgans=shift(@ARGV); chomp $rec_rate_Morgans;
my $recombination_map=shift(@ARGV); chomp $recombination_map;
my $cross_contam=shift(@ARGV); chomp $cross_contam;
my $program_path=shift(@ARGV); chomp $program_path;

open IN, $id_list or die "can't open list of indivs\n";

while(my $line=<IN>){
    chomp $line;
    $counter=$line;
    
    my $chr1="$chr"."_select_par1.fa";
    my $chr2="$chr"."_select_par2.fa";

    if($use_macs eq 0){

    system("perl $program_path/generate_shared_private_poly.pl $counter $aims $snp_freqs $poly_par1 $poly_par2 $chr_length $chr $chr1 $chr2");
    
    my $current_haps="admix_simulation_demography_output_results_"."$counter";

    #set the haplotypes to the just-generated parental haplotypes
    my $hap1="indiv"."$counter"."_"."hap1.fa";
    my $hap2="indiv"."$counter"."_"."hap2.fa";

    #generate and stitch together tracts
    system("perl $program_path/selam_to_ancestry_tracts.pl $current_haps 0 $hap1 $hap2 $counter $total_par1 $total_par2 $use_map $rec_rate_Morgans $recombination_map $program_path");
    
    print "perl $program_path/selam_to_ancestry_tracts.pl $current_haps 0 $hap1 $hap2 $counter $total_par1 $total_par2 $use_map $rec_rate_Morgans $recombination_map $program_path\n";

    }#only generate polymorphic versions if not using macs
    elsif($use_macs eq 1){
	my $current_haps="admix_simulation_demography_output_results_"."$counter";

    #generate and stitch together tracts                                                                                                             
    system("perl $program_path/selam_to_ancestry_tracts.pl $current_haps 1 $genome1 $genome2 $counter $total_par1 $total_par2 $use_map $rec_rate_Morgans $recombination_map $program_path");
    
    print "perl $program_path/selam_to_ancestry_tracts.pl $current_haps 0 $genome1 $genome2 $counter $total_par1 $total_par2 $use_map $rec_rate_Morgans $recombination_map$program_path\n";

    }#otherwise use macs
 
    #generate indels
    #generate reads

	my $name="indiv"."$counter".".fa";
	my $r1="$reads_folder"."/"."indiv"."$counter"."_read1.fq"; my $r2="$reads_folder"."/"."indiv"."$counter"."_read2.fq";
       
    system("wgsim -N $number_reads -1 $read_length -2 $read_length -S $counter -e $sequence_error -r $per_bp_indel -R 1 $name $r1 $r2");

    if($cross_contam > 0){

	my $number_contam_par1=$mixture_proportion*$number_reads*$cross_contam;
	my $number_contam_par2=(1-$mixture_proportion)*$number_reads*$cross_contam;
	my $contam_seed1=$counter+1; my $contam1_par1=""; my $contam2_par1="";
	my $contam_seed2=$counter+2; my $contam1_par2=""; my $contam2_par2="";
	print "generating $number_contam_par1 contamination reads for parent1\n";
	print "generating $number_contam_par2 contamination reads for parent2\n";

	if($use_macs eq 0){
	    $contam1_par1="$reads_folder"."/"."indiv"."$counter"."_par1contam_read1.fq"; $contam2_par1="$reads_folder"."/"."indiv"."$counter"."_par1contam_read2.fq";
	    $contam1_par2="$reads_folder"."/"."indiv"."$counter"."_par2contam_read1.fq"; $contam2_par2="$reads_folder"."/"."indiv"."$counter"."_par2contam_read2.fq";
	    system("wgsim -N $number_contam_par1 -1 $read_length -2 $read_length -S $contam_seed1 -e $sequence_error -r $per_bp_indel -R 1 $chr1 $contam1_par1 $contam2_par1");
	    system("wgsim -N $number_contam_par2 -1 $read_length -2 $read_length -S $contam_seed2 -e $sequence_error -r $per_bp_indel -R 1 $chr2 $contam1_par2 $contam2_par2"); 
	}else{
	    $contam1_par1="$reads_folder"."/"."indiv"."$counter"."_par1contam_read1.fq"; $contam2_par1="$reads_folder"."/"."indiv"."$counter"."_par1contam_read2.fq";
            $contam1_par2="$reads_folder"."/"."indiv"."$counter"."_par2contam_read1.fq"; $contam2_par2="$reads_folder"."/"."indiv"."$counter"."_par2contam_read2.fq";
            system("wgsim -N $number_contam_par1 -1 $read_length -2 $read_length -S $contam_seed1 -e $sequence_error -r $per_bp_indel -R 1 macs_simulation_results_trees.par1.fa $contam1_par1 $contam2_par1");
            system("wgsim -N $number_contam_par2 -1 $read_length -2 $read_length -S $contam_seed2 -e $sequence_error -r $per_bp_indel -R 1 macs_simulation_results_trees.par2.fa $contam1_par2 $contam2_par2");
	}#generate appropriate contamination source

	#combine contamination and non-contamination reads

	my $z1="$r1".".gz"; my $z2="$r2".".gz";
       
	system("cat $r1 $contam1_par1 $contam1_par2 | gzip > $z1");
	system("cat $r2 $contam2_par1 $contam2_par2 | gzip > $z2");
	system("mv $name $reads_folder");

	system("rm $r1 $contam1_par1 $contam2_par1 $r2 $contam1_par2 $contam2_par2");

    }#add cross contamination
    if($cross_contam eq 0){

	system("gzip $r1 $r2");
	system("mv $name $reads_folder");

    }#otherwise just zip those files

    print LIST "$r1".".gz"."\t$r2".".gz"."\n";
    
    my $combined_bed="indiv"."$counter".".bed";
    my $bed1="indiv"."$counter"."_hap1.bed";
    my $bed2="indiv"."$counter"."_hap2.bed";

    system("bedtools intersect -a $bed1 -b $bed2 -wo | sort -n -k2,2 -k7,7  > $combined_bed");
    system("mv $combined_bed $reads_folder");

    #cleanup
    my $cleanup="indiv"."$counter"."_hap"."*";
    system("rm $cleanup");
}
