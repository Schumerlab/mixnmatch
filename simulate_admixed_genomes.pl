#perl! -w

my $infile=shift(@ARGV); chomp $infile;
open IN, $infile or die "cannot open configuration file\n";

my $genome1=""; my $genome2=""; my $mixture_proportion=""; my $rec_rate_Morgans=""; my $num_indivs=""; my $gens_since_admix=""; my $rate_shared_poly=""; my $per_bp_indel=""; my $number_reads=""; my $gens_drift_parental=""; my $sequence_error=""; my $read_type=""; my $read_length=""; my $poly_par1=""; my $poly_par2=""; my $chr="";

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
    if($line =~ /gens_drift_parental=/g){
	$gens_drift_parental=$elements[1]; chomp $gens_drift_parental;
    }#generations drift
    if($line =~ /sequencing_error=/g){
	$sequence_error=$elements[1]; chomp $sequence_error;
    }#error rate

}#for all lines in the infile

my $counter=0; my $track=0;
my $snp_freqs="shared_polymorphism_distribution";
my $start=1;
my $stop=1;
my $freq_diff=0.99;
my $error=1-$freq_diff;

#catalog AIMs, extract focal chromosome
my $chr1="$chr"."_select_par1.fa";
my $chr2="$chr"."_select_par2.fa";
system("perl getScaffold_samtools.pl $genome1 $chr > $chr1");
system("perl getScaffold_samtools.pl $genome2 $chr > $chr2");

my $chr_length=qx(cat $chr1 | tail -n +2 | perl -p -e 's/\n//g' | wc -c | perl -p -e 's/ +/\t/g' | cut -f 1); chomp $chr_length;

print "number of basepairs in $chr is $chr_length\n";

my $aims="simulation_ancestry_informative_sites_"."$genome1"."_"."$genome2"."_"."$chr";
system("perl identify_AIMs_two_genomes.pl $chr1 $chr2 > $aims");

#generate frequency file based on rate of shared polymorphism                                                                                     
print "Rscript poly_dist.R $rate_shared_poly $poly_par1 $poly_par2 $chr_length $num_indivs $aims $error\n";
system("Rscript poly_dist.R $rate_shared_poly $poly_par1 $poly_par2 $chr_length $num_indivs $aims $error");

open LIST, ">simulated_hybrids_readlist_"."gen"."$gens_since_admix"."_prop_par1_"."$mixture_proportion";
while($counter<=$num_indivs){
    $counter=$counter+1; $track=$track+1;
    
    my $chr1="$chr"."_select_par1.fa";
    my $chr2="$chr"."_select_par2.fa";

    #print "perl generate_shared_private_poly.pl $counter $aims $snp_freqs 0.001 0.0002 $chr_length $chr $chr1 $chr2\n";
    system("perl generate_shared_private_poly.pl $counter $aims $snp_freqs 0.001 0.0002 $chr_length $chr $chr1 $chr2");

    $chr1="indiv"."$counter"."_hap1.fa";
    $chr2="indiv"."$counter"."_hap2.fa";

    #generate and stitch together tracts
    my $chr1_tracts=qx(Rscript tract_lengths.R $mixture_proportion $gens_since_admix $rec_rate_Morgans $chr_length);
    #print "Rscript tract_lengths.R $mixture_proportion $gens_since_admix $rec_rate_Morgans $chr_length\n";
    my @chrtracts1=split(/\n/,$chr1_tracts);
    my $chr2_tracts=qx(Rscript tract_lengths.R $mixture_proportion $gens_since_admix $rec_rate_Morgans $chr_length);
    my @chrtracts2=split(/\n/,$chr2_tracts);

    my $bed1="indiv"."$counter"."_tracts_hap1.bed"; my $bed2="indiv"."$counter"."_tracts_hap2.bed";
    open BED1, ">$bed1";
    open BED2, ">$bed2";

    my $start_bed=1; my $stop_bed=1;
    for my $m (0..scalar(@chrtracts1)-1){
	if ($m % 2 == 0){
	   $stop_bed=$chrtracts1[$m]+$start_bed-1;
	    print BED1 "$chr\t$start_bed\t$stop_bed\tpar1\n";
	} else{
	   $stop_bed=$chrtracts1[$m]+$start_bed-1;
	    print BED1 "$chr\t$start_bed\t$stop_bed\tpar2\n";
      	}#check odd or even for ancestry purposes
	$start_bed=$stop_bed+1;
    }#for all tracts for haplotype 1

    my $start_bed=1; my $stop_bed=1;
    for my $n (0..scalar(@chrtracts2)-1){
        if ($n % 2 == 0){
	    $stop_bed=$chrtracts2[$n]+$start_bed-1;
            print BED2 "$chr\t$start_bed\t$stop_bed\tpar1\n";
        } else{
	    $stop_bed=$chrtracts2[$n]+$start_bed-1;
            print BED2 "$chr\t$start_bed\t$stop_bed\tpar2\n";
        }#check odd or even for ancestry purposes                                                                                
        $start_bed=$stop_bed+1;
    }#for all tracts for haplotype 2   

    my $seq1=""; my $seq2=""; my $total_length_chr1=0; my $total_length_chr2=0;
    open OUT, ">"."indiv"."$counter".".fa";
    for my $j (0..scalar(@chrtracts1)-1){
	my $current_tract=$chrtracts1[$j];

       	if($j eq 0){$start=1;}

	$stop=$current_tract+$total_length_chr1;

	if (($j % 2 == 0) && ($current_tract != 0)) {
	    my $par2 = qx(fastahack $chr2 -r $chr:$start..$stop); chomp $par2; 
	    $seq1="$seq1"."$par2";
	} elsif($current_tract != 0){
	    my $par1 = qx(fastahack $chr1 -r $chr:$start..$stop); chomp $par1;
	    $seq1="$seq1"."$par1";
	}#odd or even?

	$start=$start+$current_tract;
        $total_length_chr1=$total_length_chr1+$current_tract;

   }#generate chromosome 1 sequence

    print OUT ">indiv"."$j"."_hap1\n"."$seq1\n";

    for my $j (0..scalar(@chrtracts2)-1){
        my $current_tract=$chrtracts2[$j];

        if($j eq 0){$start=1;}

	$stop=$current_tract+$total_length_chr2;
	#print "$start\t$stop\n";
	if (($j % 2 == 0) && ($current_tract != 0)) {
            my $par2 = qx(fastahack $chr2 -r $chr:$start..$stop); chomp $par2;
            $seq2="$seq2"."$par2";
        } elsif($current_tract != 0){
            my $par1 = qx(fastahack $chr1 -r $chr:$start..$stop); chomp $par1;
            $seq2="$seq2"."$par1";
        }#odd or even?                                                                                                                                    
	$start=$start+$current_tract;
        $total_length_chr2=$total_length_chr2+$current_tract;

    }#generate chromosome 2 sequence 

    print OUT ">indiv"."$j"."_hap2\n"."$seq2\n";

    #generate indels
    #generate reads

	my $name="indiv"."$counter".".fa";
	my $r1="indiv"."$counter"."_read1.fq"; my $r2="indiv"."$counter"."_read2.fq";
       
	system("wgsim -N $number_reads -1 $read_length -2 $read_length -S $counter -e $sequence_error -r $per_bp_indel -R 1 $name $r1 $r2");

    print LIST "$r1\t$r2\n";

    #cleanup
    my $cleanup="indiv"."$counter"."_hap"."*";
    system("rm $cleanup");

}#for all jobs


