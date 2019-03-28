#perl! -w

if(@ARGV<6){
print ""; exit;
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

open IN, $id_list or die "can't open list of indivs\n";

while(my $line=<IN>){
    chomp $line;
    $counter=$line;
    
    my $chr1="$chr"."_select_par1.fa";
    my $chr2="$chr"."_select_par2.fa";

    #print "perl generate_shared_private_poly.pl $counter $aims $snp_freqs 0.001 0.0002 $chr_length $chr $chr1 $chr2\n";
    system("perl generate_shared_private_poly.pl $counter $aims $snp_freqs 0.001 0.0002 $chr_length $chr $chr1 $chr2");

    my $chr1_mix="indiv"."$counter"."_hap1.fa";
    my $chr2_mix="indiv"."$counter"."_hap2.fa";

    #generate and stitch together tracts
    my $chr1_tracts=qx(Rscript tract_lengths.R $mixture_proportion $gens_since_admix $rec_rate_Morgans $chr_length);
    #print "Rscript tract_lengths.R $mixture_proportion $gens_since_admix $rec_rate_Morgans $chr_length\n";
    my @chrtracts1=split(/\n/,$chr1_tracts);
    my $chr2_tracts=qx(Rscript tract_lengths.R $mixture_proportion $gens_since_admix $rec_rate_Morgans $chr_length);
    my @chrtracts2=split(/\n/,$chr2_tracts);

    my $bed1="indiv"."$counter"."_tracts_hap1.bed"; my $bed2="indiv"."$counter"."_tracts_hap2.bed";
    open BED1, ">$bed1";
    open BED2, ">$bed2";

    my $seq1=""; my $seq2=""; my $par1=""; my $par2="";

    my $start_bed=1; my $stop_bed=1;
    for my $m (0..scalar(@chrtracts1)-1){
	if ($m % 2 == 0){
	   $stop_bed=$chrtracts1[$m]+$start_bed;
	   print BED1 "$chr\t$start_bed\t$stop_bed\tpar1\n";
	   $par1 = qx(fastahack $chr1_mix -r $chr:$start_bed..$stop_bed); chomp $par1;
	   #print "fastahack $chr1_mix -r $chr:$start_bed..$stop_bed\n";
	   $seq1="$seq1"."$par1";
	} else{
	   $stop_bed=$chrtracts1[$m]+$start_bed;
	   print BED1 "$chr\t$start_bed\t$stop_bed\tpar2\n";
	   $par2 = qx(fastahack $chr2_mix -r $chr:$start_bed..$stop_bed); chomp $par2;
	   #print "fastahack $chr2_mix -r $chr:$start_bed..$stop_bed\n"; 
	   $seq1="$seq1"."$par2";
      	}#check odd or even for ancestry purposes
	$start_bed=$stop_bed+1;
    }#for all tracts for haplotype 1

    my $start_bed=1; my $stop_bed=1;
    for my $n (0..scalar(@chrtracts2)-1){
        if ($n % 2 == 0){
	    $stop_bed=$chrtracts2[$n]+$start_bed;
            print BED2 "$chr\t$start_bed\t$stop_bed\tpar1\n";
	    $par1 = qx(fastahack $chr1_mix -r $chr:$start_bed..$stop_bed); chomp $par1;
	    #print "fastahack $chr1_mix -r $chr:$start_bed..$stop_bed\n";
	    $seq2="$seq2"."$par1";
        } else{
	    $stop_bed=$chrtracts2[$n]+$start_bed;
            print BED2 "$chr\t$start_bed\t$stop_bed\tpar2\n";
	    $par2 = qx(fastahack $chr2_mix -r $chr:$start_bed..$stop_bed); chomp $par2;
	    #print "fastahack $chr2_mix -r $chr:$start_bed..$stop_bed\n";
	    $seq2="$seq2"."$par2";
        }#check odd or even for ancestry purposes                                                                                
        $start_bed=$stop_bed+1;
    }#for all tracts for haplotype 2   

    my $combined_bed="indiv"."$counter"."_tracts.bed";
    system("bedtools intersect -a $bed1 -b $bed2 -wo > $combined_bed");
    system("rm $bed1 $bed2");

    open OUT, ">"."indiv"."$counter".".fa";
    print OUT ">indiv"."$j"."_hap1\n"."$seq1\n";
    print OUT ">indiv"."$q"."_hap2\n"."$seq2\n";

    #generate indels
    #generate reads

	my $name="indiv"."$counter".".fa";
	my $r1="$reads_folder"."/"."indiv"."$counter"."_read1.fq"; my $r2="$reads_folder"."/"."indiv"."$counter"."_read2.fq";
       
	system("wgsim -N $number_reads -1 $read_length -2 $read_length -S $counter -e $sequence_error -r $per_bp_indel -R 1 $name $r1 $r2");

    print LIST "$r1".".gz"."\t$r2".".gz"."\n";
    
    system("gzip $r1 $r2");
    system("mv $combined_bed $reads_folder");
    system("mv $name $reads_folder");

    #cleanup
    my $cleanup="indiv"."$counter"."_hap"."*";
    system("rm $cleanup");

}

