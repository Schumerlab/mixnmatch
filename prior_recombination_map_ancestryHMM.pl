#perl! -w

if(@ARGV<3){
    print "perl prior_recombination_map_ancestryHMM.pl ancestry_hmm_file rec_map_macs rate_M_per_kb chrom_length program_path\n";
}#print usage and exit

my $infile=shift(@ARGV); chomp $infile;
open IN, $infile or die "cannot open $infile\n";

my $map=shift(@ARGV); chomp $map;

my $rate=shift(@ARGV); chomp $rate;

my $chrom=shift(@ARGV); chomp $chrom;

my $path=shift(@ARGV); chomp $path;

open OUT, ">$infile"."_recombination_prior.txt";

my $start_prev=1; my $counter=0; my $local_rate=0;
while(my $line=<IN>){

    $counter++;

    my @elements=split(/\t/,$line);
    my $stop=$elements[1];

    if($counter ne 1){
	#print "Rscript bp_to_morgans_recmap.R $map $start_prev $stop $rate $chrom\n";
    $local_rate=qx(Rscript $path/bp_to_morgans_recmap.R $map $start_prev $stop $rate $chrom); chomp $local_rate;
    }

    print OUT "$elements[0]\t$elements[1]\t$elements[2]\t$elements[3]\t$elements[4]\t$elements[5]\t$local_rate\n";

    $start_prev=$stop;

}#go through all lines
