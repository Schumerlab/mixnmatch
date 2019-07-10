#perl! -w 

use Math::Round;

if(@ARGV<4){
    print "perl seq-gen_sequences_to_aim_definition.pl shared_polymorphism par1_seqs.fa par2_seqs.fa diff_thresh\n"; exit;
}#print usage

my $shared_poly=shift(@ARGV); chomp $shared_poly;
open IN, $shared_poly or die "cannot open shared polymorphism file\n";

my $par1=shift(@ARGV); chomp $par1;

my $par2=shift(@ARGV); chomp $par2;

my $thresh=shift(@ARGV); chomp $thresh;

my $par1num=qx(grep '>' $par1 | wc -l | perl -p -e 's/ +/\t/g' | cut -f 1); chomp $par1num;
my $par1seqids=qx(grep '>' $par1 | perl -p -e 's/>//g'); chomp $par1seqids;
my @par1seqs=split(/\n/,$par1seqids);

my $par2num=qx(grep '>' $par2 | wc -l | perl -p -e 's/ +/\t/g' | cut -f 1); chomp $par2num;
my $par2seqids=qx(grep '>' $par2 | perl -p -e 's/>//g'); chomp $par2seqids;
my @par2seqs=split(/\n/,$par2seqids);

my $total=$par1num + $par2num;

$par1num=$par1num-1;
$par2num=$par2num;

open INSNP, ">macs_simulation_results_trees.mask.insnp";
open OUT, ">simulated_parental_counts_for_AncestryHMM_sharedpoly";
open AIMs, ">simulated_AIMs_for_AncestryHMM_sharedpoly";
#print "$par1num\t$par2num\n";

while(my $line=<IN>){

    chomp $line;
    my $coord=$line;

    #print "$coord\n";
    open OUT1, ">par1_coordinates_fastahack";
    open OUT2, ">par2_coordinates_fastahack";

    for my $i (0..scalar(@par1seqs)-1){
	my $currseq=$par1seqs[$i];
	print OUT1 "$currseq:$coord\n";
	#print "$i:$line\n";
    }#go through all par1 sequences
    for my $j (0..scalar(@par1seqs)-1){
	my $currseq=$par2seqs[$j];
	print OUT2 "$currseq:$coord\n";
	#print "$j:$line\n";
    }#go through all par2 sequences

    my @combined=();

    my $par1seqs=qx(cat par1_coordinates_fastahack | fastahack macs_simulation_results_trees.par1.fa -c); chomp $par1seqs;
    #print "$par1seqs\n";
    my @par1bps=split(/\n/,$par1seqs);

    my $par2seqs=qx(cat par2_coordinates_fastahack | fastahack macs_simulation_results_trees.par2.fa -c); chomp $par2seqs;
    #print "$par2seqs\n\n";
    my @par2bps=split(/\n/,$par2seqs);

    push(@combined, @par1bps); push(@combined,@par2bps);

    my @uniquepar= uniq(@combined);
    #print scalar(@uniquepar1),"\n";

    my $num_alleles=scalar(@uniquepar);

    my $par1freq=0; my $par2freq=0; my $par1count=0; my $par2count=0; my $a1=""; my $a2="";
    if($num_alleles eq 2){

	$a1=$uniquepar[0];
	$a2=$uniquepar[1];

    for my $k (0..scalar(@par1bps)-1){

	if($par1bps[$k] eq $a1){
	    $par1count++;
	}#count a1

    }#count all bases

    for my $l (0..scalar(@par2bps)-1){

	if($par2bps[$l] eq $a1){
	    $par2count++;
	}#count a1

    }#count all bases

	if(scalar(@par1bps) >0 & scalar(@par2bps) >0){
	$par1freq=$par1count/scalar(@par1bps);
	$par2freq=$par2count/scalar(@par2bps);
	}#count first
	#print "$par1freq\t$par2freq\n";

 }#only evaluate biallelic sites

    my $thresh_inverse=1-$thresh;
    if((($par1freq < $thresh) && ($par1freq > $thresh_inverse)) or (($par2freq < $thresh) &&($par2freq > $thresh_inverse))){
	print INSNP "chr1\t$coord\tX\tN\n";
	#!print "MASK\t$par1freq\t$par2freq\t$coord\n";
    }#site does not pass due to frequency, mask
    elsif((($par1freq < $thresh) && ($par2freq < $thresh)) or (($par1freq < $thresh_inverse) && ($par2freq < $thresh_inverse))){
	print INSNP "chr1\t$coord\tX\tN\n";
	#!print "MASK\t$par1freq\t$par2freq\t$coord\n";
    }#site does not pass due to sharing mask
    elsif((($par1freq>= $thresh) && ($par2freq <= $thresh_inverse)) or (($par2freq>= $thresh) && ($par1freq <= $thresh_inverse))){
	my $a1_1=round($par1freq*$par1num*2); my $a2_1=round((1-$par1freq)*$par1num*2);
	my $a1_2=round($par2freq*$par2num*2); my $a2_2=round((1-$par2freq)*$par2num*2);
        print "$coord\t$par1freq\t$par2freq\n"; 
	print OUT "chr1\t$coord\t$a1_1\t$a2_1\t$a1_2\t$a2_2\n";
	print AIMs "chr1\t$coord\t$a1\t$a2\n";
	#need to take out the dummy rate variable
    }#passes all thresholds

}#process all shared polymorphisms


sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}
