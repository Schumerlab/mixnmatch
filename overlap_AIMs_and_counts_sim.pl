#perl! -w

my $poly=shift(@ARGV); chomp $poly;

my $aims=shift(@ARGV); chomp $aims;
open IN, $aims or die "cannot open aims file\n";

my $rate=shift(@ARGV); chomp $rate;
my $par1=shift(@ARGV); chomp $par1;
my $par2=shift(@ARGV); chomp $par2;
my $read_length=shift(@ARGV); chomp $read_length;
$rate=$rate/1000;

open OUT, ">simulated_parental_counts_for_AncestryHMM";
open AIMS, ">simulated_AIMs_for_AncestryHMM";

my $currbp=""; my $bpprev=1; my $curr_rate="";
while(my $line=<IN>){

    chomp $line;

    my @elements=split(/\t/,$line); $currbp=$elements[1];

    $curr_rate=$rate*($currbp-$bpprev);

    my $distance=$currbp-$bpprev;

    open POLY, $poly or die "cannot open shared poly file\n"; 
    my $tracker=0;
    while (my $polyline=<POLY>){
	chomp $polyline; my @polybp=split(/\t/,$polyline); $polyfocal=$polybp[1]; 
	#print "$currbp\t$polyline\n";
	if($currbp eq $polyfocal){
	    if($distance >= $read_length){
	    print OUT "$polyline\t$curr_rate\n"; $tracker=1;
	    print AIMS "$line\n"; $bpprev=$currbp;
	    }#one snp per read
	} elsif(($polyfocal < $currbp) && ($polyfocal > $bpprev)){
	    #!print OUT "$polyline\n";
	    #not currently printing out because need to add parental bps to the counts file
	}#fell between prev and current aim
    }#run through poly lines

    if($tracker eq 0){
	if($distance >= $read_length){
	print AIMS "$line\n";
	print OUT "$elements[0]\t$elements[1]\t$par1\t0\t0\t$par2\t$curr_rate\n";
	$bpprev=$currbp;
	}#one snp per read 
    }#line hasn't already been printed

    #print "$currbp\n";

}#run through aims
