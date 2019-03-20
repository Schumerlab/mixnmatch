#perl! -w

use Math::Random qw(:all);

my $indiv=shift(@ARGV); chomp $indiv;

my $aims=shift(@ARGV); chomp $aims;
open AIMS, $aims or die "cannot open parental AIMs file\n";

my $frequencies=shift(@ARGV); chomp $frequencies; 
open FREQS, $frequencies or die "cannot open parental AIM frequencies\n";

my $poly_par1=shift(@ARGV); chomp $poly_par1;

my $poly_par2=shift(@ARGV); chomp $poly_par2;

my $chr_length=shift(@ARGV); chomp $chr_length;

my $chr=shift(@ARGV); chomp $chr;

my $genome1=shift(@ARGV); chomp $genome1;
my $genome2=shift(@ARGV); chomp $genome2;

my $insnp1="indiv"."$indiv"."_hap1".".insnp";
open INSNP1, ">$insnp1";
my $insnp2="indiv"."$indiv"."_hap2".".insnp";
open INSNP2, ">$insnp2";
while((my $line1=<AIMS>) && (my $line2=<FREQS>)){

    chomp $line1;
    chomp $line2; 

    my @elements=split(/\t/,$line1);
    my $chrom=$elements[0]; my $bp=$elements[1]; my $species1=$elements[2]; my $species2=$elements[3];

    my @freqs=split(/\t/,$line2);
    
    my $par1_freq=$freqs[0]; chomp $par1_freq;
    my $par2_freq=$freqs[1]; chomp $par2_freq;
    #print "$par1_freq\t$par2_freq\n";

    my $focal_freq1=random_binomial(1,1,$par1_freq);

    if($focal_freq1 eq 1){
	print INSNP1 "$chrom\t$bp\t$species1\t$species2\n";
    }#update site hap1

    my $focal_freq2=random_binomial(1,1,$par2_freq);

    if($focal_freq2 eq 1){
	print INSNP2 "$chrom\t$bp\t$species2\t$species1\n";
    }#update site hap2

}#for all aims

my $add_poly_par1=int($poly_par1*$chr_length); print "adding $add_poly_par1 polymorphisms to parent 1 sequence\n";

my $add_poly_par2=int($poly_par2*$chr_length); print "adding $add_poly_par2 polymorphisms to parent 2 sequence\n";

my @polys_par1=random_uniform_integer($add_poly_par1,1,$chr_length);
my @polys_par2=random_uniform_integer($add_poly_par2,1,$chr_length);

my $file1poly="$insnp1".".poly"; chomp $file1poly;
my $file2poly="$insnp2".".poly"; chomp $file2poly;
open POLY1, ">"."$file1poly";
open POLY2, ">"."$file2poly";

for my $i (0..scalar(@polys_par1)-1){

    print POLY1 "$chr".":"."$polys_par1[$i]\n";

}#print out polymorphisms parent1

for my $j (0..scalar(@polys_par2)-1){

    print POLY2 "$chr".":"."$polys_par2[$j]\n";

}#print out polymorphisms parent1  

my $file1curr="$insnp1".".anc";
my $file2curr="$insnp2".".anc";
system("cat $file1poly | fastahack $genome1 -c > $file1curr");
system("cat $file2poly | fastahack $genome2 -c > $file2curr");

#print "cat $file1poly | fastahack $genome1 -c > $file1curr\n";
#print "cat $file2poly | fastahack $genome2 -c > $file2curr\n";

open POLYSITES1, "$file1curr" or die "cannot open sites to update for parent1\n";
open POLYSITES2, "$file2curr" or die "cannot open sites to update for parent2\n";

my @bp_array= qw(A T G C);

my $counter1=0;
while (my $poly1=<POLYSITES1>){
    chomp $poly1; my $shuffled_base=$poly1;
    while($poly1 eq $shuffled_base){
    my $rand_base=random_uniform_integer(1,0,3);
    $shuffled_base = $bp_array[$rand_base];
    #print "$shuffled_base\n";
    }#find another basepair
#    print "$chr\t$polys_par1[$counter]\t$poly1\t$shuffled_base\n";

    #sample the frequency of this polymorphism, constrained to be consistent across samples:
    my $freq_poly1=qx(head -n $counter1 polymorphism_distribution_par1 | tail -n 1); chomp $freq_poly1;
   
    my $change_base=random_binomial(1,1,$freq_poly1);
    if($change_base eq 1){
	print INSNP1 "$chr\t$polys_par1[$counter1]\t$poly1\t$shuffled_base\n";
        $counter1=$counter1+1;
    }#change this base in this individual

}#for all random polymorphic sites species1

my $counter2=0;
while (my $poly2=<POLYSITES2>){
    chomp $poly2; my $shuffled_base=$poly2;
    while($poly2 eq $shuffled_base){
	my $rand_base=random_uniform_integer(1,0,3);
	$shuffled_base = $bp_array[$rand_base];
    #print "$shuffled_base\n";                                                                                                     
    }#find another basepair

    my $freq_poly2=qx(head -n $counter2 polymorphism_distribution_par2 | tail -n 1); chomp $freq_poly2;

    my $change_base=random_binomial(1,1,$freq_poly2);
    if($change_base eq 1){
	print INSNP2 "$chr\t$polys_par2[$counter2]\t$poly2\t$shuffled_base\n";
	$counter2=$counter2+1;
    }#change this base in this individual 
   
}#for all random polymorphic sites species2

my $output1="indiv"."$indiv"."_hap1.fa";
my $output2="indiv"."$indiv"."_hap2.fa";
system("seqtk mutfa $genome1 $insnp1 > $output1");
system("seqtk mutfa $genome2 $insnp2 > $output2");

sub rand_sample {
    my ($n,@n) = (shift,@_);
    return 0 unless ($n < scalar @n);# see note below  
    my %seen = ();
    until (keys %seen == $samples) {
	$seen{$pop[rand @pop]}=1;
    }    
    return(keys %seen);
}
