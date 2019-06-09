#perl! -w

my $file1=shift(@ARGV); chomp $file1;
my $file2=shift(@ARGV); chomp $file2;

open OUT, ">macs_simulation_results_trees.aims.sharedpoly";

open F2, $file2 or die $!;
while (<F2>) {
    $h2{$_}++
};
open F1, $file1 or die;
$total=$.;
$printed=0;
while (<F1>) {
    $total++;
    if ($h2{$_}) {
	print OUT $_;
	$h2{$_} = "";
	$printed++;
	
    }
}
