#perl! -w

my $file1=shift(@ARGV); chomp $file1;
my $file1=shift(@ARGV); chomp $file2;

($file1, $file2) = @ARGV;
$printed = 0;
$total = 0;
open F2, $file2 or die "$file2: $!\n";
@lines2 = <F2>;
$total += $.;
foreach (@lines2) {
    $h2{$_}++;
}
open F1, $file1 or die "$file1: $!\n";
while (<F1>) {
    if (exists $h2{$_}) {
	$h2{$_} = 0;
	
    }
    else {
	print $_;
	$printed++;
	
    }
}
$total += $.;
foreach (@lines2) {
    if ($h2{$_}) {
	print $_;
	$printed ++;
	
    }
}
