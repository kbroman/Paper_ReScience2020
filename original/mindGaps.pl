#!/usr/local/bin/perl

# mindGaps.pl
# 
# Look at each pair of adjacent genes in Data/mygeneloc.csv 
#     Same orientation?
#     Length of gap?

$ifile = "Data/mygeneloc.csv";
$ofile = "Data/gaps.csv";

open(IN, $ifile) or die("Cannot read from $ifile");
open(OUT, ">$ofile") or die("Cannot write to $ofile");
print OUT ("MT1,MT2,same_orient,gap,orient1,orient2\n");

$line = <IN>; # header row
$line = <IN>; chomp($line); # first gene
($mt1, $start1, $end1, $orient1) = split(/,/, $line);
while($line = <IN>) {
    chomp($line);
    ($mt2,$start2,$end2,$orient2) = split(/,/, $line);

    printf OUT ("%s,%s,%d,%d,%d,%d\n", $mt1, $mt2, ($orient1 == $orient2),
	       $start2-$end1-1, $orient1, $orient2);

    ($mt1,$start1,$end1,$orient1) = ($mt2,$start2,$end2,$orient2);
}
close(IN);
close(OUT);
