#!/usr/local/bin/perl

# exploreSeq.pl
#
# Parse the sequence data, identifying the genes' starts 
# and ends and the location of TAs, and then allow the 
# user to view the sequence of particular genes.


$ifile1 = "Rawdata/TIGR/GMT.1con";
    # complete sequence (73400 lines; ignore the first line)

$ifile2 = "Rawdata/TBCDC1551_rev.csv"; 
    # gene name, start, end, number of TAs

$ifile3 = "Rawdata/MTCoords_rev.csv";
    # gene name, start, end, number of TAs

# if genome not unzipped, unzip it
unless(-e $ifile1) {
    print(" -Unzipping\n");
    system("gunzip Rawdata/TIGR/GMT.1con.gz");
    $rungzip=1;
}

# read genome into a single big string
open(IN, $ifile1) or die("Cannot read from $ifile1");
$line = <IN>;
$data = "";
while($line = <IN>) {
    chomp($line);
    $data .= $line;
}
close(IN);
$n = length($data);
print("Length of genome (bp): $n\n");

# read information on genes' starts and ends + orientation
open(IN, $ifile3) or die("Cannot read from $ifile3");
while($line = <IN>) {
    chomp($line);
    ($mt, $start, $end) = split(/,/, $line);
    $mt = substr($mt, 2);
    $mt =~ s/^0+//;
    if($start > $end) {  
	$orient = -1;
	($start,$end) = ($end,$start);
    }
    else { $orient = 1; }
    $start2{$mt} = $start;
    $end2{$mt} = $end;
    $len2{$mt} = $end - $start + 1;
    $orient{$mt} = $orient;
}
close(IN);

# Read the other version of this file (w/o orientation info)
open(IN, $ifile2) or die("Cannot read from $ifile2");
$line = <IN>;
while($line = <IN>) {
    chomp($line);
    ($mt, $start, $end, $numta) = split(/,/, $line);
    $mt = substr($mt, 2);
    $mt =~ s/^0+//;
    # make sure the two files agree
    if($start != $start2{$mt} or $end != $end2{$mt}) {
	printf("Problem: %-10s %10d %10d    %10d %10d\n", $mt, 
	       $start, $start2{$mt}, $end, $end2{$mt});
    }

    # include stop codon
    if($orient{$mt} == -1) { 
	$start = $start - 3; 
    }
    else { 
	$end = $end + 3; 
    }

    $start{$mt} = $start;
    $end{$mt} = $end;
    $len{$mt} = $end - $start + 1;
    $numta{$mt} = $numta;
    push(@mt, $mt);
}
close(IN);
$nmt = @mt;
print("Number of genes: $nmt\n");

# if we unzipped the genome here, re-zip it
if($rungzip == 1) {
    print(" -Zipping\n");
    system("gzip $ifile1");
}

while(0 == 0) {
    print("\n\n");
    print("Enter gene: ");
    $input = <STDIN>; chomp($input);
    if($input eq "quit") { die("\n"); }

    print("\n    Orientation: $orient{$input}");
    print("\n    Length:      $len{$input}\n");
    print(substr($data, $start{$input}-21, 20), "\n\n");
    print(substr($data,$start{$input}-1,$len{$input}), "\n\n");
    print(substr($data, $end{$input}, 20), "\n");
}
