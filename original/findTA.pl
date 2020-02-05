#!/usr/local/bin/perl

# findTA.pl
#
# Use the genomic sequence and the gene start and end
# to identify the location of each TA within each gene.


$ifile1 = "Rawdata/TIGR/GMT.1con";
    # complete sequence (73400 lines; ignore the first line)

$ifile2 = "Rawdata/TBCDC1551_rev.csv"; 
    # gene name, start, end, number of TAs

$ifile3 = "Rawdata/MTCoords_rev.csv";
    # gene name, start, end, number of TAs

$ifile4 = "Rawdata/MT-RvConversion_rev.csv";
    # RV#, Product, Gene, Sub-classification, MT, MT#, #TA, TotalTA/Group

$ifile5 = "Rawdata/MtbGeneClassification.csv";
    # RV#, Product, Gene, Sub-classification

# if genome not unzipped, unzip it
unless(-e $ifile1) {
    print(" -Unzipping\n");
    system("gunzip Rawdata/TIGR/GMT.1con.gz");
    $rungzip=1;
}

# header rows in output files
$ofile = "Data/TAloc.csv";
$ofileb = "Data/geneinfo.csv";
#$ofile2 = "Data/diff_num_ta.txt";
#$ofile3 = "Data/diff_genes.txt";
$ofile4 = "Data/noTAs.txt";
$ofile5 = "Data/start_n_end.txt";
$ofile6 = "Data/mygeneloc.csv";
$ofile7 = "Data/overlaps.txt";
$ofile8 = "Data/doubleTA.csv";

open(OUT, ">$ofile") or die("Cannot write to $ofile");
open(OUTb, ">$ofileb") or die("Cannot write to $ofileb");
#open(OUT2, ">$ofile2") or die("Cannot write to $ofile2");
#open(OUT3, ">$ofile3") or die("Cannot write to $ofile3");
open(OUT5, ">$ofile5") or die("Cannot write to $ofile5");
open(OUT6, ">$ofile6") or die("Cannot write to $ofile6");
open(OUT7, ">$ofile7") or die("Cannot write to $ofile7");
open(OUT8, ">$ofile8") or die("Cannot write to $ofile8");

print OUT ("gene,TAloc,ORFsize\n");
print OUTb ("gene,length,class\n");
#print OUT2 ("gene  original  karl  diff\n");
print OUT5 ("MTnum   orient     start    stop\n");
print OUT6 ("gene,start,end,orient\n");
print OUT8 ("gene1,taloc1,gene2,taloc2\n");

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

# Read subclass information (for each RV number)
print(" -Reading classifications\n");
open(IN, $ifile5) or die("Cannot read from $ifile5");
$line = <IN>;
while($line = <IN>) {
    chomp($line);
    ($rv,$prod,$gene,$subclass) = split(/,/, $line);

    # if this rv number has been seen already, make sure it has the same class
    if($rv_subclass{$rv} ne "" and $rv_subclass{$rv} != $subclass) {
	print("REPEAT CLASS: rv: $rv \t");
	print("class:$subclass \told class: $rv_subclass{$rv}\n");
    }

    # save subclass
    $rv_subclass{$rv} = $subclass;
}
close(IN);

# read the MT <-> RV conversion file
open(IN, $ifile4) or die("Cannot read from $ifile4");
$line = <IN>;
while($line = <IN>) {
    chomp($line);
    ($rv,$prod,$gene,$subclass,$junk,$mt,$ta,$junk2) = split(/,/, $line);
    $junk =~ s/^MT(0*)//; # take mt number to be the column with MT####
    $mt = $junk;

    if($mt eq "") { next; }     # if not mt number, skip this line

    # if this MT was already assigned to an RV, make sure it's the same class
    if($rv{$mt} ne "" and $rv_subclass{$rv} ne $rv_subclass{$rv{$mt}}) {
	print("REPEAT: mt: $mt \trv:$rv  \trvold:$rv{$mt}\n");
    }
    # if this RV number was already seen, report that.
    if($mt{$rv} ne "") {
	print("REPEAT: rv: $rv \tmt:$mt  \tmtold:$mt{$rv}\n");
    }

    # save info
    $rv{$mt} = $rv;
    $mt{$rv} = $mt;
    $subclass{$mt} = $rv_subclass{$rv};
}
close(IN);
print(" -End classifications\n");

# look through genome for all TA locations
#     save location of the T in the TA
foreach $i (1..$n-1) {
    if(substr($data,$i-1,2) eq "TA") {
	push(@taloc, $i);
    }
}
$N = @taloc; 
print("Number of TA's: $N\n");

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
    if($orient{$mt} == -1) { $start -= 3; }
    else { $end += 3; }

    $start{$mt} = $start;
    $end{$mt} = $end;
    $len{$mt} = $end - $start + 1;
    $numta{$mt} = $numta;
    push(@mt, $mt);
}
close(IN);
$nmt = @mt;
print("Number of genes: $nmt\n");

# go through each gene and figure out the locations of any TAs
$next = 0;
foreach $gene (@mt) {
    if($subclass{$gene} eq "") {
	$subclass{$gene} = "noRV";
    }
    print OUTb ("$gene,$len{$gene},$subclass{$gene}\n");

    print OUT6 ("$gene,$start{$gene},$end{$gene},$orient{$gene}\n");

    # look at start and end codons
    if($orient{$gene} == 1) {
	$a = substr($data,$start{$gene}-1,3);
	$b = substr($data,$end{$gene}-3,3);
    }
    else {
	$a = join("", complement(substr($data, $end{$gene}-1, 1),
				 substr($data, $end{$gene}-2, 1),
				 substr($data, $end{$gene}-3, 1)));
	$b = join("", complement(substr($data, $start{$gene}+1, 1),
				 substr($data, $start{$gene}, 1),
				 substr($data, $start{$gene}-1, 1)));
    }
    if(($a eq "TTG" or $a eq "ATG" or $a eq "GTG") and 
       ($b eq "TAA" or $b eq "TAG" or $b eq "TGA")) { 
	($keepstart{$a})++; # keep track of start and stop codon usage
	($keepend{$b})++;
    }
    elsif($a eq "TTG" or $a eq "ATG" or $a eq "GTG") { # stop codon different than expected
	printf OUT5 ("%7s    %2d               %3s\n", $gene, 
		    $orient{$gene}, $b);
	($keepstart{$a})++;
    }
    elsif($b eq "TAA" or $b eq "TAG" or $b eq "TGA") { # start codon diff't than expected
	printf OUT5 ("%7s    %2d        %3s       \n", $gene, 
		    $orient{$gene}, $a);
	($keepend{$b})++;
    }
    else { # start and stop codons different than expected
	printf OUT5 ("%7s    %2d        %3s    %3s\n", $gene,
		     $orient{$gene}, $a, $b);
    }


    $n = 0;
    if($next < 0) { $next = 0; }
    # look through TAs, starting 50 back
    foreach $i ($next..($N-1)) {
	$taloc = @taloc[$i];
	$next = $i-50; 

	# is this TA within the gene?
	if(($orient{$gene} == 1 and $taloc >= $start{$gene}-1 and
	    $taloc < $end{$gene}-1) or
	   ($orient{$gene} == -1 and $taloc >= $start{$gene}+1 and 
	    $taloc < $end{$gene}+1)) {

	    if($orient{$gene} == -1) {
		$taloc = $end{$gene} - $taloc+1;
	    }
	    else {
		$taloc = $taloc-$start{$gene} + 2;
	    }
	    print OUT ("$gene,$taloc,$len{$gene}\n"); 
	    $n++;
	}
	elsif($taloc > $end{$gene}) { last; }
    }
    $overallTA += $n;

    # print info on discrepancies with old TA counts
#    if($n != $numta{$gene}) {
#	printf OUT2 ("%-10s %3d   %3d    %2d\n", $gene, 
#		     $numta{$gene}, $n, $n-$numta{$gene});
#
#	print OUT3 ("$gene :: $numta{$gene} $n    $orient{$gene}\n");
#	@geneseq = split(//, substr($data, $start{$gene}-1, $len{$gene}));
#	if($orient{$gene} == -1) {
#	    @geneseq = (complement(substr($data, $end{$gene}, 1)), " ",
#			complement(reverse(@geneseq))); 
#	}
#	else {
#	    @geneseq = (substr($data, $start{$gene}-2, 1), " ", @geneseq);
#	}
#	foreach $i (0..(@geneseq-1)) {
#	    print OUT3 ("$geneseq[$i]");
#	}
#	print OUT3 ("\n\n");
#
#    }
    
    if($n == 0) { # no TAs: print out
	print OUT4 ("$gene\n");
    }
	
}
print(" -TA's in genes: $overallTA\n");


# look for overlapping genes
print(" -Looking at overlapping regions\n");
foreach $i (1..(@mt-1)) {

    if($start{$mt[$i]} < $end{$mt[$i-1]}) { # overlap?
	if($end{$mt[$i-1]} > $end{$mt[$i]}) { $a = "*"; }
	else { $a = " "; }

	# look for TAs
	if($orient{$mt[$i-1]} == 1 and $orient{$mt[$i]} == 1) {
	    $overlap = substr($data, $start{$mt[$i]}-2, 
			      $end{$mt[$i-1]} - $start{$mt[$i]}+1);
	    $orientation = "pp";
	}
	elsif($orient{$mt[$i-1]} == -1 and $orient{$mt[$i]} == -1) {
	    $overlap = substr($data, $start{$mt[$i]}, 
			      $end{$mt[$i-1]} - $start{$mt[$i]}+1);
	    $orientation = "mm";
	}
	elsif($orient{$mt[$i-1]} == 1 and $orient{$mt[$i]} == -1) {
	    $overlap = substr($data, $start{$mt[$i]}, 
			      $end{$mt[$i-1]} - $start{$mt[$i]}-1);
	    $orientation = "pm";
	}
	else {
	    $overlap = substr($data, $start{$mt[$i]}-2, 
			      $end{$mt[$i-1]} - $start{$mt[$i]}+3);
	    $orientation = "mp";
	}

	$numta = 0;
	foreach $j (0..(length($overlap)-1)) {
	    if(substr($overlap,$j,2) eq "TA") { # a double-counted TA!
		$numta++;

		if($orientation eq "pp") {
		    $temp = $j + $start{$mt[$i]}-1;
		    print OUT8 ($mt[$i-1],",",$temp-$start{$mt[$i-1]}+2,",",
				$mt[$i],",",$temp-$start{$mt[$i]}+2,"\n");
		}
		elsif($orientation eq "mm") {
		    $temp = $j + $start{$mt[$i]}+1;
		    print OUT8 ($mt[$i-1],",",$end{$mt[$i-1]}-$temp+1,",",
				$mt[$i],",",$end{$mt[$i]}-$temp+1,"\n");
		}
		elsif($orientation eq "pm") {
		    $temp = $j + $start{$mt[$i]}+1;
		    print OUT8 ($mt[$i-1],",",$temp-$start{$mt[$i-1]}+2,",",
				$mt[$i],",",$end{$mt[$i]}-$temp+1,"\n");
		}
		else {
		    $temp = $j + $start{$mt[$i]}-1;
		    print OUT8 ($mt[$i-1],",",$end{$mt[$i-1]}-$temp+1,",",
				$mt[$i],",",$temp-$start{$mt[$i]}+2,"\n");
		}
	    }
	}
	$doubleTA += $numta;

	printf OUT7 ("%7s %10d %10d %4d %2d        %5d%1s\n",
		     $mt[$i-1], $start{$mt[$i-1]}, $end{$mt[$i-1]}, 
		     $len{$mt[$i-1]}, $orient{$mt[$i-1]},
		     $end{$mt[$i-1]} - $start{$mt[$i]} + 1, $a);
	printf OUT7 ("%7s %10d %10d %4d %2d        %5d\n",
		     $mt[$i], $start{$mt[$i]}, $end{$mt[$i]}, 
		     $len{$mt[$i]}, $orient{$mt[$i]}, $numta);
	print OUT7 (substr($data, $start{$mt[$i]}-2, 1), " ");
	print OUT7 (substr($data, $start{$mt[$i]}-1, 
			   $end{$mt[$i-1]} - $start{$mt[$i]}+1), " ");
	print OUT7 (substr($data, $end{$mt[$i-1]}, 1), "\n\n");


    }
}
print(" -Double-counted TA's: $doubleTA\n");

print OUT5 ("\n\nStarts:\n");
$num = 0;
foreach $start (keys %keepstart) {
    print OUT5 ("    $start    $keepstart{$start}\n");
    $num += $keepstart{$start};
}
print OUT5 ("    TOTAL    $num\n");
print OUT5 ("\n\nStops:\n");
$num = 0;
foreach $end (keys %keepend) {
    print OUT5 ("    $end    $keepend{$end}\n");
    $num += $keepend{$end};
}
print OUT5 ("    TOTAL    $num\n");

close(OUT);
#close(OUT2);
#close(OUT3);
close(OUT4);
close(OUT5);
close(OUT6);
close(OUT7);
close(OUT8);

# if we unzipped the genome here, re-zip it
if($rungzip == 1) {
    print(" -Zipping\n");
    system("gzip $ifile1");
}

sub complement {
    @a = @_;
    foreach $i (0..(@a-1)) {
	if($a[$i] eq "A") { $a[$i] = "T"; }
	elsif($a[$i] eq "T") { $a[$i] = "A"; }
	elsif($a[$i] eq "C") { $a[$i] = "G"; }
	else { $a[$i] = "C"; }
    }
    @a;
}
