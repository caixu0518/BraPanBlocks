#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- pixy_fst.txt

my %index = ();
open IN0, $in0;
<IN0>;
while(<IN0>){
  chomp;
  my @temp = split(/\t/, $_);
  my $midPos = ($temp[3]+$temp[4])/2;
  $index{$temp[0]}{$temp[1]}{$temp[2]}{$midPos} = $temp[5];
}
close IN0;

my @results = ();
for my $g1(sort keys %index){
    for my $g2(keys %{$index{$g1}}){
        
        my $tmpOut = $g1."vs".$g2.".out";
        open OUTX, ">$tmpOut";        
        for my $chr(sort keys %{$index{$g1}{$g2}}){
            for my $pos(sort {$a<=>$b} keys %{$index{$g1}{$g2}{$chr}}){
                print OUTX join("\t", $chr, $pos, $index{$g1}{$g2}{$chr}{$pos}, $g1."vs".$g2), "\n";
            }
        }     
        close OUTX;

        `perl /data/mg1/caix/scripts/XP-EHH/process_zscore.pl $tmpOut`;

        open INX, "$tmpOut.zscore";
        while(<INX>){
          chomp;
          push(@results, $_);
        }
        close INX;

        `rm $tmpOut $tmpOut.zscore`;  
    }
}


open OUT, ">pixy_fst.groups.lst";
print OUT join("\t", "Chr", "Pos", "zFst", "Fst", "Groups"), "\n";
print OUT join("\n", @results), "\n";
close OUT;

