#!/usr/bin/perl -w
use strict;
use threads;

    ##- Xu Cai
    ##- Bug report: caixu0518@163.com

my $input = $ARGV[0]; ##- final.plk

my @thr = ();
   for(my $i=1;$i<=8;$i+=1) {
        $thr[$i] = threads->create(\&run_structure, $input, $i+1);
   }

   for(my $i=1;$i<=8;$i+=1) {
       $thr[$i]->join;
   }


##-----------------------------------------------------------------------
sub run_structure {
 
    my ($input, $k) = @_;

    my $cmdString = "/usr/bin/python    /fastStructure/structure.py  --input=$input  --output=final  -K $k";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
}
