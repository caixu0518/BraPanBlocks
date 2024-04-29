#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- Pangenome.blocks.txt.add


my $outfileDir = "BraPan_intervals_vcfs_file";
my $outLog     = $in0.".index";
   system("mkdir $outfileDir"), if(not -d $outfileDir);
   #system("rm -rf $outfileDir"), if(-d  $outfileDir);
   #system("mkdir  $outfileDir");
   &output($in0, $outLog);

##------------------------------------------------------------------------
sub output {

    my ($fileIn, $logFile) = @_;

    my $count = 0;
    open OUTX, ">$logFile";
    print OUTX join("\t", "PanIntervalName", "VarCount", "Chr", "Start", "End", "Reference", "GenesInInterval"), "\n";
    open IN0, $fileIn;
    while(<IN0>){
      chomp;
      $count += 1;
      my @temp         = split(/\t/, $_);
      my $newnum       = &numFormat($count);
      my $index        = "Brapan".$newnum;
      my $intervalVCFs = $index.".vcf"; 
      my $region       = $temp[0].":".$temp[1]."-".$temp[2];
     
      my $snpcount = 0;

         my $vcffile      = $temp[3].".pan.snp.vcf.gz";
         #print  "waining(s): cannot find: $vcffile", if(not -e $vcffile);

         if(-e $vcffile){
            system("bcftools view  $vcffile -r $region  > $intervalVCFs");

            if(-e $intervalVCFs){

               my $fileLine = &countFileline($intervalVCFs);
               next, if($fileLine == 0);
               system("bgzip $intervalVCFs");
               system("tabix -p vcf $intervalVCFs.gz");
               system("mv  $intervalVCFs.gz  $intervalVCFs.gz.tbi  $outfileDir");

               ##- print Log
               $snpcount = &varcount("$outfileDir/$intervalVCFs.gz");
               print OUTX join("\t", $index, $snpcount, $temp[0], $temp[1], $temp[2], $temp[3], $temp[4]), "\n";  ##- index, snp count , ....

            }
         }

    }
    close IN0;
    close OUTX;

}


sub countFileline {

    my ($fileIn) = @_;
 
    my $tmpCount = 0;
    open INX, $fileIn;
    while(<INX>){
      chomp;
      $tmpCount += 1;
    }
    close INX;    
    return($tmpCount);

}

sub varcount {

    my ($fileIn) = @_;

    my $count = 0;
    open IN1, "pigz -dc $fileIn | ";
    while(<IN1>){
      next, if(/^#/);
      $count += 1;     
    }
    close IN1;
    return($count);
 
}

sub numFormat {
 
    my ($n) = @_;
    
    my $m = $n;
    if($n >= 0 && $n < 10){
       $m = "0000".$n; 
    }
    elsif($n >= 10 && $n < 100){
       $m = "000".$n;
    }
    elsif($n >= 100 && $n < 1000){
       $m = "00".$n;
    }
    elsif($n >= 1000 && $n < 10000){
       $m = "0".$n;
    }
    else{
      $m = $n;
    }
    return($m); 

}
