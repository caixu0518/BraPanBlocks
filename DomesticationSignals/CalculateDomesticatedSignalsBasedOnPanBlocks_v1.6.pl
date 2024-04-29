#!/usr/bin/perl -w
use strict;
use Cwd;

my $in0       = $ARGV[0]; ##- Pangenome.blocks.txt.add.index
my $in1       = $ARGV[1]; ##- group.lst

my $scripts   = "/10t/caix/Bra_2023/PanVcfs/snp/domestication_analysis_2024/domestication_analysis_2024_new/scripts";
my $DataDir   = "/10t/caix/Bra_2023/PanVcfs/snp/domestication_analysis_2024/domestication_analysis_2024_new/BraPan_intervals_vcfs_file";

my $threadNum = 150;

my @indexArray = ();
   &file2Array($in0, \@indexArray);
   
   ##- calculate ref/missing alleles frequencies
   &generateCMD(\@indexArray, $threadNum);

   ##- calculate raw p values
   &calcutePvalue(\@indexArray);

   ##- merged together
   my $refAlleleFileOut = "Merged.".$in1.".RefAllele.Pvalue.txt";
   my $missAleleFileOut = "Merged.".$in1.".MissingAllele.Pvalue.txt";
      &mergedINVsPvalues(\@indexArray, $refAlleleFileOut, $missAleleFileOut);   

   my $RefPvalueOut  = "Merged."."$in1.RefAllele.adjustPvalue.txt";
   my $MissPvalueOut = "Merged."."$in1.MissAllele.adjustPvalue.txt";
      &AdjustPvalue($refAlleleFileOut, $RefPvalueOut);
      &AdjustPvalue($missAleleFileOut, $MissPvalueOut);
  
      #system("rm -rf  $refAlleleFileOut  $missAleleFileOut");

##--------  all subroutines -------------------------------
sub AdjustPvalue {

    my ($fileIn, $fileOut) = @_;

    my $tmpRawPvalueFile = $fileIn.".tmpRawPvalue.txt";

    open OUTPVALUE, ">$tmpRawPvalueFile";
    print OUTPVALUE "RawPvalue", "\n";
    open IN0, $fileIn;
    while(<IN0>){
      chomp;
      my @temp = split(/\t/, $_);
         print OUTPVALUE  $temp[7], "\n";
    }
    close IN0;
    close OUTPVALUE;

    ##- run p.adjust
    my $rTmpOut = $tmpRawPvalueFile.".tmp.R";
    my $fileOuttMP = $fileOut.".tmp";
    open OUTR, ">$rTmpOut";
    print OUTR  "data <- read.table(\"$tmpRawPvalueFile\", header = T, sep=\"\\t\")", "\n";
    print OUTR  "data\$adjustP <- p.adjust(data\$RawPvalue,  method =\"BH\")", "\n";
    print OUTR  "write.table(data, \"$fileOuttMP\", sep = \"\\t\", quote = FALSE, row.names = FALSE)", "\n";
    close OUTR;
    system("Rscript  $rTmpOut");
    system("rm -rf $rTmpOut");

    my @rawP = ();
    my @adjustP = ();
       &Pvaluefile2Array($fileOuttMP, \@rawP, \@adjustP);

    my $count = 0;
    open IN11, $fileIn;
    open OUT11, ">$fileOut";
    print OUT11 join("\t", "PanBlockId", "Chr", "Start", "End", "Group1", "Group2", "Frequency1", "Frequency2", "Pvalue", "adjustPvalue", "NegativeLog10"), "\n";
    while(<IN11>){
      chomp;
      my @temp = split(/\t/, $_);
      my $rawPvalue = $rawP[$count];
      my $adjustPvalue = $adjustP[$count];
      my $lgvalue = "NaN";
      if($adjustPvalue >0){
         $lgvalue = -1 * log($adjustPvalue)/log(10);
      }
      
      if($temp[7] == $rawPvalue){
         print OUT11  join("\t", $temp[8], $temp[0], $temp[1], $temp[2], $temp[3], $temp[4], $temp[5], $temp[6], $temp[7], $adjustPvalue, $lgvalue), "\n";
      }
      else{
         print OUT11  join("\t", @temp, "NaN", "NaN"), "\n";
      }
      $count += 1;
    }
    close IN11;
    close OUT11;
    system("rm -rf $fileOuttMP  $tmpRawPvalueFile");

}

sub Pvaluefile2Array {

    my ($fileIn, $rawP, $adjustP) = @_;

    open IN00, $fileIn;
    <IN00>;
    while(<IN00>){
      chomp;
      my @temp = split(/\t/, $_);
      push(@{$rawP}, $temp[0]);
      push(@{$adjustP}, $temp[1]);
    }
    close IN00;

}
    
sub mergedINVsPvalues {

    my ($PanIndexArray, $refAlleleFileOut, $missAleleFileOut) = @_;

    my $refAlleleFileOutTmp = $refAlleleFileOut.".tmpxx";
    my $missAleleFileOutTmp = $missAleleFileOut.".tmpxx";

    open OUTrefAlleleFileOutTmp, ">$refAlleleFileOutTmp";
    open OUTmissAleleFileOutTmp, ">$missAleleFileOutTmp";
    for my $tmpkey(@{$PanIndexArray}){
        my $tmpfileref = $tmpkey.".groups.alleles.freq.refAllele.Pvalue";    
        my $tmpfileMissing = $tmpkey.".groups.alleles.freq.missAllele.Pvalue"; 
       
        if(-e $tmpfileref){ 
           ##- process ref    
           open INtmpfileref, $tmpfileref;
           <INtmpfileref>;
           while(<INtmpfileref>){
             chomp;
             my @temp = split(/\t/, $_);
             if(defined $temp[0] && defined $temp[1] && defined $temp[2] && defined $temp[4] && defined $temp[5] && defined  $temp[6] &&  defined  $temp[7] && defined $temp[8]){
                next, if($temp[8] eq "nan");
                print OUTrefAlleleFileOutTmp  join("\t", $temp[0], $temp[1], $temp[2], $temp[4], $temp[5], $temp[6], $temp[7], $temp[8], $tmpkey), "\n";
             }
             else{
               #print "Warining(s): undefined values in the file: $tmpfileref", "\n";
             }
           }
           close INtmpfileref;              
        }
        if(-e $tmpfileMissing){
           ##- process missing
           open INtmpfileMissing, $tmpfileMissing;  
           <INtmpfileMissing>;
           while(<INtmpfileMissing>){
             chomp;
             my @temp = split(/\t/, $_);
             if(defined $temp[0] && defined $temp[1] && defined $temp[2] && defined $temp[4] && defined $temp[5] && defined  $temp[6] &&  defined  $temp[7] && defined $temp[8]){
                next, if($temp[8] eq "nan");
                print OUTmissAleleFileOutTmp  join("\t", $temp[0], $temp[1], $temp[2], $temp[4], $temp[5], $temp[6], $temp[7], $temp[8], $tmpkey), "\n";
             }
             else{
                #print "Warining(s): undefined values in the file: $tmpfileMissing", "\n";
             }
           }
           close INtmpfileMissing; 
        } 
    }
    close OUTrefAlleleFileOutTmp;
    close OUTmissAleleFileOutTmp;

    ##- sort
    system("sort -k8,8g $refAlleleFileOutTmp > $refAlleleFileOut");     
    system("sort -k8,8g $missAleleFileOutTmp > $missAleleFileOut");
    
    system("rm -rf $refAlleleFileOutTmp  $missAleleFileOutTmp");
    
    my $tmpOut = "TmpOutput";
    my $currentPath = getcwd();
    
    my $freq  = "$currentPath/$tmpOut/00.freq";
    my $piDir = "$currentPath/$tmpOut/01.pi";
    my $pvalues = "$currentPath/$tmpOut/02.RawPvalues";

       system("mkdir -p  $freq"),    if(not -e $freq);
       system("mkdir -p  $piDir"),   if(not -e $piDir);
       system("mkdir -p  $pvalues"), if(not -e $pvalues);
       for my $tmpKey(@{$PanIndexArray}){
           my $tmpfreqFile       = $tmpKey.".groups.alleles.freq";
           my $tmpfreqMissAllele = $tmpKey.".groups.alleles.freq.missAllele";
           my $freqrefAllele     = $tmpKey.".groups.alleles.freq.refAllele";
           if(-e $tmpfreqFile){
              system("pigz $tmpfreqFile");
              system("mv  $tmpfreqFile.gz  $freq");
           }
           if(-e $tmpfreqMissAllele){
              system("pigz $tmpfreqMissAllele");
              system("mv $tmpfreqMissAllele.gz  $freq");
           }
           if(-e $freqrefAllele){
              system("pigz $freqrefAllele");
              system("mv  $freqrefAllele.gz $freq");
           }
           my $pifile = $tmpKey.".groups.pi"; 
           if(-e  $pifile){
              system("pigz  $pifile");
              system("mv  $pifile.gz  $piDir");
           }

           my $pvalueRef = $tmpKey.".groups.alleles.freq.refAllele.Pvalue";
           my $pvalueMissing = $tmpKey.".groups.alleles.freq.missAllele.Pvalue";
           system("mv  $pvalueRef  $pvalues"), if(-e $pvalueRef);
           system("mv  $pvalueMissing  $pvalues"), if(-e $pvalueMissing);
       } 

}

sub calcutePvalue {

    my ($PanIndexArray) = @_;

    my $tmpRunPvalueCMDs = "Temp.Run.Pvalues.cmds";
    system("rm -rf $tmpRunPvalueCMDs"), if(-e $tmpRunPvalueCMDs);

    open OUTCMDPvalue, ">$tmpRunPvalueCMDs";
    for my $tmpKey(@{$PanIndexArray}){
        my $tmpFile = $tmpKey.".groups.alleles.freq";
  
        if(-e $tmpFile){
           my $samLineCount = 0;
              $samLineCount = &samLine($tmpFile);          
              if($samLineCount > 1){
                 print OUTCMDPvalue "perl  $scripts/calculate_RawPvalue.pl  $tmpFile  $in0", "\n";
              }
        }
    }
    close OUTCMDPvalue;
    system("ParaFly -c $tmpRunPvalueCMDs -CPU  $threadNum");
    system("rm -rf $tmpRunPvalueCMDs*");

}

sub samLine {

    my ($fileIn) = @_;

    my $count = 0;
    open INXXX, $fileIn;
    while(<INXXX>){
      $count += 1; 
    }
    close INXXX;
    return($count);

}

sub generateCMD {

    my ($fileIndexArray, $threadNum) = @_;

    my %finished = ();
    my $count = scalar(@{$fileIndexArray});
    my $round = 1;
    for($round=1; ($round*$threadNum)<=($count-1); $round++){
               
        my ($start, $end) = (($round-1)*$threadNum, $round*$threadNum);  ##- not include end

        ##- generate ParaFly cmommands
        my $cmds = "batch.generate_Pi_AlleleFreq.round.$round.cmds";
        system("rm -rf  $cmds"), if(-e $cmds);
        system("rm -rf  $cmds.compelated"), if(-e "$cmds.compelated");

        my $count = 0; 
        open OUT, ">$cmds";
        for(my $m=$start; $m<$end; $m++){
            my $indexTmp  = $fileIndexArray ->[$m];
               $finished{$indexTmp} = "Y";
            my $fileName  = $indexTmp.".vcf.gz";
            my $fileIndex = $indexTmp.".vcf.gz.tbi";
               system("ln -s  $DataDir/$fileName  .");
               system("ln -s  $DataDir/$fileIndex .");
               print OUT  "perl  $scripts/generate_allele_frequency_and_Pi_v1.1.pl   $in1   $fileName", "\n";         
               $count += 1;
        }
        close OUT; 
        
        ##- paraFly
        system("ParaFly -c $cmds -CPU $count"); 

        ##- clean
        for(my $m=$start; $m<=$end; $m++){
            my $fileName  = ($fileIndexArray ->[$m]).".vcf.gz";
            my $fileIndex = ($fileIndexArray ->[$m]).".vcf.gz.tbi";
            system("rm  $fileName  $fileIndex");
        } 
        system("rm -rf $cmds"), if(-e $cmds);
        system("rm -rf $cmds.completed"), if(-e "$cmds.completed");

     }
     ##- add supplymentary
     my %left = ();
     $round += 1;
     my $cmd1s = "batch.generate_Pi_AlleleFreq.round.$round.cmds";
     my $countTmp = 0;
     open OUTX, ">$cmd1s";
     for my $eachtmp(@{$fileIndexArray}){
         if(not exists $finished{$eachtmp}){
                $left{$eachtmp} = "Y";
             my $fileName1  = $eachtmp.".vcf.gz";
             my $fileIndex1 = $eachtmp.".vcf.gz.tbi";
                system("ln -s     $DataDir/$fileName1  .");
                system("ln -s     $DataDir/$fileIndex1 .");
                print OUTX  "perl  $scripts/generate_allele_frequency_and_Pi.pl  $in1  $fileName1", "\n";                   
                $countTmp += 1;
          }
      }
      close OUTX;  
             
      ##- paraFly
      system("ParaFly -c $cmd1s -CPU $countTmp");

      ##- clean
      for my $eachTmp(keys %left){
          system("rm -rf $eachTmp.vcf.gz $eachTmp.vcf.gz.tbi");
      }
      system("rm -rf $cmd1s"), if(-e  $cmd1s);
      system("rm -rf $cmd1s.completed"), if(-e  "$cmd1s.completed");
}

sub file2Array {

    my ($fileIn, $array) = @_;

    open IN0, $fileIn;
    while(<IN0>){
      chomp;
      my @temp = split(/\t/, $_);
      next, if($temp[1] <= 10);  ##- skip blocks with snp number less than 10 -----
      push(@{$array}, $temp[0]);
    }
    close IN0;

}
