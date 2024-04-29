#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- Merged.group.lst.RefAllele.adjustPvalue.txt
my $in1 = $ARGV[1]; ##- Merged.group.lst.MissAllele.adjustPvalue.txt
my $in2 = $ARGV[2]; ##- Pangenome.blocks.txt.add.index
my $in3 = $ARGV[3]; ##- Pan.blocks.gff.coords
my $in4 = $ARGV[4]; ##- PanGene_Tair10.SynOrths.results.txt
my $in5 = $ARGV[5]; ##- Pangenome.blocks.txt.add.index.valid.Pos.txt

my $threshold = 1e-5;
my $flankLen = 5000;

##---------------------------------------------------------------------------------
my %panIndex2Sam   = ();
my %coord2geneName = ();
   &readCoords($in2, $in3, \%panIndex2Sam, \%coord2geneName);

my %genename2info = ();
   &gene2Description($in4, \%genename2info);

my %ValidPos = ();
   &readValidPos($in5, \%ValidPos);

my $outfile = "merged.outliers.signals.txt";
   &output($in0, $in1, \%panIndex2Sam, \%coord2geneName, \%genename2info, \%ValidPos, $outfile);

##------------all subroutines---------------------------------------------------------
sub output {

    my ($refPvalue, $missingPvalue, $panIndex2Sam, $coord2geneName, $genename2info, $ValidPos, $out) = @_;


    my $refPvalueCandidate = $refPvalue.".tmp";
       &extractCandidate($refPvalue, $threshold, $refPvalueCandidate);
    
    my $missingPvalueCandidate = $missingPvalue.".tmp";
       &extractCandidate($missingPvalue, $threshold, $missingPvalueCandidate);

         
    open OUT0, ">$out.tmp";
    ##---------------------------------------------------------------------------------------------------------------------------------------------------------------------
    my %recordgenes = ();
    open IN2, $refPvalueCandidate;
    while(<IN2>){
      chomp;
      my @temp  = split(/\t/, $_);

      my $samid = $panIndex2Sam ->{$temp[0]};
      my ($chr, $s, $e) = ($temp[1], $temp[2]-$flankLen, $temp[3]+$flankLen);
      my @geneArray = ();
      for my $tmps(keys %{$coord2geneName ->{$samid} ->{$chr}}){
          for my $tmpe(keys %{$coord2geneName ->{$samid} ->{$chr} ->{$tmps}}){
              
              my $flag   = &checkOverlap($ValidPos, $temp[0], $chr, $tmps, $tmpe, $s, $e);              
              my $geneid = $coord2geneName ->{$samid} ->{$chr} ->{$tmps} ->{$tmpe};

              if($flag == 1 && not exists $recordgenes{$geneid}){
                 push(@geneArray, $geneid);
              }
          }
      }
      if(scalar(@geneArray) > 0){
         for my $eachgene(@geneArray){
             my $info = join("\t", "/", "/");
             if(exists $genename2info ->{$eachgene}){
                $info = $genename2info ->{$eachgene};
             }
             print OUT0  join("\t", $temp[0], $samid, $chr, $temp[2], $temp[3], $temp[4], $temp[5], $temp[6], $temp[7], $temp[8], $temp[9], "RefAllele", $eachgene, $info), "\n";
             $recordgenes{$eachgene} = "Y";
         }
      }
      else{
         next; ##- fragments without gene models
         #print  OUT0  join("\t", $temp[0], $samid, $chr, $s, $e, $temp[4], $temp[5], $temp[6], $temp[7], $temp[8], $temp[9], "RefAllele"), "\n";   
      }
    }
    close IN2;    
    %recordgenes = ();

    ##------------------------------------------------------------------------------------------------------------------------------------------------------------------
    open IN3, $missingPvalueCandidate;
    while(<IN3>){
      chomp;
      my @temp  = split(/\t/, $_);
     
      my $samid = $panIndex2Sam ->{$temp[0]};
      my ($chr, $s, $e) = ($temp[1], $temp[2]-$flankLen, $temp[3]+$flankLen);
      my @geneArray = ();
      for my $tmps(keys %{$coord2geneName ->{$samid} ->{$chr}}){
          for my $tmpe(keys %{$coord2geneName ->{$samid} ->{$chr} ->{$tmps}}){

              my $flag = &checkOverlap($ValidPos, $temp[0], $chr, $tmps, $tmpe, $s, $e);
              my $geneid = $coord2geneName ->{$samid} ->{$chr} ->{$tmps} ->{$tmpe};

              if($flag == 1 && not exists $recordgenes{$geneid}){
                 push(@geneArray, $geneid);
              }
          }
      }
      if(scalar(@geneArray) > 0){
         for my $eachgene(@geneArray){
             my $info = join("\t", "/", "/");
             if(exists $genename2info ->{$eachgene}){
                $info = $genename2info ->{$eachgene};
             }
             print OUT0  join("\t", $temp[0], $samid, $chr, $temp[2], $temp[3], $temp[4], $temp[5], $temp[6], $temp[7], $temp[8], $temp[9], "MissingAllele", $eachgene, $info), "\n";
         }
      }
      else{
         next; ##- fragments without gene models
      }
    }
    close IN2;
    ##---------------------------------------------------------------------------------------------------------------------------------------------------------------------
    close OUT0;
     
    system("sort -k11,11g   $out.tmp  >  $out");
    system("rm -rf $out.tmp  $refPvalueCandidate  $missingPvalueCandidate");

}


sub extractCandidate {

    my ($fileIn, $cutoff, $fileOut) = @_;

    open IN0X, $fileIn;
    <IN0X>;
    open OUTX, ">$fileOut";
    while(<IN0X>){
      chomp;
      my @temp = split(/\t/, $_);

      if($temp[9] <= $cutoff){
         my ($Aaverage, $Asd) = split(/:/, $temp[6]);
         my ($Baverage, $Bsd) = split(/:/, $temp[7]); 
         
         my ($minAverage, $minSD) = ($Aaverage, $Asd); 
         my ($maxAverage, $maxSD) = ($Baverage, $Bsd);  
         if($Aaverage > $Baverage) {
            ($minAverage, $minSD) = ($Baverage, $Bsd);
            ($maxAverage, $maxSD) = ($Aaverage, $Asd);
         }    
         my $min3sigama = $minAverage+3*$minSD;     ##- it depends_ u 3 xigama
         my $max3sigama = $maxAverage-3*$maxSD;     ##- it depends_ u 3 xigama 
         if($max3sigama > $min3sigama){
            print OUTX $_, "\n"; 
         } 
      }
    }
    close IN0X;
    close OUTX;

}

sub checkOverlap {

    my ($validpos, $panindexName, $chr, $s1, $s2, $ns1, $ns2) = @_;

    my $flag = 0;

    my @data = sort {$a<=>$b} ($s1, $s2, $ns1, $ns2);
    my $length = ($s2-$s1+1) + ($ns2-$ns1+1);

    if(($data[3]-$data[0]+1) >= $length){
       $flag = 0;
    }
    else{
       for my $posTmp(keys %{$validpos ->{$panindexName} ->{$chr}}){
           if($posTmp >= $s1 && $posTmp <= $s2){
              $flag = 1;
              last;
           }
       }  
    }
    return($flag);

}

sub gene2Description {

    my ($fileIn, $gene2info) = @_;
 
    open IN1, $fileIn;
    <IN1>;
    while(<IN1>){
      chomp;
      my @temp = split(/\t/, $_);
   
      my $info = $temp[1];
      if(defined $temp[2]){
         $info .= "\t".$temp[2];
      }
      $gene2info ->{$temp[0]} = $info;
    }
    close IN1;

}

sub readCoords {

    my ($panindexFile, $genecoords, $panIndex2Sam, $coord2geneName) = @_;
 
    open IN0, $panindexFile;
    <IN0>;
    while(<IN0>){
      chomp;
      my @temp = split(/\t/, $_);
         $panIndex2Sam ->{$temp[0]} = $temp[5];
    }
    close IN0;

    open IN1, $genecoords;
    while(<IN1>){
      chomp;
      my @temp = split(/\t/, $_);
         $coord2geneName ->{$temp[1]} ->{$temp[2]} ->{$temp[3]} ->{$temp[4]} = $temp[0];
    }
    close IN1;

}

sub readValidPos {

    my ($fileIn, $posindex) = @_;
 
    open IN2XXXX, $fileIn;
    while(<IN2XXXX>){
      chomp;
      my @temp = split(/\t/, $_);
         $posindex ->{$temp[0]} ->{$temp[1]} ->{$temp[2]} = "Y";
    }
    close IN2XXXX;

}
