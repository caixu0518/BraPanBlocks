#!/usr/bin/perl -w
use strict;
use Statistics::Descriptive;

my $in0 = $ARGV[0]; ##- Brapan03936.groups.alleles.freq
my $in1 = $ARGV[1]; ##- Pangenome.blocks.txt.add.index

my $winsize  = 10000;
my $stepsize = 5000;
my $minMaf   = 0.1;


my $refAlleleFile   = $in0.".refAllele";
   &generateRefAlleleFreq($in0, $refAlleleFile);

my $missAlleleFile  = $in0.".missAllele";
   &generateMissAlleleFreq($in0, $missAlleleFile);

my %RefratioIndex   = ();
my @groupArrayR     = ();
   &recordRatio($refAlleleFile, \@groupArrayR, \%RefratioIndex);

my %missratioIndex  = ();
my @groupArrayM     = ();
   &recordRatio($missAlleleFile, \@groupArrayM, \%missratioIndex);

my $panId  = "nan";
my $region = "nan";
   &extractRegion($refAlleleFile, $in1, \$panId, \$region);

my $PvalueOut = $refAlleleFile.".Pvalue";
   &calculateTest($region, $refAlleleFile, $PvalueOut, $winsize, $stepsize);

my $missOut   = $missAlleleFile.".Pvalue";
   &calculateTest($region, $missAlleleFile, $missOut, $winsize, $stepsize);

  # system("rm -rf $refAlleleFile  $missAlleleFile");

##--------------------------------------------------------------------

sub extractRegion {

    my ($fileIn, $paninfo, $panID, $regionInfo) = @_;
  
    my $indexName = $fileIn;
       $indexName =~ s/\.groups\.alleles\.freq\.refAllele//; 

    $$panID = $indexName;

    open IN5, $paninfo;
    while(<IN5>){
      chomp;
      my @temp = split(/\t/, $_);
      if($temp[0] eq $indexName){
         $$regionInfo = $temp[2].":".$temp[3]."-".$temp[4];
      }
    }
    close IN5;

}

sub calculateTest {

    my ($region, $fileIn, $filePvalueout, $stepSize, $winSize) = @_;

    open OUTPvalue, ">$filePvalueout";
    print OUTPvalue join("\t", "Chr", "Start", "End", "MidPos", "Control", "Derive", "AlleleFreqControl", "AlleleFreqDerived", "Pvalue"), "\n";
    ##- record allele
    my @groupsid    = ();
    my %alleleIndex = ();
       &recordAleleRatio($fileIn, \%alleleIndex, \@groupsid);

    if($region =~ /^(\S+?):(\S+?)-(\S+)$/){
       my ($chr, $start, $end) = ($1, $2, $3);

       ##- the length of intervals less than a window size
       my $s     = "nan";
       my $e     = "nan";
       my @resultArray = ();
       my $count = 0;
       if(($end-$start+1) <= $winSize){
           ($s, $e) = ($start, $end); 
           ##- record vars in the target region
           @resultArray = ();
           &calcualteTtestforEach(\%alleleIndex, \@groupsid, $chr, $s, $e, \@resultArray);
           if(scalar(@resultArray) > 0){
              print OUTPvalue join("\n", @resultArray), "\n";
           }
       }
       elsif(($end-$start+1) > $winSize){
           for($count=0; ($start+$stepSize*$count+$winSize)<=$end; $count++){
               ($s, $e) = ($start+$stepSize*$count, $start+$stepSize*$count+$winSize); 
               ##-                  
               @resultArray = (); 
               &calcualteTtestforEach(\%alleleIndex, \@groupsid, $chr, $s, $e, \@resultArray);
               print OUTPvalue join("\n", @resultArray), "\n";
           } 
       }
       else{
           die "error 2: xxxxx.\n";
       }
    } 
    else{
       die "error 1: the input region format is wrong.\n";
    }
    close OUTPvalue;

}

sub calcualteTtestforEach {
    
    my ($AlleleIndex, $groupsid, $chr, $startPos, $endPos, $results) = @_;
    
    my $groupCount = scalar(@{$groupsid});
    my $flag = 0; 
    my %tmpAllele = ();
    for my $each(@{$groupsid}){
        for my $pos(sort {$a<=>$b} keys %{$AlleleIndex ->{$each} ->{$chr}}){
            if($pos >= $startPos && $pos <= $endPos){
               if(not exists $tmpAllele{$each}){
                  $tmpAllele{$each} = $AlleleIndex ->{$each} ->{$chr} ->{$pos};
               }
               else{
                  $tmpAllele{$each} .= "\t".($AlleleIndex ->{$each} ->{$chr} ->{$pos});
               }
            }
        }
    }
    
    ##- skip ....
    my $controlId = $groupsid ->[($groupCount-1)]; 
    if(exists $tmpAllele{$controlId}){
    my @control   = split(/\t/, $tmpAllele{$controlId});
    if(scalar(@control) <= 10){
       @{$results} = ();
    }   
    else{ 
      my $rscriptPrefix = $chr.":".$startPos."_".$endPos;
      if(scalar(@control) > 10){
         my $statcontrol    = Statistics::Descriptive::Full->new();
            $statcontrol->add_data(@control); 
         my $averagecontrol = $statcontrol->mean();
         my $sdcontrol      = $statcontrol->standard_deviation();

         for(my $m=0; $m<=($groupCount-2); $m++){
             my $tmpid = $groupsid ->[$m]; 
             my @alleleArray = split(/\t/, $tmpAllele{$tmpid});
             if(scalar(@alleleArray) > 0){

                my $statQ        = Statistics::Descriptive::Full->new(); 
                   $statQ->add_data(@alleleArray);
                my $qveragegroup = $statQ->mean();
                my $sdqgroup     = $statQ->standard_deviation();

                my $derived2control = "nan";
                next, if(scalar(@control) < 10);
                next, if(scalar(@alleleArray) < 10);
 
                if($qveragegroup >= $averagecontrol){
                   $derived2control = &ttest(\@alleleArray, \@control, $rscriptPrefix, \$qveragegroup, \$sdqgroup, \$averagecontrol, \$sdcontrol);
                }
                else{
                   $derived2control = &ttest(\@control, \@alleleArray, $rscriptPrefix, \$averagecontrol, \$sdcontrol, \$qveragegroup, \$sdqgroup);
                }
                $averagecontrol = sprintf("%.2f", $averagecontrol);
                $qveragegroup   = sprintf("%.2f", $qveragegroup);

                my $tmpresult = join("\t", $chr, $startPos, $endPos, int(($startPos+$endPos)/2), $groupsid ->[-1], $groupsid ->[$m], $averagecontrol.":".$sdcontrol, $qveragegroup.":".$sdqgroup, $derived2control);
                push(@{$results}, $tmpresult);

             }
         }
      }
      else{
         $flag = 1;
      }
   }
   }
}

sub recordAleleRatio {
    
    my ($fileIn, $alleleIndex, $groupids) = @_;
   
       open IN2, $fileIn;
       my $firstLine = <IN2>;
          chomp($firstLine);
       my @first     = split(/\t/, $firstLine);
       for(my $m=2; $m<=$#first; $m++){
           push(@{$groupids}, $first[$m]);
       }
       while(<IN2>){
         chomp;
         my @temp = split(/\t/, $_);
         for(my $m=2; $m<=$#temp; $m++){
             my $id = $first[$m];
                $alleleIndex ->{$id} ->{$temp[0]} ->{$temp[1]} = $temp[$m];
         }
       }
       close IN2;
      
}

sub recordRatio {

    my ($fileIn, $groupArray, $ratioIndex) = @_;
 
    open IN1, $fileIn; 
    my $first = <IN1>;
       chomp($first);
    my @first = split(/\t/, $first);
    for(my $m=2; $m<=$#first; $m++){
        push(@{$groupArray}, $first[$m]);
    } 

    while(<IN1>){
      chomp;
      my @temp = split(/\t/, $_);
      next, if(scalar(@temp) != scalar(@first));
      for(my $n=2; $n<=$#temp; $n++){
          my $id = $first[$n];
          $ratioIndex ->{$id} ->{$temp[0]} ->{$temp[1]} = $temp[$n];
      }
    }
    close IN1;
    
}

sub generateMissAlleleFreq {

    my ($fileIn, $fileOut) = @_;
    
    open IN0, $fileIn;
    open OUT0, ">$fileOut";
    my $firstLine = <IN0>;
    print OUT0 $firstLine;
    while(<IN0>){
      chomp;
      my @temp = split(/\t/, $_);
 
      my @values = (); 
      for(my $m=2; $m<=$#temp; $m++){
          my $tmpValue = "nan";
          if($temp[$m] =~ /^(\S+?)\/(\S+?);(\S+?):(\S+?):(\S+?):(\S+)$/){
             my ($ref, $alt, $heter, $miss) = ($3, $4, $5, $6);
             
             next, if($ref+$alt+$heter+$miss <= 0);            
             $tmpValue = $miss/($ref+$alt+$heter+$miss);     ##- it depends_
             push(@values, $tmpValue); 
          }  
          else{
             next;
             #die "Error 1: cannot match $temp[$m].\n";
          }  
      }      
      print OUT0  join("\t", $temp[0], $temp[1], @values), "\n";
    }        
    close IN0;
    close OUT0;

}

sub generateRefAlleleFreq {

    my ($fileIn, $fileOut) = @_;

    open IN0, $fileIn;
    open OUT0, ">$fileOut";
    my $firstLine = <IN0>;
    print OUT0 $firstLine;
    while(<IN0>){
      chomp;
      my @temp = split(/\t/, $_);
      my @values = ();

      for(my $m=2; $m<=$#temp; $m++){
          my $tmpValue = "nan";
          if($temp[$m] =~ /^(\S+?)\/(\S+?);(\S+?):(\S+?):(\S+?):(\S+)$/){
             my ($ref, $alt, $heter, $miss) = ($3, $4, $5, $6);
             
             next, if($ref+$alt+$heter+$miss <= 0);             
             $tmpValue = $ref/($ref+$alt+$heter+$miss);     ##- it depends_
             push(@values, $tmpValue); 
          }
          else{
             next;
             #die "Error 1: cannot match $temp[$m].\n";
          }
      }
      print OUT0  join("\t", $temp[0], $temp[1], @values), "\n";
    } 
    close IN0;
    close OUT0;

}

sub ttest {

    my ($groupA, $groupB, $prefix, $averageA, $sdA, $averageB, $sdB) = @_;
 
    my $pvalue = "nan";

    ##- average (A) > average (B)
    my $Aelements = scalar(@{$groupA}); 
    my $Belements = scalar(@{$groupB});  

    if($Aelements != $Belements){
       $pvalue = "nan";
    }      
    else{
     
       my $fileTmp = $prefix.".ttest.r";
       my $Rout = $fileTmp.".r2.r.out";

       open OUT, ">$fileTmp";
       print OUT "sink(\"$Rout\")", "\n";
     
       my @Avalues = ();
       my @Bvalues = ();
       for(my $m=0; $m<=($Aelements-1); $m++){
           my $elementA = $groupA ->[$m];
           my $elementB = $groupB ->[$m];
 
           if($elementA <= 0.1 && $elementB <= 0.1){
              next;
           }
           elsif($elementA >= 0.9 && $elementB >= 0.9){
              next;
           }
           else{
              my $min = $elementA;
              my $max = $elementB;
              if($elementA > $elementB){
                 $max = $elementA;
                 $min = $elementB;
              }    
              $min = 0.0001,   if($min == 0);              
              $max = 0.0001,   if($max == 0);

              if($max/$min >= 2){
                 push(@Avalues, $elementA);
                 push(@Bvalues, $elementB);
              }
           }

       }
       print OUT  "groupA <- c(", join(",", @Avalues), ")", "\n";
       print OUT  "groupB <- c(", join(",", @Bvalues), ")", "\n";
       print OUT  "t.test(groupA, groupB, alternative = \"greater\", paired=T, conf.level=0.99)\$p.value", "\n";  
       print OUT  "sink()", "\n";
       
       if(scalar(@Avalues) >= 3) {
          system("Rscript $fileTmp");
           
          my $statA = Statistics::Descriptive::Full->new();
             $statA->add_data(@Avalues);
          my $statB = Statistics::Descriptive::Full->new(); 
             $statB->add_data(@Bvalues);

             $$averageA = $statA->mean(); 
             $$sdA      = $statA->standard_deviation(); 
             $$averageB = $statB->mean();
             $$sdB      = $statB->standard_deviation(); 

             open INX, $Rout;
             while(<INX>){
               chomp;
               my @temp = split(/\s+/, $_);
                  $pvalue = $temp[1];
             }
             close INX;
       }
       system("rm $fileTmp"), if(-e $fileTmp);
       system("rm $Rout"), if(-e $Rout);
    }
    return($pvalue);

}
