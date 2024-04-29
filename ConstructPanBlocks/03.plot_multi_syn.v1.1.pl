#!/usr/bin/perl -w
use SVG; 

my $in0 = $ARGV[0]; ##- micro_for_eachGenome.bed  
my $out = $ARGV[1]; ##- A01

my $coordDir = "/10t/caix/Pan-genome/coords/";   

#my @color = ("#F8766D", "#00BFC4", "#B79F00", "#619CFF", "#00BA38", "#F564E3");
my @color = &colorSet();
my @samInfo = ();
my @sams = ();
my %sam2Len = ();
my %sam2Start = ();
my %sam2End = ();
   &getLongestChr($in0, \@samInfo, \@sams, \%sam2Len, \%sam2Start, \%sam2End);

my %sam2color = ();
   &sam2color(\@sams, \@color, \%sam2color);


#--Initialize picture----------------------------
my $width  = 1200;
my $height = 100 + 50*$samInfo[1] + 250;
my $svg    = SVG->new(width => $width, height => $height);
my $xstart = $width - 1100; 
my $ystart = 80;

my $y0     = 25;
   for(my $m=0; $m<=($#sams-1); $m++){
       my $samname   = $sams[$m];
       my $colorname = $sam2color{$samname};
       $svg -> rect(x => $xstart+50*$m, y => $y0, width => 20, height => 20,
                         fill => $colorname, "stroke-width"=>"0.1", 'stroke-opacity' => 0.2, 'fill-opacity' => 1);
       $svg -> text(x => $xstart+50*$m+10, y => $y0+30, 'font-family' => 'Arial', 'font-size' => 10,                  
                         'font-weight' => 'normal', stroke => 'black', 'stroke-width' => 0.2, 'text-anchor' => 'middle', -cdata => $samname);
   }   
 
#--draw scales-----------------------------------
my $numtmp  = int($samInfo[0]/1000000/5)+1;
my $largest = 5*$numtmp;
my $xratio  = ($width-200)/$largest/1000000;
   $svg -> line(x1 => $xstart, y1 => $ystart, x2 => $xstart+1000, y2 => $ystart, 'stroke' => 'black', 'stroke-width' => 2);

   my $Rcanvas = ($width-100-100)/$largest;
   for (my $i=0;$i<=$largest;$i++){                                                                                                      
                if ($i%5 == 0 ){                                                                                                              
                        $svg -> line(x1 => $xstart+$i*$Rcanvas, y1 => $ystart, x2 => $xstart+$i*$Rcanvas, y2 => $ystart+8,                    
                                 'stroke' => 'black', 'stroke-width' => 1.2);                                                                 
                        $svg -> text(x => $xstart+$i*$Rcanvas, y => $ystart+30, 'font-family' => 'Arial', 'font-size' => 15,                  
                                'font-weight' => 'normal', stroke => 'black', 'stroke-width' => 0.2, 'text-anchor' => 'middle', -cdata => $i);
                }else{                                                                                                                        
                        $svg -> line(x1 => $xstart+$i*$Rcanvas, y1 => $ystart, x2 => $xstart+$i*$Rcanvas, y2 => $ystart+5,                    
                                'stroke' => 'black', 'stroke-width' => 1.2);                                                                  
                }                                                                                                                             
   }                                                                                                                                     
   $svg -> text(x => $xstart+$width-100-100+25+20, y => $ystart+30, 'font-family' => 'Arial', 'font-size' => 15,                         
                        'font-weight' => 'normal', stroke => 'black', 'stroke-width' => 0.2, 'text-anchor' => 'middle', -cdata => 'Mb');

##- make the plot in the middle
my $largestScale = $largest*1000000;
my %sam2PlotStart = ();
   for my $eachSam(@sams){
       my $fragmentLen = $sam2Len{$eachSam};
       my $plotStart = ($largestScale/2-$fragmentLen/2)*$xratio;        
          $sam2PlotStart{$eachSam} = $plotStart;
   }

#--draw each chromosomes--------------------------
   my $chrXstart = $xstart;
   my $chrYstart = $ystart+30;
   my $chrnum = scalar(@sams); 
   #my $interval = ($height-50-50-100)/$chrnum;
   my $interval = 30;   

   my $flag = 0;
   for(my $m=0; $m<=$#sams; $m++){
       $eachSam = $sams[$m];
       $flag++;
       my $samLen = $sam2Len{$eachSam}; 

       my $samStartx = 0;
       my $samEndx   = $samLen;
          if(exists $sam2Start{$eachSam}){
             $samStartx = $sam2Start{$eachSam};
          } 
          if(exists  $sam2End{$eachSam}){
             $samEndx = $sam2End{$eachSam};
          }

        my $sf = sprintf "%.2f", ($samStartx/1000000);
        my $ef = sprintf "%.2f", ($samEndx/1000000); 
        my $sef = $sf."-".$ef." Mb";

        ##- set plot start
        my $plotStartTmp = $sam2PlotStart{$eachSam};

        my $high = 20;
           $svg -> rect(x => $chrXstart+$plotStartTmp, y => $chrYstart+$interval*(2*$flag-1), width => $samLen*$xratio, height => $high,               
                        stroke => 'black', fill => '#f6f5ec', "stroke-width"=>"0.5", 'stroke-opacity' => 0.2, 'fill-opacity' => 1);
           $svg -> text(x => $chrXstart-50+$plotStartTmp, y => ($chrYstart+$interval*(2*$flag-1)+$high+$chrYstart+$interval*2*$flag+$high)/2-25,
                       'font-family' => 'Arial', 'font-size' => 10, 'font-weight' => 'normal', stroke => 'black', 
                       'stroke-width' => 0.2, -cdata => $eachSam);
           
           $svg -> text(x => $chrXstart-70+$plotStartTmp, y => ($chrYstart+$interval*(2*$flag-1)+$high+$chrYstart+$interval*2*$flag+$high)/2-15,
                       'font-family' => 'Arial', 'font-size' => 8, 'font-weight' => 'normal', stroke => 'black', 
                       'stroke-width' => 0.2, -cdata => $sef);  ##- draw fragment region 
 


          ##- plot align
          my $fileName = $eachSam.".filtered.merged";    
          open INX, $fileName || die "cannot find $fileName";
          while(<INX>){
            chomp;
            my @tmp       = split(/\t/, $_);
            
             if($samInfo[-1] eq $tmp[0]){
            if($tmp[3] eq "CFv4" || $tmp[3] eq "ECD04" || $tmp[3] eq "A03" || $tmp[3] eq  "PCE" || $tmp[3] eq  "Z1v2"){
              if($tmp[1] >= $samStartx && $tmp[2] <= $samEndx){  
         
                 if(exists $sam2color{$tmp[3]}){
                    my $colorname = $sam2color{$tmp[3]};

                    ##- manually give color to Chiifu and ECD04/OIB haplotypes                   
                    #my $colorname = "#F8766D";    ##- "#F8766D", "#00BFC4" 
                    #if($tmp[3] eq "ECD04"){
                    #   $colorname = "#00BFC4";
                    #}
                    #if($tmp[3] eq "OIB"){
                    #   $colorname = "#B79F00";
                    #}

                    $svg -> rect(x => $chrXstart+($tmp[1]-$samStartx)*$xratio+$plotStartTmp, y => $chrYstart+$interval*(2*$flag-1), width => ($tmp[2]-$tmp[1]+1)*$xratio, height => $high,
                             stroke => 'black', fill => $colorname, "stroke-width"=>"0", 'stroke-opacity' => 0, 'fill-opacity' => 1);
                 }
                 else{
                     print "cannot find the color captured by the sample: $tmp[3] in file $fileName", "\n";
                 }   
              }
            }
           }
          }
          close INX;
         
        
          ##- plot synteny links
          if(($m+1) <= $#sams){
              my $ref   = $sams[$m];
              my $query = $sams[$m+1];
              my $file  = $coordDir.$ref."_".$query.".filtered.coords";
          
              ##- record plot start
              my $rplotStart = $sam2PlotStart{$ref};
              my $qplotStart = $sam2PlotStart{$query};

              my $rS = 0;
              my $rE = $sam2Len{$ref};
              my $qS = 0;
              my $qE = $sam2Len{$query};

              if(exists $sam2Start{$ref}){
                 $rS = $sam2Start{$ref};
              }
              if(exists $sam2End{$ref}) {
                 $rE = $sam2End{$ref};
              }        

              if(exists $sam2Start{$query}){
                 $qS = $sam2Start{$query};
              }
              if(exists $sam2End{$query}){
                 $qE = $sam2End{$query};
              }

              open IN2X, $file || die "cannot find  $file ";
              while(<IN2X>){
                chomp;
                my @temp = split(/\t/, $_);

                next, if($temp[11] ne $samInfo[-1]);
                next, if($temp[12] ne $samInfo[-1]);

                next, if($temp[0] < $rS);
                next, if($temp[1] > $rE);
                
                next, if($temp[2] < $qS);
                next, if($temp[3] > $qE);
 
 
                my ($refStart, $refEnd)     = ($temp[0]-$rS, $temp[1]-$rS);
                my ($QueryStart, $QueryEnd) = ($temp[2]-$qS, $temp[3]-$qS);
           
                my $x1 = $chrXstart+$refStart*$xratio+$rplotStart;
                my $x2 = $chrXstart+$refEnd*$xratio+$rplotStart;
                my $x3 = $chrXstart+$QueryStart*$xratio+$qplotStart;
                my $x4 =  $chrXstart+$QueryEnd*$xratio+$qplotStart;

                my $y1 = $chrYstart+$interval*(2*$flag-1)+$high+4;
                my $y2 = $chrYstart+$interval*(2*$flag-1)+$high+4;
                my $y3 = $chrYstart+$interval*(2*($flag+1)-1)-4;
                my $y4 = $chrYstart+$interval*(2*($flag+1)-1)-4;

                my $xv = [$x1, $x2, $x4, $x3];
                my $yv = [$y1, $y2, $y4, $y3];       
 
=pop
             my $xv = [$chrXstart+$refStart*$xratio, $chrXstart+$refEnd*$xratio, $chrXstart+$QueryEnd*$xratio, $chrXstart+$QueryStart*$xratio];
             my $yv = [$chrYstart+$interval*(2*$flag-1)+$high+4, $chrYstart+$interval*(2*$flag-1)+$high+4, $chrYstart+$interval*(2*($flag+1)-1)-4, $chrYstart+$interval*(2*($flag+1)-1)-4];
            
             my $points = $svg->get_path(                                                                                                                                       
                           x     =>  $xv,                                                                                                                                                  
                           y     =>  $yv,                                                                                                                                                  
                           -type =>'polygon'                                                                                                                                               
                        );

             $svg -> polygon(%$points, 'fill' => '#d3d7d4', "stroke-width"=>"0.1", "stroke-opacity" => "0.2", "fill-opacity" => "0.8") if($refStart < $refEnd && $QueryStart < $QueryEnd);  ##- gray 
             $svg -> polygon(%$points, 'fill' => '#00EE00', "stroke-width"=>"0.1", "stroke-opacity" => "0.2", "fill-opacity" => "0.8") if($refStart > $refEnd || $QueryStart > $QueryEnd);  ##- green

=cut
 
                my $x24= $x2;
                my $y24=($y2+$y4)/2;
                my $x42=$x4;
                my $y42=$y24;
                my $x31=$x3;
                my $y31=($y1+$y3)/2;
                my $x13=$x1;
                my $y13=($y1+$y3)/2;
 
             
                my $path = "M$x1,$y1 L$x2,$y2 C$x24,$y24 $x42,$y42 $x4,$y4 L$x3,$y3 C$x31,$y31 $x13,$y13 $x1,$y1 Z";

                ##- gray
                if($refStart < $refEnd && $QueryStart < $QueryEnd){
                   $svg->path(
                             d => $path,
                             fill=> "#d3d7d4",
                             stroke=>"#d3d7d4",
                             'stroke-width'=>'0.1',
                             'stroke-opacity'=>'0.2',
                             'fill-opacity'=>'0.8',
                             );
 
                } 
 
                ##- green
                if($refStart > $refEnd || $QueryStart > $QueryEnd){
                   $svg->path(
                             d => $path,
                             fill=> "#00EE00",
                             stroke=>"#00EE00",
                             'stroke-width'=>'0.5',
                             'stroke-opacity'=>'0.2',
                             'fill-opacity'=>'0.8',
                             );

                }

              }
              close IN2X;
       }    

   }

open OUT,">$out.svg";
print OUT $svg->xmlify;
system ("/10t/caix/bin/distributing_svg_4.74/svg2xxx_release/svg2xxx   $out.svg -t png -dpi 500");



#--subs------------------------------------------
sub getLongestChr {

    my ($FileIn, $samInfo, $sams, $sam2Len, $sam2startPos, $sam2endPos) = @_;

    my $tmp = $FileIn.".tmp";
    my $tmpchr = "";
 
    system("awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$4-\$3}'  $FileIn   | sort -k5,5nr  >  $tmp");   

    my $flag = 0;
    my $longest = 0;
    open IN1, $tmp;
    while(<IN1>){
      chomp;
      my @temp = split(/\t/, $_);
         $flag += 1;
         $longest = $temp[4], if($flag == 1); 
         $tmpchr  = $temp[1];
    }
    close IN1; 
    system("rm $tmp");  
      
    push(@{$samInfo}, $longest, $flag, $tmpchr);

    open INs, $FileIn;
    while(<INs>){
      chomp;
      my @temp = split(/\t/, $_);
      $sam2Len ->{$temp[0]} = $temp[3]-$temp[2]; 
      $sam2startPos ->{$temp[0]} = $temp[2];
      $sam2endPos ->{$temp[0]} = $temp[3];
      push(@{sams}, $temp[0]);
    }
    close INs;

}


sub colorSet {

    my @defaultColors = ("#F8766D", "#00BFC4", "#B79F00", "#619CFF", "#00BA38", "#F564E3", "#694d9f", "#f36c21", "#007d65", "#fdb933", "#ADFF2F", "#BC8F8F", "#A0522D", "#DEB887", "#BDB76B", "#F08080", "#5F9EA0", "#EE82EE", "#9370DB", "#DB7093");
    return(@defaultColors); 

}

sub sam2color {

    my ($sam, $color, $sam2color) = @_;
   
    my $count = scalar(@{$sam});
    my $colorCount = scalar(@{$color});
    for(my $m=0; $m<=($count-1); $m++){
        if($m<=($colorCount -1)){
           my $colorName = $color ->[$m];
              $sam2color ->{$sam ->[$m]} = $colorName;
        }
        else{
           my $colorName = $color ->[$colorCount-1]; 
              $sam2color ->{$sam ->[$m]} = $colorName;   
        }
    }

}
