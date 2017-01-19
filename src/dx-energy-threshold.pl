#!/usr/bin/perl 

# File: dx-combine.pl
# File: densityfinder_dxcontour_userinput.pl
# Original Script by Steven Ramsey
# Heavely modifed by Trent Balius, in B. Shoichet Lab, 2014
# Purpose: To find the high density water locations and print in a dx file the values (entropy/energy/whatever) to be contoured as wanted. 
# Will print to a single column .dat file to convert easily back to dx and also translate directly to dx format.
# 
#specifically uses user input to determine everything basically



###################################################
# This function reads in a dx file
# input: file name
#
# output: a vector containing the dx values. 
###################################################
sub read_dx_file {
  print "in side sub rotine read_dx_file\n";
  my($file, $IN, $OUT, $name, @values); # this local varible. 
  #my($file, $name, @values); # this local varible. 
  my($tempvalues, $line); # this local varible. 
  my($countx,$county,$countz,$originx,$originy,$originz,$incrementx,$incrementy,$incrementz);

  $file = $_[0] ;# this is argv 
  $name = $_[1] ;# this is argv 

  print "file = $file \n";

  open (IN, "< $file" || die "cannot open $file1 for data extraction \n");
  
  $line = <IN>; #first line is usually a comment
  print substr( $line, 0 , 1 ) . "\n"; # this is the first char of the string
  if ((substr( $line, 0 , 1 )) eq "#"){ # if comment read in next line.
     print $line . "\n";
     $line = <IN>; #line 2 contains the count for each variable: "object 1 class gridpositions counts ## ## ##"
  }
  # else use the first line 
  @count = split(/ /, $line);
 
  print $line, "\n"; 
  $countx = $count[5]; 
  $county = $count[6];
  #$countz =  $count[7];
  #$countz =  ($count[7] =~ s/\n//g);
  $countz = (substr( $count[7] , 0, 2 )); # removes last char which is the new line
  #$countz = (substr( $countz, 0, 3 )); #$count[7];

  #print $countz;
  #print ($countx =~ s/\n//g);
  #print ($countz =~ s/\n//g);
  #print (substr( $countz, 0 , $#countz )); # removes last char which is the new line

  print $countx,$county,$countz,"\n"; 
  
  $line = <IN>; #line 3 contains origins for each variable: "origin ## ## ##"
  @origin = split(/\s+/, $line);
  print $line, "\n";
  
  $originx = $origin[1];
  $originy = $origin[2];
  $originz = $origin[3];
  
  
  $line = <IN>; #line 4 contains increment value for first variable (x): "delta  (#) # #"
  $incrementx = substr($line, 7, 3);
  
  $line = <IN>; #line 5 contains increment value for second variable (y): "delta  # (#) #"
  $incrementy = substr($line, 9, 3);
  
  $line = <IN>; #line 6 contains increment value for third variable (z): "delta  # # (#)"
  $incrementz = substr($line, 11, 3);
  
  
  $line = <IN>; #line 7: object2
  $line = <IN>; #line 8: object3
  # $line = <IN>; #line 9: sometimes there is a blank line
  #
  $count = 0;
  @values = ();
  open (OUT, "> $name.dat" || die "Cannot create $name.dat to write data to" );  #printing all the data to a .dat file of a single column, this intermediate step may reduce time necessary to read in and evaluate the data.
  while (<IN>) {
          $line = $_;
          #print $line;
          ($tempvalues = $line) =~ s/\s+/\n/g;
          #@tmpvec = split(/\n/,$tempvalues),"\n";
          @values = (@values, split(/\n/,$tempvalues));
          #print $#values,"\n" ; 
          #@values = (@values, split($tempvalues));
          print OUT $tempvalues,"\n";
          print OUT $values[$count+0]," ",$values[$count+1]," ",$values[$count+2],"\n";
          #print $tempvalues,"\n";
          #print $values[$count+0]," ",$values[$count+1]," ",$values[$count+2],"\n";
          #print $tempvalues,"\n";
          $count = $count + 3;
  }

  print $#values,"\n" ; 

  close (IN);
  close (OUT);

  $total = int($countx)*int($county)*int($countz);
  #print $count," ", $countx," ", $county," ", $countz," ", $total,"\n";
  print $countx,"\n";
  print $county,"\n";
  print $countz,"\n";
  print $originx,"\n";
  print $originy,"\n";
  print $originz,"\n";
  print $incrementx,"\n";
  print $incrementy,"\n";
  print $incrementz,"\n";
  print $total,"\n";

  return ($countx,$county,$countz,$originx,$originy,$originz,$incrementx,$incrementy,$incrementz,$total,@values);# return value
} # end Subroutine read_dx_file

######################################################################
# this function writes a dx file. 
#
######################################################################
sub write_dx_file {
  print "in side sub rotine write_dx_file\n";
  my($file, $IN, $OUT, $name, @values); # this local varible. 
  #my($file, $name, @values); # this local varible. 
  my($countx, $county, $countz);
  my($originx,$originy,$originz); # this local varible. 
  my($incrementx,$incrementy,$incrementz);
  my($total);
  my(@values);
  # input arguments
  $file       = $_[0]; # this is argv 
  $countx     = $_[1];
  $county     = $_[2];
  $countz     = $_[3];
  $originx    = $_[4];
  $originy    = $_[5];
  $originz    = $_[6];
  $incrementx = $_[7];
  $incrementy = $_[8];
  $incrementz = $_[9];
  $total      = $_[10];
  @values     = @{$_[11]};

  print $file,"\n";
  print $countx,"\n";
  print $county,"\n";
  print $countz,"\n";
  print $originx,"\n";
  print $originy,"\n";
  print $originz,"\n";
  print $incrementx,"\n";
  print $incrementy,"\n";
  print $incrementz,"\n";
  print $total,"\n";

  #start dx format print modelled by code from Crystal Nguyen
  open (DX, ">$file" || "Can't open dx file to write into" );
  
  #print DX "#water thermodynamics on a grid\n"; # dock does like this coment line write now. 
  print DX "object 1 class gridpositions counts $countx $county $countz\n";
  print DX "origin $originx $originy $originz\n";
  #print DX "delta $incrementx 0 0\n";
  #print DX "delta 0 $incrementy 0\n";
  #print DX "delta 0 0 $incrementz\n";
  printf DX "delta %3.1f 0 0\n",$incrementx;
  printf DX "delta 0 %3.1f 0\n",$incrementy;
  printf DX "delta 0 0 %3.1f\n",$incrementz;
  print DX "object 2 class gridconnections counts $countx $county $countz\n";
  #print DX "object 3 class array type double rank 0 items $total data follows\n\n";
  print DX "object 3 class array type double rank 0 items $total data follows\n";
  
  
  $a = 0;   
  # We can replace this with a single loop. TEB 2014/04/02
  #$count = 0;   
  #for ($i = 0; $i < $countx; $i++ ) {
  #        for ($j = 0; $j < $county; $j++ ) {
  #        	for ($k = 0; $k < $countz; $k++ ) {
  #        		if ($a == 2) {
  #        			#print ("$values[$count-2] $values[$count-1] $values[$count]\n");
  #        			#print DX ("$values[$count-2] $values[$count-1] $values[$count]\n");
  #        			#printf DX "%6.3f %6.3f %6.3f\n", $values[$count-2] $values[$count-1] $values[$count];
  #        			printf DX "%g %g %g\n", $values[$count-2], $values[$count-1], $values[$count];
  #        			$a = -1;
  #        		}
  #        		$a++;
  #        		$count++;
  #        	}
  #        	$k = 0;
  #        }
  #        $j = 0;
  #}
  #if ($a > 0 && $a <= 2) {
  #        print "a is $a\n count = $count\n";
  #        print "size of values = ", $#values, "\n";
  #        for ($p = $a; $p > 0; $p--) {
  #              print (" p=",$p," index = ",($count-$p-1)," value=",$values[($count-$p-1)],"\n");
  #        	print DX ($values[($count-$p-1)]);
  #        }
  #        print DX "\n";
  #}

  $count_xyz = $countx*$county*$countz;
  for ($i = 0; $i < $count_xyz; $i++ ) {
      if ($a == 2) {
          printf DX "%g %g %g\n", $values[$i-2], $values[$i-1], $values[$i];
          $a = -1;
      }
      $a++;
  }

  if ($a > 0 && $a <= 2) {
          print "a is $a\n count = $count_xyz\n";
          print "size of values = ", $#values, "\n";
          for ($p = $a; $p > 0; $p--) {
                print (" p=",$p," index = ",($count_xyz-$p-1)," value=",$values[($countxyz-$p-1)],"\n");
          	print DX ($values[($count_xyz-$p-1)]);
          }
          print DX "\n";
  }
} # end Subroutine write_dx_file

######################################################################
# this function combines values file from the dx files. 
# function v[i] = w1*v1[i] + w2*v2[i]
######################################################################
sub combine_values {
  print "in side sub rotine combine values\n";
  my(@values1,$weight1); # this local varible. 
  my(@values2,$weight2);
  my(@values);
  # input arguments
  @values1    = @{$_[0]};
  $weight1    = $_[1];
  @values2    = @{$_[2]};
  $weight2    = $_[3];

  if ($#values1 != $#values2) {
     print "dx files do not have the same number of values\n";
     exit(-1);
  }

  print "size of values = $#values1 \n";

  for ($i = 0; $i < $#values1; $i++ ) {
       $values[$i] = $weight1 * $values1[$i] + $weight2 * $values2[$i];
       print "$values[$i] = $weight1 * $values1[$i] + $weight2 * $values2[$i] \n";
  }
  return (@values);# return value

} # end Subroutine combine_values

######################################################################
# this function reads in energy values and densities as vectors.
# output vn[i] = v[i], if d[i] > cutoff
#     or vn[i] = 0 otherwise
######################################################################
sub apply_density_threshold_values {
  print "in side sub rotine combine values\n";
  my(@values1); # this local varible.
  my(@densities);
  my(@values);
  my($threshold);
  # input arguments
  @values1    = @{$_[0]};
  @densities  = @{$_[1]};
  $threshold  = $_[2];

  if ($#values1 != $#densities) {
     print "dx files do not have the same number of values\n";
     exit(-1);
  }

  print "size of values = $#values1 \n";

  for ($i = 0; $i < $#values1; $i++ ) {
       if ( $densities[$i] > $threshold ){
          $values[$i] = $values1[$i]
       }
       else {
          $values[$i] = 0.0  
       }
       print "$values[$i] --- $values1[$i] , $densities[$i] , $threshold \n";
  }
  return (@values);# return value

} # end Subroutine combine_values

######################################################################
# this function apply a threshold to the dx file. 
# This will remove point that are close to zero 
# note the use of aboslute values.
# output vn[i] = v[i], if |v[i]| > cutoff
#     or vn[i] = 0 otherwise
#
######################################################################
sub apply_threshold_values {
  print "in side sub rotine combine values\n";
  my(@values1); # this local varible.
  my(@values);
  my($threshold);
  # input arguments
  @values1    = @{$_[0]};
  $threshold  = $_[2];

  print "size of values = $#values1 \n";

  for ($i = 0; $i < $#values1; $i++ ) {
       if ( abs($values1[$i]) > $threshold ){
          $values[$i] = $values1[$i]
       }
       else {
          $values[$i] = 0.0  
       }
       print "$values[$i] ---> $values1[$i] , $threshold \n";
  }
  return (@values);# return value

} # end Subroutine combine_values



sub main {
  # difine as local varables
  my($energyfile1,$theshold,$name1,@values_cut,$fileout_cut);
  my($countx1,$county1,$countz1,$originx1,$originy1,$originz1,$incrementx1,$incrementy1,$incrementz1,$total1,@energies);

  print "runing dx-threshold.pl\n";
  if (@ARGV != 3) {
        print "There are not the right number of inputs.\n";
        print "Please enter [dx energy filename] [threshold] [output name prfix]\n";
        exit(-1);
  }
  $energyfile1   = $ARGV[0];
  $threshold = $ARGV[1];
  $name1   = $ARGV[2];

  print "input: $energyfile1, $densityfile2, $threshold, $name1 \n";

  $fileout_cut = $name1 . "_cut.dx";
  $fileout1 = $name1 . "_energy1.dx";

  ($countx1,$county1,$countz1,$originx1,$originy1,$originz1,$incrementx1,$incrementy1,$incrementz1,$total1,@energies) = &read_dx_file($energyfile1,$fileout1);

  (@values_cut) = &apply_threshold_values(\@energies, $threshold);

  &write_dx_file($fileout_cut,$countx1,$county1,$countz1,$originx1,$originy1,$originz1,$incrementx1,$incrementy1,$incrementz1,$total1,\@values_cut);
} #end sub main

# call function main
&main

