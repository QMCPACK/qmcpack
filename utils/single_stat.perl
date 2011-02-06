#!/usr/bin/env perl
#################################################
# evaluate the average, error, sigma and autocorellation time of multiple files
# Each file is an isolated data and no grouped average is taken
#################################################
# script is "./g"
use Getopt::Std;
%options=();
getopts("i:n:f:",\%options);

my $target_index=1;
my $target_name="LE";
my $first_data=0;
$target_index=$options{i} if defined $options{i};
$target_name=$options{n} if defined $options{n};
$first_data=$options{f} if defined $options{f};

#print "Select index $options{i}\n" if defined $options{i};

foreach my $fname (@ARGV) {
  &processAfile($fname);
}

sub processAfile {
  my @data_array;
  my($a)=@_;
  my $irow=0;
  open(MFILE,"$a");
  while(<MFILE>) {
    chop($_);
    @w=split(' ');
    if($w[0] ne "#" && $irow>$first_data) 
    {
      push(@data_array,$w[$target_index]); 
    }
    $irow++;
  }
  close(MFILE);

  my $lower_index=0;
  my $upper_index=@data_array-1;

  if($upper_index>5)
  {
    my @res = &calcAutoCorr($lower_index,$upper_index,@data_array);
    printf("%s %s %.6e %.1e %.6e %.6e\n",$a, $target_name, $res[0], $res[1], $res[2], $res[3]);
  }
}

###########################################
# below : Common subroutines
###########################################
#calculate mean and variance of a data 
sub calcMeanAndVariance {
  my(@data) = @_;
  my $sum=0.0;
  my $var=0.0;
  foreach my $val (@data) {
    $sum += $val;
    $var += $val*$val;
  }

  my $n = @data;
  my $avg = $sum/$n;
  return ($avg,$var/$n-$avg*$avg);
}

#calculate mean, error, sigma=sqrt(variance) and autocorrelation time  of a data
sub calcAutoCorr {
  my ($imin,$imax,@data) = @_;
  my ($mean,$var) = &calcMeanAndVariance(@data);
  my $auto_min=1;
  my $auto_max=1;
  my $n=$imax-$imin+1;
  my $cutoff = $n;
  my @auto;
  my @corrTime;
  for (my $j=0; $j<$cutoff; $j++) {
    push(@auto,0.0);
  }
  for (my $j=0; $j<$cutoff; $j++) {
    for (my $i=$imin; $i+$j<=$imax; $i++) {
      $auto[$j]+=($data[$i]-$mean)*($data[$i+$j]-$mean);
    }
    $auto[$j]=$auto[$j]/($var*($n-$j));
    if ( $auto[$j] < sqrt( 2.0/($n-$j))) {
      if($cutoff>5*$j) {
        $cutoff = 5*$j;
      }
    }
    if ($auto[$j]>$auto_max) {
      $auto_max=$auto[$j];
    }
    if ($auto[$j]<$auto_min) {
      $auto_min=$auto[$j];
    }
  }

  #remove data from an array
  my $autolen=@auto;
  if ($autolen!=$cutoff) {
    my $dn=$autolen-$cutoff;
    while($dn>0) {
      pop(@auto); 
      $dn--;
    }
    $autolen-=$dn;
  }
  if($autolen == 0) {
    $corrTime=1.0;
  } else {
    $corrTime=$auto[0];
    for (my $j=1; $j*5 < 2*$cutoff; $j++) {
      $corrTime += 2*$auto[$j];
    }
  }
  my $sigma = sqrt($var);
  my $err = $sigma*sqrt($corrTime/$n);
  return ($mean,$err,$sigma,$corrTime);
}
