#! /usr/bin/perl
# only works for csf's right now

if($#ARGV != 1)
{
  print "Usage: script filename cutoff \n";
  exit(1)
}

$filename = $ARGV[0];
$dbl1 = $ARGV[1];
# check file exists
open(INN,"< $filename");
@alldata = <INN>;
close(INN);

$dcnt = 0;
$int1 = 0;
$ndet = 0;

while($dcnt < $#alldata+1)
{
 if($alldata[$dcnt] =~ /<csf/)
 {
   $alldata[$dcnt] =~ /qchem_coeff.\"(\S*)\"/;
   $somevar = $1;
   #print "my num $somevar\n";
   if(abs($somevar) > abs($dbl1))
   { 
    $ndet++;
    $alldata[$dcnt] =~ /occ.\"(.*)[12](0*)"/;
    $mylength = length($1)+1;
    #print "dol1 $1 dol2 $2 \n";
    #print "my length is $mylength\n";
    if($mylength > $int1) 
    {
      $int1 = $mylength;
    } 
   }
 } 
 $dcnt++;
}

# must add number of core orbitals to this
print "Number of determinants: $ndet \n"; 
print "Required number of orbitals: $int1 \n"; 
