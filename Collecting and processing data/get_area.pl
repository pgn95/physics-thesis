#!/usr/bin/perl
$mocfile=$ARGV[0];
open(MOC,"$mocfile") || die "Cannot open $mocfile\n";
while(<MOC>){
    chomp;
    if($_ =~ /\"/){
# Linea anterior	
	$npix=@pixels;
#	print "$order $npix\n";
	$Omega=$Omega+$npix*(($pix/3600)**2);
# Nueva linea	
	undef @pixels;	$npix=0;
	/\"(.*)\"\:\[(.*)\,/;
	$order=$1;
	@pixels=split(",",$2);
	$linepix=`grep \" $order \" MOCorder2pixel.topcat`;chomp($linepix);
	($order,$pix)=split(" ",$linepix);
    }else{
	/(.*)\,/;
	@data=split(",",$1);
	@pixels=(@pixels,@data);
    }
}
# Ultima linea
$npix=@pixels;
#print "$order $npix\n";
$Omega=$Omega+$npix*(($pix/3600)**2);
#print "Omega: $Omega\n";
print "$Omega\n";
