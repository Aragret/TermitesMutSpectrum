#!/usr/bin/perl -w

%H=();
$h=$s="";
open (I, "$ARGV[0]");
while (<I>)
{
chomp;
if ($_=~/^>/) 
    {
    $t=$_;
    if ($s ne "")
	{
	$s=~s/\s+//g;
	$s=~tr/a-z/A-Z/;
	$H{$s}=$h;
	$s="";
	}
    $h=$t;
    }
else
    {
    $s.=$_;
    }
}
close I;

$H{$s}=$h;

foreach $k (keys %H)
{
print "$H{$k}\n$k\n";
}
