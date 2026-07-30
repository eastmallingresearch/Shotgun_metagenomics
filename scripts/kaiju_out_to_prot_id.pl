#!/usr/bin/perl -s -w
use List::Util 'sum';
my %taxon_hash; 
while(<>) {
  chomp;
  my @array=split /\t/,$_;
  shift(@array);
  #print "$array[0]";
  
  my $al = scalar @array ;
  #print "$al\t";
  if( $al > 3) {  #- edited 15/01/24 to prevent warnings
  #if( $array[0] eq "C") {
	my @prot=split /,/,$array[4];
	
	$taxon_hash{$prot[0]}++;
#  pop(@array);
#  foreach(@array) {
 #   $taxon_hash{$_}++;
 # }
  }
}

my $taxon_sum = sum values %taxon_hash;

foreach (keys %taxon_hash) {
  print "$_\t$taxon_hash{$_}\n";
}

