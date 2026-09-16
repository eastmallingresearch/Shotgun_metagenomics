#!/usr/bin/perl -s -w
use List::Util 'sum';
my %taxon_hash; 
while(<STDIN>) {
  chomp;
  my @array=split /,/,$_;
  foreach(@array) {
    $taxon_hash{$_}++;
  }
}

my $taxon_sum = sum values %taxon_hash;

foreach (keys %taxon_hash) {
  print "$_\t$taxon_hash{$_}\n";
}