#!/usr/bin/perl -s -w
use autodie;
# use Number::Range;
use lib '/home/deakig/MyPerlLib/lib/perl5';
use Set::IntervalTree;

my (%ranges_hash,%objects_hash) = create_hash_from_gff();


while (my $sam_line = <STDIN>) {
	my @proccesedSam = process_sam($sam_line);
	if (exists $ranges_hash{$proccesedSam[0]}) {
		my $results_ref = $ranges_hash{$proccesedSam[0]}->fetch($proccesedSam[1],$proccesedSam[2]);
		foreach (@{$results_ref}) {
			#print"$proccesedSam[0], $_, $proccesedSam[1],  $proccesedSam[2],$proccesedSam[3],$proccesedSam[4]\n";   
			$objects_hash{$_}++;
		}
	}
}

foreach (keys %objects_hash) {
	print "$_ $objects_hash{$_}\n" if  $objects_hash{$_}>0;
}

sub process_sam {
	my @line = split /\t/,$_[0];
	my $x = $line[2];
	($line[2])= split / /,$line[2];
	return($line[2],$line[3],length($line[9])+$line[3]-1); #,$line[0],$x);
}

sub create_hash_from_gff {
	# this reads the gff file three times - it's faster than checking if keys already exist

	# get unique values from first column of gff file 
	open my $fh, "cut -f1 $ARGV[0]|sort|uniq|";

	# read column into hash and chomp
	my %gff_hash=map{chomp; $_ => 1 }<$fh>;

	# create an empty interval tree for each hash key
	foreach (keys %gff_hash) {
		$gff_hash{$_}=Set::IntervalTree->new;
	}

	close $fh;
	
	# get unique values from ID column of gff file 
	open $fh, "cut -f9 $ARGV[0]|sort|uniq|";

	# read column into hash and chomp
	#my %id_hash=map{chomp; substr($_,3) => 1 }<$fh>;
	my %id_hash=map{s/ID=//;$_}map{chomp; $_ => 1}<$fh>;

	# set value to 0 for each ID key
	foreach (keys %id_hash) {
		$id_hash{$_}=0;
		#print"$_\n";
	}

	close $fh;	

	# reopen gff file
	open $fh, '<', $ARGV[0];
	my @gff = <$fh>;
	chomp(@gff);
	#print(pop @gff);
	shift @gff;

	# go through the gff and for each unique hash key assign each gff ID and range to value
	foreach(@gff) {
		my @line = split /\t/, $_;
		$line[8]=~s/ID=//;
		$gff_hash{$line[0]}->insert($line[8],$line[3],$line[4])
	}
	close $fh;
	return (%gff_hash,%id_hash);
}

# debug
# foreach (keys %hash) {
#	my $results_ref = $hash{$_}->fetch(500,510);
#	print "$_  @{$results_ref}\n";
#}




