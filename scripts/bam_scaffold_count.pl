#!/usr/bin/perl -s -w
use autodie;
use lib '/home/deakig/MyPerlLib/lib/perl5';
use Set::IntervalTree;
#use Parallel::ForkManager;

my $gff_out=0;

$gff_out=1 if scalar @ARGV>1;

# testing parallel processing
#my $MAX_PROCESS = shift;

#my $pm = new Parallel::ForkManager($MAX_PROCESS); 

#my $INP;

#foreach (@ARGV) {
#	$pm->start and next;
#	open $INP,'<',$_;
#	while (<$INP>) {
#		print"$_";
#	}
#	close $INP;
#	$pm->finish;
#}
#exit;

my (%ranges_hash,%objects_hash) = create_hash_from_gff($gff_out);
print STDERR ("Created hash from GFF\n");
my %pe;

while (my $sam_line = <STDIN>) {
	my @proccesedSam = process_sam($sam_line);
	if (exists $ranges_hash{$proccesedSam[0]}) {
		my $results_ref = $ranges_hash{$proccesedSam[0]}->fetch($proccesedSam[1],$proccesedSam[2]);
		# check if paired end - only count new matches
		if(exists $pe{$proccesedSam[3]}) {
			foreach (@{$results_ref}) {
				$objects_hash{$_}++ if !exists $pe{$proccesedSam[3]}{$_};
			}
			delete $pe{$proccesedSam[3]};	
		} else {
			#$pe{$proccesedSam[3]} = $results_ref;		
			foreach (@{$results_ref}) {
				#print"$proccesedSam[0], $_, $proccesedSam[1],  $proccesedSam[2],$proccesedSam[3],$proccesedSam[4]\n";   
				$objects_hash{$_}++;
				$pe{$proccesedSam[3]}{$_}=0
			}
		}
	}
}

foreach (keys %objects_hash) {
	print "$_\t$objects_hash{$_}\n" if  $objects_hash{$_}>0;
}

sub process_sam {
	my @line = split /\t/,$_[0];
	my $x = $line[2];
	($line[2])= split / /,$line[2];
	return($line[2],$line[3],length($line[9])+$line[3]-1),$line[0]; #,$line[0],$x);
}

sub create_hash_from_gff {
	
	# get first argument
	my ($gfo) = @_;

	# this reads the gff file three times - it's faster than checking if keys already exist

	# get unique values from first column of gff file 
	open my $fh, "cut -f1 $ARGV[0]|sort|uniq|";

	# read column into hash and chomp
	my %gff_hash=map{chomp; $_ => Set::IntervalTree->new }<$fh>;

	# create an empty interval tree for each hash key
#	foreach (keys %gff_hash) {
#		$gff_hash{$_}=Set::IntervalTree->new;
#	}
	close $fh;

	my %id_hash;

	if(!$gfo) {

		# get unique values from ID column of gff file (domain names)
		open $fh, "cut -f9 $ARGV[0]|sort|uniq|";

		# read column into hash and chomp
		%id_hash=map{s/ID=//;$_}map{chomp; $_ => 0}<$fh>;

		close $fh;
	}

	# reopen gff file
	open $fh, '<', $ARGV[0];
	my @gff = <$fh>;
	close $fh;
	chomp(@gff);
	shift @gff;
	
	# go through the gff and for each unique hash key assign each gff ID and range to value
	
	foreach(@gff) {
		#$id_hash{$_}=0;
		my @line = split /\t/, $_;
		$line[8]=~s/ID=//;
		$gff_hash{$line[0]}->insert($_,$line[3],$line[4])
	}
	%id_hash=map{$_ => 0}@gff if $gfo;

	return (%gff_hash,%id_hash);
}

# debug
# foreach (keys %hash) {
#	my $results_ref = $hash{$_}->fetch(500,510);
#	print "$_  @{$results_ref}\n";
#}




