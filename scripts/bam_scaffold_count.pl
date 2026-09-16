#!/usr/bin/perl -s -w
use autodie;
use lib '/home/deakig/MyPerlLib/lib/perl5';
use Set::IntervalTree;
#use Parallel::ForkManager;

################################################################################
#																																							 #
# Counts feature overlaps between a gff file and headerless SAM input.				 #
#																																							 #
#	bam_scaffold_count.pl GFF_FILE [BEDTOOLS] < SAM_LINE			      						 #
#																																							 #
# IF BEDTOOLS = True, counts will be as per bedtools with features 					   #
# overlapping forward and reverse reads counted twice. Default counts once.    #
#																																							 #
################################################################################

# get command line variables (if set)
my $gff = shift or die "Usage: $0 FILENAME\n";
my $gff_out = 1; # shift // 0; Always BEDTOOLS cov output
my $bed_compatible = shift // 0;
#$gff_out=1 if scalar @ARGV>1;

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

# get refs to hashes
my($ranges_ref,$objects_ref) = &create_hash_from_gff($gff_out);

print STDERR ("Created hash from GFF\n");
my %pe;
print STDERR "BEDTOOLS counts - PE matches counted twice\n" if $bed_compatible;

# Process stdin (SAM)
while (my $sam_line = <STDIN>) {
	my @proccesedSam = process_sam($sam_line);
	if (exists $ranges_ref->{$proccesedSam[0]}) {
		my $results_ref = $ranges_ref->{$proccesedSam[0]}->fetch($proccesedSam[1],$proccesedSam[2]);
		# check if paired end - count twice if BEDTOOLS is true
		if(exists $pe{$proccesedSam[3]}) {
			foreach (@{$results_ref}) {
				if(!$bed_compatible) {
					$objects_ref->{$_}++ if !exists $pe{$proccesedSam[3]}{$_};
				} else {
					$objects_ref->{$_}++;
				}
			}
			delete $pe{$proccesedSam[3]};
		} else {
			#$pe{$proccesedSam[3]} = $results_ref;
			foreach (@{$results_ref}) {
				#print"$proccesedSam[0], $_, $proccesedSam[1],  $proccesedSam[2],$proccesedSam[3],$proccesedSam[4]\n";
				$objects_ref->{$_}++;
				$pe{$proccesedSam[3]}{$_}=0
			}
		}
	}
}

foreach (keys %$objects_ref) {
	if ($objects_ref->{$_}>0) {
		my $out = $objects_ref->{$_};
		print "$_\t$out\n"
	}
}

sub process_sam {
	# splits up tab seperated sam line and returns columns we're interested in
	my @line = split /\t/,$_[0];
	my $x = $line[2];
	($line[2])= split / /,$line[2];
	return($line[2],$line[3],(length($line[9])+$line[3]-1),$line[0]); #,$line[0],$x);
}

sub create_hash_from_gff {

	# get first argument
	my ($gfo) = @_;

  # reopen gff file
	open $fh, '<', $gff;
	my @gff = <$fh>;
	close $fh;
	chomp(@gff);
	shift @gff;

  my %gff_hash;
  my %id_hash;

	# go through the gff and for each unique hash key assign each gff ID and range to value
	foreach(@gff) {
		my @line = split /\t/, $_;
		#$line[8]=~s/ID=//;
    if (exists $gff_hash{$line[0]}) {
		    $gff_hash{$line[0]}->insert($_,$line[3],$line[4]);
    } else{
      $gff_hash{$line[0]}=Set::IntervalTree->new;
      $gff_hash{$line[0]}->insert($_,$line[3],$line[4]);
    }
    # $id_hash{$line[8]}=0 if (!exists $id_hash{$line[8]});
	}
  %id_hash=map{$_ => 0}@gff if $gfo;
	return (\%gff_hash,\%id_hash);
}
