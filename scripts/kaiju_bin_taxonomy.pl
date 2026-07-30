#!/usr/bin/perl -s -w
use autodie;
use lib '/home/deakig/MyPerlLib/lib/perl5';

open (my $kf, '<', shift);# or die "Need kaiju file\n");
open (my $outfile, '>', shift);# or die "Need output file";
my $PROJECT_FOLDER=shift;

my $bin  = <$kf>;

my @splitline = split /\t/,$bin;

my $binname = (split /\./,$splitline[1])[1];


while (my $kline = <$kf>) {
	@splitline = split /\t/,$kline;
	my $newbinname = (split /\./,$splitline[1])[1];
	if($binname!=$newbinname) {
		open(my $t,">","temp.out");
		print $t $bin;
		my $taxon = parse_bin();
		print $outfile "bin.$binname\t$taxon";
		$bin="";
	} else {
		$bin = $bin.$kline;
	}
	
	$binname = $newbinname;
}

sub parse_bin {
	my $p = kaiju2table("phylum");
	my $c = kaiju2table("class");
	my $o = kaiju2table("order");
	my $f = kaiju2table("family");
	my $g = kaiju2table("genus");
	my $s = kaiju2table("species");

	return "$p\t$c\t$o\t$f\t$g\t$s\n"; 

}

sub kaiju2table {

	# set first argument to rank
	my ($rank) = @_;

	# system call to kaiju2table
	my $x = system("kaiju2table -t $PROJECT_FOLDER/data/kaiju/nodes.dmp -n $PROJECT_FOLDER/data/kaiju/names.dmp	-c 2 -r $rank -o temp.tsv temp.out");
	
	open (my $f, '<', "temp.tsv");
	my $header = <$f>;
	my @tophit = split /\t/, <$f>;
	my @tophit2 = split /\t/, <$f>; 
	close $f;
	chomp @tophit;
	
	my $tot = $tophit[2];#sprintf("%.0f",($tophit[2]*100)/$tophit[1]);
	my $pval = sprintf("%.2f",$tophit[1]);
	my $hitname = $tophit[4];
	$hitname="$hitname,$tophit2[4]" if $tot==$tophit2[2];
	
	return "$hitname\t$pval\t$tot";

}

open(my $t,">","temp.out");
print $t $bin;
my $taxon = parse_bin();
print $outfile "bin.$binname\t$taxon";
close $t;
close $outfile;
close $kf;


