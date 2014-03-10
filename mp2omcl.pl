use strict;
use warnings;

my $inFile=shift;
my %clusters;
open IN, $inFile;
my $ln=0;
while (my $line=<IN>) {
$ln=$ln+1;
next if $line=~/^cluster/;
chomp $line;
my @cols=split /\t/,$line;
#my $clusterName="cluster".$cols[0];
my $clusterName=$cols[0];
$clusters{$clusterName}=() unless defined $clusters{$clusterName};
my $geneName=$cols[1]."##".$cols[2];
push @{$clusters{$clusterName}},$geneName;
#last if $ln==10;
}
for my $clusterName (keys %clusters) {
print "$clusterName";
foreach my $geneName (@{$clusters{$clusterName}}) {
	print "\t$geneName";
}
print "\n";
}

