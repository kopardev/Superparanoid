use strict;
use DBI;
use lib qw(/usr/global/blp/perllib);
use Vutil;
use Bio::PrimarySeq;
#use Data::Dumper;
use Bio::SeqIO;

sub usage;

usage unless (@ARGV==2 or @ARGV==3 or @ARGV==4);

#INPUT
my $dbfile=shift;
#OUTPUT
my $geneidList=shift;
my $geneid_w_sequences_fasta=shift;
my $outputSpeciesName=shift; # only output SCOs from this species
chomp $outputSpeciesName;
my $dbh = DBI->connect(          
    "dbi:SQLite:dbname=$dbfile", 
    "",                          
    "",                          
    { RaiseError => 1 },         
) or die $DBI::errstr;

# get the species

my $sql_cmd="SELECT sname FROM Species ORDER BY sid";
my $sth = $dbh->prepare( $sql_cmd );
my $rv = $sth->execute() or die $DBI::errstr;
if($rv < 0){
   print $DBI::errstr;
}
my %sid2sname;
my $sid=0;
my $outputSpeciesFound=0;
my @fastas;
while (my @row = $sth->fetchrow_array()) {
	$sid+=1;
	my $speciesName="Species_".$sid;
	my $sname=$row[0];
	chomp $sname;
	push @fastas,$sname;
	Vutil::fileCheck($sname,"Fasta required to extract sequences by seqid!") if (defined $geneid_w_sequences_fasta) ;
	$sid2sname{$speciesName}=$sname;
	$outputSpeciesFound=1 if $sname eq $outputSpeciesName;
}
$sth->finish;
exit "No such species found=$outputSpeciesName\n!" if ($outputSpeciesFound==0 and defined $outputSpeciesName);

my %allsequences;
if (defined $geneid_w_sequences_fasta) {
	for my $fasta (@fastas){
		next if (defined $outputSpeciesName and $outputSpeciesName ne $fasta);
		$allsequences{$fasta}={};
		my $f=Bio::SeqIO->new(-file=>$fasta,-format=>'fasta');
		while (my $s=$f->next_seq){
			$allsequences{$fasta}{$s->id}=$s;
		}
	}	
}

my $nfastas = scalar @fastas;

my $sql_cmd="SELECT * FROM Clusters WHERE";
my $nwheres=0;
for (my $i=1;$i<$nfastas+1;$i++) {
	$sql_cmd.=" Species_".$i."=1" if $nwheres==0;
	$sql_cmd.=" AND Species_".$i."=1" if $nwheres!=0;
	$nwheres+=1;
}
my $sth = $dbh->prepare( $sql_cmd );
my $rv = $sth->execute() or die $DBI::errstr;
if($rv < 0){
   print $DBI::errstr;
}
my @cids;
open ID, ">$geneidList";
my $fo=Bio::SeqIO->new(-file=>">$geneid_w_sequences_fasta",-format=>'fasta') if (defined $geneid_w_sequences_fasta);
while (my @row = $sth->fetchrow_array()) {
	my $cid=shift @row;
	chomp $cid;
	push @cids,$cid;
	for my $x (1..$nfastas) {
		shift @row;
	}
	for my $i (1..$nfastas) {
		my @gids=split/:/,$row[$i-1];
		next if @gids == 0;
		my $sname=$sid2sname{"Species_$i"};
		for my $gid (@gids) {
			if (defined $geneid_w_sequences_fasta) {
				next if ((defined $outputSpeciesName) and ($outputSpeciesName ne $sname));
				my $id=sprintf "%s|%s|%s",$gid,$sname,$cid;
				my $seq=$allsequences{$sname}{$gid};
				my $s=Bio::PrimarySeq->new(-id=>$id,-seq=>$seq->seq,-desc=>$seq->desc);
				$fo->write_seq($s);
				my $desc=$s->desc;
				print ID "$cid\t$sname\t$gid\t$desc\n";
			} else {
				print ID "$cid\t$sname\t$gid\n";
			}
		}
	}
}
$sth->finish;
$dbh->disconnect();
my $nscos=scalar @cids;
print "$nscos single copy orthologs found.\n";
exit;


sub usage {
	print "Usage: perl $0 <db file> <output file with geneids> <output fasta with sequences(optional)> <species Name if only one species should be in the output(optional)>\n";
	exit;
}
