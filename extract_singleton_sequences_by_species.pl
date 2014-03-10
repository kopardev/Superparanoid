use strict;
use DBI;
use lib qw(/usr/global/blp/perllib);
use Vutil;
use Bio::SeqIO;

sub usage;

usage unless (@ARGV==3 or @ARGV==4);

#INPUT
my $dbfile=shift;
my $speciesName=shift;
#OUTPUT
my $geneidList=shift;
my $geneid_w_sequences_fasta=shift;

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
my $ctr=0;
my $sid=0;
while (my @row = $sth->fetchrow_array()) {
	$ctr+=1;
	my $sname=$row[0];
	chomp $sname;
	$sid=$ctr if $sname eq $speciesName;
}
$sth->finish;
exit "Species $speciesName does not exist!\n" if $sid==0;

my $fi=Bio::SeqIO->new(-file=>$speciesName,-format=>'fasta');
my %allsequences;
while (my $s=$fi->next_seq) {
	$allsequences{$s->id}=$s;
}
my $fo=Bio::SeqIO->new(-file=>">$geneid_w_sequences_fasta",-format=>'fasta') if (defined $geneid_w_sequences_fasta);
open FO, ">$geneidList";

my $sql_cmd="select geneid from Singletons where sid = $sid";
my $sth = $dbh->prepare( $sql_cmd );
my $rv = $sth->execute() or die $DBI::errstr;
if($rv < 0){
   print $DBI::errstr;
}
while (my @row = $sth->fetchrow_array()) {
	my $geneid=$row[0];
	chomp $geneid;
	my $desc=$allsequences{$geneid}->desc;
	print FO "$geneid\t$desc\n";
	$fo->write_seq($allsequences{$geneid}) if (defined $geneid_w_sequences_fasta);
};
$sth->finish;
close FO;

$dbh->disconnect();
exit;


sub usage {
	print "Usage: perl $0 <db file> <species fasta> <output file with geneids> <output fasta with sequences(optional)>\n";
	exit;
}
