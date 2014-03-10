use strict;
use warnings;
use Bio::SeqIO;
use Bio::PrimarySeq;
my $inFile=shift;
my $outFile=shift;
my $in=new Bio::SeqIO(-file=>$inFile,-format=>"genbank");
my $out=new Bio::SeqIO(-file=>">$outFile",-format=>"fasta");
while (my $seq=$in->next_seq) {
	#$out->write_seq($seq);
for my $feat_object ($seq->get_SeqFeatures) {          
	next unless ($feat_object->primary_tag eq "CDS");
   my $start=$feat_object->start;
   my $end=$feat_object->end;
   my $contigname=$feat_object->spliced_seq->id;
   my $newSeqName="${contigname}_${start}_${end}";
   my $newSeq=$feat_object->spliced_seq->seq;
   #my $newSeq=$feat_object->get_tag_values("translation");
   my $newSeqObj=new Bio::PrimarySeq(-id=>$newSeqName,-seq=>$newSeq);
   my $newProtSeqObj=$newSeqObj->translate();
   my $newProtSeq=$newProtSeqObj->seq();
   $newProtSeq=~s/\*//g;
   my $anotherNewProtSeqObj=new Bio::PrimarySeq(-id=>$newSeqName,-seq=>$newProtSeq);
   #print $feat_object->get_tag_values("translation"),"\n";
   #for my $tag ($feat_object->get_all_tags) {            
    #  print "  tag: ", $tag, "\n";            
     # for my $value ($feat_object->get_tag_values($tag)) {
#
 #       print "    value: ", $value, "\n";            
  #  }          
   #}      
	#$out->write_seq($newSeqObj);
	#$out->write_seq($newProtSeqObj);
	$out->write_seq($anotherNewProtSeqObj);
}

}
