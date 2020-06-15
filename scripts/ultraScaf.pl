#!/usr/bin/perl -w
use strict;

=head1 NAME

ultraScat.pl


=head1 VERSION

v1.0

=head1 PURPOSE

Used for generating ultra scaffolds/chromosomes from a fasta file containing scaffolds or contings

=head1 USAGE

ultraScat.pl <fasta file> <list> <prefix>

=head1 AUTHOR

Christina Toft, christina.toft@uv.es

=cut

use Bio::Perl;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

open(IN,$ARGV[0]) || die;

my($seq,$name);
my(%dataIN);
my($prefix) = $ARGV[2];

while(my $line = <IN>){
	if($line =~ /^>(\S+)/){
		$name = $1;
		#print "$name ";
	}else{
		chomp $line;
		$dataIN{$name} .= $line;
	}
}

close IN;

open(COORDS,$ARGV[1]) || die;
my(%dataOUT);
my($refChrm) = "";
my($contig) = "";
my($cov,$end1,$start1,$start2,$end2,$length),$name;
my($first) = 1;
my(@used);
my($ultraScaf);
while(my $line = <COORDS>){
	if($line =~ /^\d/){
		chomp $line;
		my(@array) = split(/\t/,$line);

		if($contig ne $array[12]){
			if("@used " !~ /\b$array[12]\s/){

				if($contig ne ""){
					if($start2 < $end2){
						print "$start1\t$end1\t$start2\t$end2\t$contig\n";
						$ultraScaf .= "$contig ";
					}else{
						print "$end1\t$end2\t$start2\t$end2\t$contig"."COM\n";
						# complement
						$ultraScaf .= $contig."COM ";
						my($temp) = $end2;
						$end2 = $start2;
						$start2 = $temp;
					}
				}
				my $temp = ($length - $end2 + $array[2] - 1);
				$temp = ($length - $end2 + $array[3] -1) if $array[3] < $array[2];
				#print "ref dif ($array[0] - $end1):" . ($array[0] - $end1) ."\t sub dif ($end2 $length $array[2]): " .($temp) . "\ndif " .($array[0] - $end1 - $temp)." \n\n";

				my($difference) = $array[0] - $end1 - $temp;
				if($refChrm eq $array[11]){
					if($difference < 1){
						$ultraScaf .= "0 " ;
					}else{
						$ultraScaf .= "$difference " ;
					}
				}
				$contig = $array[12];
				#print "$array[0]\t";
				push(@used,$contig);
				$start2 = $array[2];
				$length=$array[8];
			}else{
				print "Check $array[12]\n$line\n";
				die;
			}

		}
		$ultraScaf .= "\n$array[11] " if $refChrm ne $array[11];
		$refChrm = $array[11];
		$end1 = $array[1];
		$end2 = $array[3];
		#print "start:$start2\tend:$end2\n";
	}
	#exit if $refChrm eq "Sbay_02";
}
#$ultraScaf .= "$contig", ($start2 < $end2) ? "" : "COM", "\n";
$ultraScaf .= "$contig";
$ultraScaf .= "COM" if $start2 > $end2;
$ultraScaf .= "\n";
#print "start:$start2\tend:$end2\n$ultraScaf";
close COORDS;
#die;
#print "$ultraScaf";
my(@usedScaffolds);

foreach my $line(split(/\n/,$ultraScaf)){
	print ">$line<\n";
	next if $line eq "";
	my(@array) = split(/\s+/,$line) ;
	my(@feats);
	my($text);

	print "chrm $array[0]\n";
	my($seq);
	for(my $i = 1;$i < scalar(@array); $i++){
		if($array[$i] =~ /(\D+\d+\S*)/){
			if($array[$i] =~ /(\D+\d+\S*)COM/){
				print "$array[$i] is the reverse comp\n";
				$temp = reverse($dataIN{$1});
				$temp =~ tr/ACGTacgt/TGCAtgca/;
				$dataIN{$1} = $temp;
				$array[$i] = $1
			}
			#$array[$i] =~ /(\D+\d+)/;
			$temp = (length($seq)+1);
			#my($tempTwo) = (length($dataIN{$1}));
			$seq .= $dataIN{$array[$i]};
			my($tempTwo) = (length($seq));
			$text .= "$1 $temp $tempTwo\n";
			push(@usedScaffolds,$1);
			#print ".." .length($seq) . "\n";
			#print "scaf: " .length($dataIN{$1}) ."\n\n";

		}else{
			$seq .= "n";
			for(my $k = 2; $k <= $array[$i];$k++){
				$seq .= "n";
			}
		}
	}

	my($num) = $array[0];
	$num =~ s/\D+(\d+)/$1/; ############################################ commend out if hybrid
	print "fears:" . scalar(@fears) ."\n";
	my $seq_obj = Bio::Seq->new(-seq => "$seq",
                           -display_id => "$prefix\_$num" );

	print "$text\n";
	my $feat = new Bio::SeqFeature::Generic(-start       => 1,
                                        -end         => length($seq),
                                        -strand      => 1,
                                        -primary_tag => 'misc_feature',
                                        -tag => {locus_tag => $prefix."_$num",
                                                 note     => 'ultra scaffold' } );
	$feat->add_tag_value("colour", 2 );
	$seq_obj->add_SeqFeature($feat);

	foreach my $mis(split(/\n/,$text)){
		my(@misArray) = split(/\s/,$mis);
		my $feat = new Bio::SeqFeature::Generic(-start       => $misArray[1],
                                        -end         => $misArray[2],
                                        -strand      => 1,
                                        -primary_tag => 'misc_feature',
                                        -tag => {locus_tag => "$misArray[0]",
                                                 note     => 'scaffold' } );

		$seq_obj->add_SeqFeature($feat);
#	$seq_obj->add_SeqFeature(@feats);
	}


	my $io = Bio::SeqIO->new(-format => "genbank", -file => ">" . $prefix . "_$num.gb" );
	$io->write_seq($seq_obj);
	#die;
}

# ALL scaffolds not used....

$seq="";
$text = "";
my(@keys) = keys %dataIN;
my(@sorted) = sort { $a <=> $b } @keys;
foreach my $scaf(sort @sorted){
	next if "@usedScaffolds" =~ /\b$scaf\b/;
	print "$scaf\n";
	my $temp = (length($seq)+1);
	$seq .= $dataIN{$scaf};
	my($tempTwo) = (length($seq));
	$text .= "$scaf $temp $tempTwo\n";
}


my $seq_obj = Bio::Seq->new(-seq => "$seq",
                           -display_id => $prefix."_unplaced" );

foreach my $mis(split(/\n/,$text)){
	my(@misArray) = split(/\s/,$mis);
	my $feat = new Bio::SeqFeature::Generic(-start       => $misArray[1],
                                       -end         => $misArray[2],
                                       -strand      => 1,
                                       -primary_tag => 'misc_feature',
                                       -tag => {locus_tag => "$misArray[0]",
                                                note     => 'scaffold' } );
		$seq_obj->add_SeqFeature($feat);
	$seq_obj->add_SeqFeature(@feats);
}
my $io = Bio::SeqIO->new(-format => "genbank", -file => ">" . $prefix . "_unplaced.gb" );
$io->write_seq($seq_obj);

exit;
