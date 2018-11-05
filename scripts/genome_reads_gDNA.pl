#!/usr/local/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Config::Simple;

my $cfg = new Config::Simple('SCpipe.conf');
my %config = $cfg->vars();

$| = 1;

&GetOptions
(
	"c=s"=>\my $coords_file,
	"b=s"=> \my $bam_dir,
	"o=s"=> \my $out_dir
);
unless($coords_file) {die "ERROR: no coords file!\nUSAGE: genome_reads.pl -c full/path/coords.bed -b /full/path/bam/ -o full/path/output/\n\n";}
unless($bam_dir) {die "ERROR: no bam files dir!\ngenome_reads.pl -c full/path/coords.bed -b /full/path/bam/ -o full/path/output/\n\n";}
unless($out_dir) {die "ERROR: no output dir!\ngenome_reads.pl -c full/path/coords.bed -b /full/path/bam/ -o full/path/output/\n\n";}
unless($bam_dir =~ m/\/$/) {$bam_dir .= '/';}
unless($out_dir =~ m/\/$/) {$out_dir .= '/';}

#

my %coords;
open (COORDS, "<$coords_file") or die "cannot open $coords_file: $!\n";
while (<COORDS>)
{
	my $line = $_;
	chomp $line;
	unless ($line){next;}
	if ($line =~ m/^Chr/) {next;}
	
	my ($chr, $from, $to, $amplicon) = split(/\t/, $line);
	unless($chr =~ m/^chr/) {$chr = 'chr' . $chr;}
	
	$coords{$chr}{$from}{to} = $to;
	$coords{$chr}{$from}{amplicon} = $amplicon;
}
close COORDS;



opendir(DIR, $bam_dir);
my @files = readdir(DIR);
closedir(DIR);
 
foreach my $file (@files)
{
	if ($file eq "." or $file eq "..") {next;}
	unless ($file =~ m/bam$/) {next;}
	
	my $genome_sam_file = $file;
	$genome_sam_file =~ s/bam$/genome.sam/;
	
	$file = $bam_dir . $file;
	$genome_sam_file = $out_dir . $genome_sam_file;
	
	my $cmd = $config{SAMTOOLS_EXECUTABLE} . ' view ' . $file;
	my %sam_data = ();
	my %genome_read_ids = ();
	open (PIPE, "$cmd |") || die "$!";
	while (<PIPE>)
	{
		my $line = $_;  
		chomp $line;
		
		my ($read_id, $bitwise_flag, $ref_seq, $start, $qual, $cigar, $read_pair_ref_seq, $read_pair_start, $insert_size, $seq, $bp_scores, $optional_field) = split(/\t/, $line);
		$sam_data{$read_id}{$ref_seq}{$start} = $line;
		
		my ($bits_ref) = &parseSAMbitwise($bitwise_flag);
		my %bits = %{$bits_ref};
		
		# are there any skipped bases (ie spliced reads arising from mRNA)?
		my $skipped = 0;
		while ($cigar =~ m/(\d+)N/g)
		{
			$skipped += $1;
		}
		#print "Skipped: $skipped\n-------------\n";
		
		# check read is mapped in a proper pair:
		unless(exists($bits{b2})) {next;}
		#Eliminate any potential reads which are spliced: CIGAR FORMAT 'N' 
		unless($skipped==0) {next;}
		
		foreach my $from_coord (keys %{$coords{$ref_seq}})
		{
			my $to_coord = $coords{$ref_seq}{$from_coord}{to};
			my $mid_posn = (($to_coord - $from_coord) / 2);
			my $amplicon = $coords{$ref_seq}{$from_coord}{amplicon};

			if(exists($bits{b10}))
			{
				# reads on reverse strand:
				my $adjusted_start = ($start - $mid_posn) + length($seq) + $skipped;
				if (($adjusted_start >= $from_coord) and ($adjusted_start <= $to_coord))
				{
					$genome_read_ids{$read_id}{$ref_seq}{$start} = '';
				}
			}
			else
			{
				# reads on fwd strand:
				my $adjusted_start = ($start + $mid_posn);
				#my $adjusted_start = ($start);
				if (($adjusted_start >= $from_coord) and ($adjusted_start <= $to_coord))
				{
					$genome_read_ids{$read_id}{$ref_seq}{$start} = '';
				}
			}
		}
	}
	close PIPE;	

	my $count = 0;
	open (OUT, ">$genome_sam_file") or die "cannot open out file: $genome_sam_file: $!\n";
	foreach my $read_id (sort keys %genome_read_ids)
	{
		foreach my $read_chr (keys %{$genome_read_ids{$read_id}})
		{
			foreach my $read_coord (sort numerically keys %{$genome_read_ids{$read_id}{$read_chr}})
			{
				print OUT $sam_data{$read_id}{$read_chr}{$read_coord} . "\n";
				$count ++;
			}
		}
	}
	close OUT;
	
	print "Reads: $count\n";

}


#----------------------------------------------------------------------------------------

sub parseSAMbitwise
{
############################################
# Interpret bitwise flag form the SAM format (http://samtools.sourceforge.net/SAM1.pdf)
# Input:  flag
# Output: reference to a hash of set bits and associated values
############################################

	my ($flag) = @_;
	
	my %bits;
	
	if($flag == 0)
	{
		$bits{'b0'} = 'mapped to forward strand, not paired-end';
		return(\%bits);
	}
	
	if ($flag & 0x0001) {$bits{b1} = 'the read is paired in sequencing, (but not necesarily mapped in a pair)';}
	if ($flag & 0x0002) {$bits{b2} = 'the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment)';}
	if ($flag & 0x0004) {$bits{b4} = 'the query sequence itself is unmapped';}
	if ($flag & 0x0008) {$bits{b8} = 'the pair mate is unmapped';}
	if ($flag & 0x0010) {$bits{b10} = 'mapped to reverse strand';}
	if ($flag & 0x0020) {$bits{b20} = 'pair mate is mapped to reverse strand';}
	if ($flag & 0x0040) {$bits{b40} = 'the read is the first read in a pair';}
	if ($flag & 0x0080) {$bits{b80} = 'the read is the second read in a pair';}
	if ($flag & 0x0100) {$bits{b100} = 'the alignment is not primary (a read having split hits may have multiple primary alignment records)';}
	if ($flag & 0x0200) {$bits{b200} = 'the read fails platform/vendor quality checks';}
	if ($flag & 0x0400) {$bits{b400} = 'the read is either a PCR duplicate or an optical duplicate';}
	
	return(\%bits);	
}



sub numerically
{
	$a <=> $b;
}
