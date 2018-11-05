#!/usr/local/bin/perl -w

use strict;
use Data::Dumper;

use Getopt::Long;
use Config::Simple;

&GetOptions (
		"c=s"=> \my $config
	     );
unless($config){die "ERROR: no config file entered\nusage example: SCgenotype.pl -c /full/path/to/SCpipe.conf\n\n";}
unless(-e $config)
{
	warn "ERROR: could not find config: $config\n";
	warn "Please give full path to config file, SCpipe.conf.\n\n";
	exit;
}


my $cfg = new Config::Simple($config);
my %config = $cfg->vars();


#--------------------------------------------------------------

# check for output dir:
my $analysis_dir = $config{ANALYSIS_DIR}; 
unless(-e $analysis_dir)
{
	warn "ERROR: could not find ANALYSIS_DIR: $analysis_dir\n";
	warn "Please create this dir and try again.\n\n";
	exit;
}
# log file:
my $log_file = $analysis_dir . 'SCpipe.log';
open (LOG, ">$log_file") or die "ERROR: cannot open log file: $log_file: $!\n\n";
my $timestamp = &time_and_date();
print LOG "SCpipe\n$timestamp\n----------\n";

# make analysis dirs:
warn "making output dirs in: $analysis_dir\n";
&make_analysis_dirs();

# reading variants file:
warn "reading variants file: $config{VARIANTS}\n";
my %vars = %{&read_variants_file()};

# read fastq_dir - find legitimate samples:
my %fastqs;
my %samples;
my $fastq_dir = $config{FASTQ_DIR}; 
warn "reading fastq data: $fastq_dir\n";
print LOG "reading fastq data: $fastq_dir\n";
opendir(DIR, $fastq_dir);
my @files = readdir(DIR);
closedir(DIR); 
my $fastq_fwd_suffix = $config{FASTQ_FWD_SUFFIX};
my $fastq_rev_suffix = $config{FASTQ_REV_SUFFIX};
foreach my $file (@files)
{
	if ($file eq "." or $file eq "..") {next;}	
	if ($file =~ m/^(.*)$fastq_fwd_suffix$/) {$fastqs{$1}{f} = '';}
	elsif ($file =~ m/^(.*)$fastq_rev_suffix$/) {$fastqs{$1}{r} = '';}
}
my $sample_pair_count = 0;
foreach my $sample_name (keys %fastqs)
{
	if ( (exists($fastqs{$sample_name}{f})) and (exists($fastqs{$sample_name}{r})) )
	{
		$samples{$sample_name} = '';
	}
}
warn 'found ' . scalar(keys %samples) . ' sample pairs' . "\n";
print LOG 'found ' . scalar(keys %samples) . ' sample pairs' . "\n";

my $amplcn_count = 0;
foreach my $amplcn (sort keys %samples)
{
	
	$amplcn_count ++;
	warn "\n---\n($amplcn_count) sample: $amplcn\n";
	print LOG "\n---\n($amplcn_count) sample: $amplcn:\n";
	
	# STAR alignment
	my $amplcn_out_dir = $analysis_dir . 'STAR/output/' . $amplcn . '/';
	my $amplcn_tmp_dir = $analysis_dir . 'STAR/tmp/' . $amplcn . '/';
	unless (-e $amplcn_out_dir)
	{
		my $dir_cmd = "mkdir $amplcn_out_dir";
		system($dir_cmd) == 0 or die "system cmd [$dir_cmd] failed ($?): $!";		
	}
	
	# check to see if this job has already been carried out
	my $job_done_file = $amplcn_out_dir . 'Log.final.out';
	if ((-e $job_done_file) and (-s $job_done_file))
	{
		warn 'using existing STAR alignment data' . "\n";
		print LOG 'using existing STAR alignment data: ' . $amplcn_out_dir . "\n";
	}
	else
	{ 
		my $fwd_fastq = $fastq_dir . $amplcn . $fastq_fwd_suffix;
		my $rev_fastq = $fastq_dir . $amplcn . $fastq_rev_suffix;
		my $STAR_cmd = $config{STAR_EXECUTABLE};
		$STAR_cmd .= ' ' . '--limitBAMsortRAM 10000000000';
		if ($amplcn_count < scalar(keys %samples)) { $STAR_cmd .= ' ' . '--genomeLoad LoadAndKeep'; }
		else { $STAR_cmd .= ' ' . '--genomeLoad LoadAndRemove'; }
		$STAR_cmd .= ' ' . '--genomeDir ' . $config{STAR_INDEX_DIR};
		$STAR_cmd .= ' ' . "--readFilesIn $fwd_fastq $rev_fastq";
		$STAR_cmd .= ' ' . '--readFilesCommand zcat';
		$STAR_cmd .= ' ' . "--outFileNamePrefix $amplcn_out_dir";
		$STAR_cmd .= ' ' . '--runThreadN ' . $config{STAR_NUM_THREADS};
		$STAR_cmd .= ' ' . "--outTmpDir $amplcn_tmp_dir";
		$STAR_cmd .= ' ' . '--outReadsUnmapped Fastx';
		$STAR_cmd .= ' ' . '--outSAMtype BAM Unsorted';
	
		warn "STAR alignment\n";
		print LOG 'running STAR using: ' . $STAR_cmd . "\n";
		system($STAR_cmd) == 0 or die "system cmd [$STAR_cmd] failed ($?): $!";		
	}
	
	
	# Sort and index bam files
	my $sort_out_dir = $analysis_dir . 'SortIndex/output/' . $amplcn . '/';
	unless (-e $sort_out_dir)
	{
		my $dir_cmd = "mkdir $sort_out_dir";
		system($dir_cmd) == 0 or die "system cmd [$dir_cmd] failed ($?): $!";		
	}

	my $input_filename = $amplcn_out_dir . 'Aligned.out.bam';
	my $output_sort_filename = $sort_out_dir . $amplcn . '.sorted.bam';
	if ((-e $output_sort_filename) and (-s $output_sort_filename))
	{
		print LOG 'using sorted bam file: ' . $output_sort_filename . "\n";
	}
	else
	{
		# sort:
		my $sort_cmd = $config{SAMTOOLS_EXECUTABLE} . ' sort'; 
		$sort_cmd .= ' ' . '-@ ' . $config{SAMTOOLS_NUM_THREADS};
		$sort_cmd .= ' ' . '-m 4G';
		$sort_cmd .= ' ' . $input_filename;
		$sort_cmd .= ' ' . '-f ' . $output_sort_filename;
	
		warn "SAMtools sort\n";
		print LOG 'sorting bam file using: ' . $sort_cmd . "\n";
		system($sort_cmd) == 0 or die "system cmd [$sort_cmd] failed ($?): $!";
		
		# then, index:
		my $index_cmd = $config{SAMTOOLS_EXECUTABLE} . " index $output_sort_filename";		
		warn "SAMtools index\n";
		print LOG 'indexing bam file using: ' . $index_cmd . "\n";
		system($index_cmd) == 0 or die "system cmd [$index_cmd] failed ($?): $!";
	}

	
	# separate and process genomic and mRNA
	my %sep_output_dirs;
	$sep_output_dirs{gDNA} = $analysis_dir . 'SEP_Amplicons/gDNA/' . $amplcn . '/';
	$sep_output_dirs{mRNA} = $analysis_dir . 'SEP_Amplicons/mRNA/' . $amplcn . '/';
	
	foreach my $sep_data_type (sort keys %sep_output_dirs)
	{
		my $sep_out_dir = $sep_output_dirs{$sep_data_type};
		unless (-e $sep_out_dir)
		{
			my $sep_dir_cmd = "mkdir $sep_out_dir";
			system($sep_dir_cmd) == 0 or die "system cmd [$sep_dir_cmd] failed ($?): $!";		
		}
		
		# separate: 
		my ($primers_bed_file, $sep_script);
		if ($sep_data_type eq 'gDNA')
		{
			$primers_bed_file = $config{gPRIMERS_BED};
			$sep_script = $config{SEPARATE_gDNA_SCRIPT};
		}
		else
		{
			$primers_bed_file = $config{mPRIMERS_BED};
			$sep_script = $config{SEPARATE_mRNA_SCRIPT};
		}
		
		my $sep_cmd = 'perl ' . $sep_script;
		$sep_cmd .= ' ' . '-c ' . $primers_bed_file;
		$sep_cmd .= ' ' . '-b ' . $sort_out_dir;
		$sep_cmd .= ' ' . '-o ' . $sep_out_dir;
		warn "separate $sep_data_type\n";
		print LOG "separate $sep_data_type using: " . $sep_cmd . "\n";
		system($sep_cmd) == 0 or die "system cmd [$sep_cmd] failed ($?): $!";
		
		# copy the header from original bam file:
		my $header_file = $sep_out_dir . 'header.txt';
		my $head_cmd = $config{SAMTOOLS_EXECUTABLE} . ' view -H ' . $output_sort_filename . ' > ' . $header_file;
		warn "$sep_data_type, generating header file\n";
		print LOG "$sep_data_type, generating header file using: " . $head_cmd . "\n";
		system($head_cmd) == 0 or die "system cmd [$head_cmd] failed ($?): $!";
		
		# concatenating the header with the .sam file generated during separation
		my $sep_sam_file = $sep_out_dir . $amplcn . '.sorted.genome.sam';
		my $sep_header_sam_file = $sep_out_dir . $amplcn . '.header.sam';
		my $cat_cmd = 'cat ' . $header_file . ' ' . $sep_sam_file . ' > ' . $sep_header_sam_file;
		warn "$sep_data_type, cat header and sam file\n";
		print LOG "$sep_data_type, cat header and sam file using: " . $cat_cmd . "\n";
		system($cat_cmd) == 0 or die "system cmd [$cat_cmd] failed ($?): $!";
	
		# converting the sam file to bam file
		my $sep_bam_file = $sep_out_dir . $amplcn . '.bam';
		my $sam_to_bam_cmd = $config{SAMTOOLS_EXECUTABLE} . ' view -Sb ' . $sep_header_sam_file . ' > ' . $sep_bam_file;
		warn "$sep_data_type, convert sam to bam\n";
		print LOG "$sep_data_type, convert sam to bam using: " . $sam_to_bam_cmd . "\n";
		system($sam_to_bam_cmd) == 0 or die "system cmd [$sam_to_bam_cmd] failed ($?): $!";
		
		# sorting bam file
		my $sep_sorted_bam_file = $sep_out_dir . $amplcn . '.sorted.bam';
		my $sep_sort_cmd = $config{SAMTOOLS_EXECUTABLE} . ' sort ' . $sep_bam_file . ' -f ' . $sep_sorted_bam_file;
		warn "$sep_data_type, sorting bam file\n";
		print LOG "$sep_data_type, sorting bam file using: " . $sep_sort_cmd . "\n";
		system($sep_sort_cmd) == 0 or die "system cmd [$sep_sort_cmd] failed ($?): $!";
	
		# indexing bam file
		my $sep_index_cmd = $config{SAMTOOLS_EXECUTABLE} . ' index ' . $sep_sorted_bam_file;
		warn "$sep_data_type, indexing bam file\n";
		print LOG "$sep_data_type, indexing bam file using: " . $sep_index_cmd . "\n";
		system($sep_index_cmd) == 0 or die "system cmd [$sep_index_cmd] failed ($?): $!";
		
		# tidying up
		my @intermediate_files = ($sep_sam_file, $header_file, $sep_header_sam_file, $sep_bam_file);
		foreach my $intermediate_file (@intermediate_files)
		{
			my $tidy_cmd = 'rm ' . $intermediate_file;
			print LOG "$sep_data_type, tidying intermediate files using: " . $tidy_cmd . "\n";
			system($tidy_cmd) == 0 or die "system cmd [$tidy_cmd] failed ($?): $!";	
		}
	}
		
		
	# mpileup:
	my %mpileup_output_dirs;
	$mpileup_output_dirs{gDNA} = $analysis_dir . 'mpileup/gDNA/' . $amplcn . '/';
	$mpileup_output_dirs{mRNA} = $analysis_dir . 'mpileup/mRNA/' . $amplcn . '/';
	
	foreach my $mpileup_data_type (sort keys %mpileup_output_dirs)
	{
		# make the output dir:
		my $mpileup_out_dir = $mpileup_output_dirs{$mpileup_data_type};
		unless (-e $mpileup_out_dir)
		{
			my $mpileup_dir_cmd = "mkdir $mpileup_out_dir";
			system($mpileup_dir_cmd) == 0 or die "system cmd [$mpileup_dir_cmd] failed ($?): $!";		
		}

		# run mpileup for each variant in the variants file:
		foreach my $var_chr (sort keys %vars)
		{
			foreach my $var_genome_start (sort numerically keys %{$vars{$var_chr}})
			{
				foreach my $var_genome_end (sort numerically keys %{$vars{$var_chr}{$var_genome_start}})
				{
					my $var_name = $vars{$var_chr}{$var_genome_start}{$var_genome_end}{name};
					my $var_type = $vars{$var_chr}{$var_genome_start}{$var_genome_end}{type};
					
					my $mpileup_posn = $var_chr . ':' . $var_genome_start . '-' . $var_genome_end;
					my $mpileup_out_file = $mpileup_out_dir . $amplcn . '.' . $var_name . '.mpileup';
					
					my $sep_sorted_bam_file = $analysis_dir . 'SEP_Amplicons/' . $mpileup_data_type . '/' . $amplcn . '/' . $amplcn . '.sorted.bam';
					
					my $mpileup_cmd = $config{SAMTOOLS_EXECUTABLE} . ' mpileup';
					$mpileup_cmd .= ' ' . '--count-orphans';
					$mpileup_cmd .= ' ' . '--max-depth 999999';
					$mpileup_cmd .= ' ' . '--min-BQ 30';
					$mpileup_cmd .= ' ' . '--ignore-overlaps';
					$mpileup_cmd .= ' ' . '--region ' . $mpileup_posn;
					$mpileup_cmd .= ' ' . '--output ' . $mpileup_out_file;
					$mpileup_cmd .= ' ' . $sep_sorted_bam_file;

					warn "$mpileup_data_type, mpileup: $var_name ($mpileup_posn)\n";
					print LOG "$mpileup_data_type, mpileup: $var_name ($mpileup_posn) using:\n\t" . $mpileup_cmd . "\n";
					system($mpileup_cmd) == 0 or die "system cmd [$mpileup_cmd] failed ($?): $!";
				}
			}
		}
	}
	
}

# summarise:
my @separation_types = ('gDNA', 'mRNA');
foreach my $var_chr (sort keys %vars)
{
	foreach my $var_genome_start (sort numerically keys %{$vars{$var_chr}})
	{
		foreach my $var_genome_end (sort numerically keys %{$vars{$var_chr}{$var_genome_start}})
		{
			my $var_name = $vars{$var_chr}{$var_genome_start}{$var_genome_end}{name};
			my $var_type = $vars{$var_chr}{$var_genome_start}{$var_genome_end}{type};
			
			foreach my $separation_type (sort @separation_types)
			{
				my $summarise_out_file = $analysis_dir . 'Summarize/' . $separation_type . '_' . $var_name . '.txt';
				open (SUMM, ">$summarise_out_file") or die "cannot open summarise out file: $summarise_out_file: $!\n";
				print SUMM "id\tbp\tref\tA\tG\tC\tT\tdel\tins\tinserted\tambiguous\n";
				close SUMM;
				
				warn "summarising: $var_name, $separation_type\n";					
				
				foreach my $amplcn (sort keys %samples)
				{
					
					my $mpileup_out_file = $analysis_dir . 'mpileup/' . $separation_type . '/' . $amplcn . '/' . $amplcn . '.' . $var_name . '.mpileup';
					my $summarise_cmd = 'python ' . $config{SUMMARISE_SCRIPT} . ' ' . $mpileup_out_file . ' >> ' . $summarise_out_file;
					
					print LOG "summarise: $var_name, $amplcn, $separation_type, using: " . $summarise_cmd . "\n";
					system($summarise_cmd) == 0 or die "system cmd [$summarise_cmd] failed ($?): $!";
				}
			}
		}
	}
}

exit;






#------------------------------------------------------------------------------



sub make_analysis_dirs
{
	my @reqrd_dirs = ('STAR', 'STAR/output', 'STAR/tmp', 'SortIndex', 'SortIndex/output',
					  'SEP_Amplicons', 'SEP_Amplicons/gDNA', 'SEP_Amplicons/mRNA',
					  'mpileup', 'mpileup/gDNA', 'mpileup/mRNA', 'Summarize');
	unless($analysis_dir =~ m/\/$/){$analysis_dir .= '/';}
	
	foreach my $dir (@reqrd_dirs)
	{
		my $new_dir = $analysis_dir . $dir;
		if (-e $new_dir)
		{
			warn "WARNING: analysis output dir exists ($new_dir) - potentially overwriting\n";
			print LOG "WARNING: analysis output dir exists ($new_dir) - potentially overwriting\n";
		}
		else
		{
			my $dir_cmd = "mkdir $new_dir";
			system($dir_cmd) == 0 or die "system cmd [$dir_cmd] failed ($?): $!";
		}
	}
}

sub read_variants_file
{
	# read variants file
	my %vars;
	open (VARS, "<$config{VARIANTS}") or die "cannot open $config{VARIANTS}: $!\n";
	while (<VARS>)
	{
		my $line = $_;
		chomp $line;
		if ($line =~ m/^#/) { next; }    # ignore the header line
		unless ($line){next;}
		
		my ($genome_acc, $genome_start, $genome_end, $name, $type) = split(/\t/, $line);
		$vars{$genome_acc}{$genome_start}{$genome_end}{name} = $name;
		$vars{$genome_acc}{$genome_start}{$genome_end}{type} = $type;
	}
	close VARS;
	
	return(\%vars);
}


sub time_and_date
{
	my ($second, $minute, $hour, $dayofmonth, $month, $year, $weekday, $dayofyeaf, $IsDST) = localtime(time);
	my %months;
	$months{'01'} = 'Jan';
	$months{'02'} = 'Feb';
	$months{'03'} = 'Mar';
	$months{'04'} = 'Apr';
	$months{'05'} = 'May';
	$months{'06'} = 'Jun';
	$months{'07'} = 'Jul';
	$months{'08'} = 'Aug';
	$months{'09'} = 'Sep';
	$months{'10'} = 'Oct';
	$months{'11'} = 'Nov';
	$months{'12'} = 'Dec';
	$dayofmonth = sprintf( "%02d", $dayofmonth);    
	$month = sprintf( "%02d", ($month + 1));  
	$second = sprintf( "%02d", ($second + 1));  
	$minute = sprintf( "%02d", ($minute + 1));  
	$hour = sprintf( "%02d", ($hour + 1));  
	my $curr_month = $months{$month}; 
	$year = $year + 1900;
	
	my $timestamp = "$dayofmonth $curr_month $year - $hour:$minute:$second";
	return($timestamp);
}


sub numerically
{
	$a <=> $b;
}
