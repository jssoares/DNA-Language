#!/software/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Storable;
use List::MoreUtils qw(uniq);
use List::Util qw(shuffle);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use MIME::Lite;
use Carp;
use lib 'Utils';
use DNAcryptic;
use EnsemblFTP;

local $| = 1;       #setautoflush, outputs immediately instead of buffering

main();

sub main {

  my $mode = undef;
  my $k =undef;
  my $rand_range = undef;
  my $total_number_of_bases = undef;
  my $sequence = undef;
  my $random_sequence = undef;
  my $host = undef;
  my $branch = undef;
  my $rmode = undef;
  my $test = undef;
  my @genetic_code = qw{a t c g};

  my %opts;
  eval {
    GetOptions( \%opts,
		'm=s' => \$mode,               #Options: 'random','test' or 'real'
		'gram=i' => \$k,               #Options: 3 4 5 6 7 8 9 10 11
		'rg=i' => \$rand_range,            #Options: range of the universe of random numbers to choose form when running under random mode. From (4 to infinity really)
		'tnb=i' => \$total_number_of_bases, #Options: to infinity and beyond
		'h=s' => \$host,    #Options: 'ensembl' or 'genomes'
		'b=s' => \$branch,     #Options: 'plants' or 'fungi' or 'protist' or 'metazoa' or 'bacteria'
		'r=s' => \$rmode,      #Options: 'ro', 'rw'
		't=i' => \$test   #Options: 1 (run the test sub_routine) or 0 (don't run it)
	      );
  };
  warn() if $@;

  ($mode,$k,$rand_range,$total_number_of_bases) = checkParameters($mode,$k,$rand_range,$total_number_of_bases);


  if ($mode eq 'test') {

    #my $test_path = '../dna_language_sequences/test/sc_XI.fa';
    #my $test_path = '../dna_language_sequences/test/chrTest/seqExp.fa';
    my $test_path = '../dna_language_sequences/test/chrTest/seqExpSimple.fa';
    open (my $fh, "<", $test_path)
      or croak "cannot open $test_path: $!";
    my @lines = <$fh>;
    close($fh);

    chomp($lines[0]);

    my $chr_label = $lines[0];
    $chr_label =~ s/^>(\d)\sdna:(chr).*/$2$1/;
    $chr_label .= '_test';

    my (@sequence,@sequences);

    for (my $i =0; $i < scalar @lines; $i++) {

      chomp($lines[$i]);

      $lines[$i] = lc($lines[$i]);

      if ($i > 0) {

	@sequence = split(//,$lines[$i]);

	push(@sequences,@sequence);

      }
    }
    undef(@lines);
    print "undefined array lines\n";

    my $base_dist = simpleBaseCounter(\@sequences);


    my %args = (
		   mode => $mode,
		   k => $k,
		   sequence => \@sequences,
		  );

    my $word_count = DNAcryptic->new(%args);

    $word_count->getKgrams;
    $word_count->convertDigitToDNA;
    print Dumper($word_count->digit_to_dna);

  }
  elsif ($mode eq 'real') {

    #Get all chromossomes
    my %args;
    $args{host} = $host; 
    $args{branch} = $branch;
    $args{rmode} = $rmode;
    $args{test} = $test;

    my $ens_ftp = EnsemblFTP->new(%args);
    $ens_ftp->getSpeciesEnsemblFTP;
    $ens_ftp->testFTP;
 
    my $configuration = configuration();
    my ($files_href,$count_files) = getAllFiles($configuration);

  FILE:
    for my $file(sort keys %{$files_href}) {
      my ($dest_href) = cacheDestination($file,$files_href,$configuration,$k);
      my ($lines) = fileLoader($file,$dest_href);
      if (defined $lines->[0]) {
	chomp($lines->[0]);

	my $chr_label = $lines->[0];
	$chr_label =~ s/^>(\d)\sdna:(chr).*/$2$1/;

	my (@sequence,@sequences);

	for (my $i =0; $i < scalar @$lines; $i++) {

	  chomp($lines->[$i]);

	  $lines->[$i] = lc($lines->[$i]);

	  if ($i > 0) {

	    @sequence = split(//,$lines->[$i]);

	    push(@sequences,@sequence);
	  }
	}

	my $base_dist = simpleBaseCounter(\@sequences);
	my %args = (
		    mode => $mode,
		    k => $k,
		    sequence => \@sequences,
		   );

	my $word_count = DNAcryptic->new(%args);

	$word_count->getKgrams;
	$word_count->convertDigitToDNA;
	print "BEGIN $chr_label\n";
	print Dumper($word_count->digit_to_dna);
	print"END $chr_label\n";
      }
      else {
	next FILE;
      }
    }
  }
}

sub configuration {

  my %configuration;
  $configuration{fa_sequences_dir} = '/nfs/cancer_translation/Jorge/dna_language_proj/new_dna_language_sequences/';
  $configuration{analysis_dir} = '/nfs/cancer_translation/Jorge/dna_language_proj/dna_language_analysis/real/';

  return (\%configuration);

}

sub getAllFiles {

  my ($conf) = @_;
  my (%files_hash);
  my (@find_list);

  push(@find_list, `find $conf->{fa_sequences_dir}\* -iname '*.gz'`);
  my $count_files = scalar @find_list;
  for my $file(@find_list) {
    $file =~ s/\n$//;
    my @dirs = split(/\//,$file);
    shift(@dirs);
    for (my $i =0; $i < scalar @dirs; $i++) {
      $files_hash{$file}{$i} = $dirs[$i];
    }
  }

 DEBUG:
  #for my $file(sort keys %files_hash) {
  #print "File: $file\n";
  #for my $dir(sort {$a <=> $b} keys %{$files_hash{$file}}) {
  #print "Dir: $dir\t$files_hash{$file}{$dir}\n";
  #}
  #}

  return (\%files_hash,$count_files);
}

sub cacheDestination {

  my ($file,$files_href,$configuration,$k) = @_;
  my %destination;

  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst ) = localtime ( time ) ;
  my $stamp = sprintf "%02d-%02d-%02d_%02d-%02d-%4d_", $hour,$min,$sec,$mday,$mon+1,1900+$year;


  my $new_dir1 = $configuration->{analysis_dir} . 'k' . $k . '/';
  $new_dir1 .= $files_href->{$file}->{5} . '/';

  my $new_dir2 = $new_dir1;
  $new_dir2 .= $files_href->{$file}->{6} . '/';

  my $analysis_label = $files_href->{$file}->{5} . '_' . $files_href->{$file}->{6} . '_analysis_' . $stamp . '.txt';

  $destination{$file}{NEW_DIR1} = $new_dir1;
  $destination{$file}{NEW_DIR2} = $new_dir2;
  $destination{$file}{OUTPUT_FILE} = $new_dir2 . $analysis_label;

  #print Dumper(\%destination);
  #exit;
  return (\%destination);
}


sub fileLoader {

  my ($file) = @_;
  my $output = $file;
  $output =~ s/\.gz$//;
  my $gun_status = gunzip $file => $output
    or warn "gunzip failed: $GunzipError\n";

  if (-e $output) {
    open (my $fh, "<", $output)  or croak "cannot open < $output: $!";
    print "loading $output\n";
    my @lines = <$fh>;
    print "output is loaded\n";
    close($fh);

    fileCleanUp($output);

    return (\@lines);
  }
  else {
    return(qw());
  }
}

sub fileCleanUp {

  my $file_to_rm = shift;
 SYSTEM_CALL:
  `rm $file_to_rm`;

  return;
}

sub simpleBaseCounter {

  my ($sequences) = @_;

  my $total_number_of_bases = scalar @$sequences;

  my ($countA,$countT,$countC,$countG,$countN)  = (0,0,0,0,0);
  for my $base (@$sequences) {
    if ($base eq 'a') {
      $countA++;
    } elsif ($base eq 't') {
      $countT++;
    } elsif ($base eq 'c') {
      $countC++;
    } elsif ($base eq 'g') {
      $countG++;
    } elsif ($base eq 'n') {
      $countN++;
    }
  }

  my %base_dist;
  $base_dist{'a'} = $countA/$total_number_of_bases;
  $base_dist{'t'} = $countT/$total_number_of_bases;
  $base_dist{'c'} = $countC/$total_number_of_bases;
  $base_dist{'g'} = $countG/$total_number_of_bases;
  $base_dist{'n'} = $countN/$total_number_of_bases;

  print (scalar @$sequences,"\n");
  print Dumper(\%base_dist);

  return(\%base_dist);
}

sub checkParameters {

  my ($mode,$k,$rand_range,$total_number_of_bases) = @_;
  my $answer = '0';
  if (!$mode) {
    print "You forgot too specify the mode you want to use.\nWould you like to specify it now? (y/n)\n";
    $answer = allowManSet();

    if ($answer =~ m/y/i) {

      $mode = setManMode();
      	print "Mode has been set to: $mode\n\n";

    } elsif ($answer =~ m/n/i){
      die "You didn't want to define the mode.\nBye";
    }
    else {
      ($mode,$k,$rand_range,$total_number_of_bases) = giveMeAnotherChance($answer,$mode,$k,$rand_range,$total_number_of_bases);
    }
  }

  if ($mode eq 'real' && $k) {
    print "I have all options I need. Will continue analysis\n\n";
  } elsif ($mode eq 'real' && !$k) {

    print "You forgot too specify the kgram you want to use.\nWould you like to specify it now? (y/n)\n";
    $answer = allowManSet();

    if ($answer =~ m/y/i) {

      $k = setManKgram();
      print "I have all options I need. Will run on $mode mode and kgram: $k\n\n";
    } elsif ($answer =~ m/n/i) {
      die "You didn't want to define the kgram analysis.\nBye";
    } else {
      ($mode,$k,$rand_range,$total_number_of_bases) = giveMeAnotherChance($answer,$mode,$k,$rand_range,$total_number_of_bases);
    }
  } elsif ($mode eq 'test' && $k) {
    print "I have all options I need. Will run on $mode mode and kgram: $k\n";
  } elsif ($mode eq 'test' && !$k) {

    print "You forgot too specify the kgram you want to use.\nWould you like to specify it now? (y/n)\n";
    $answer = allowManSet();

    if ($answer =~ m/y/i) {

      $k = setManKgram();
      print "I have all options I need. Will run on $mode mode and kgram: $k\n\n";
    } elsif ($answer =~ m/n/i) {
      die "You didn't want to define the kgram analysis.\nBye";
    } else {
      ($mode,$k,$rand_range,$total_number_of_bases) = giveMeAnotherChance($answer);
    }
  } elsif ($mode eq 'random') {

    if ($k) {
      print "Kgram: $k\n";
    } else {
      print "You forgot too specify the kgram you want to use.\nWould you like to specify it now? (y/n)\n";
      $answer = allowManSet();

      if ($answer =~ m/y/i) {

	$k = setManKgram();
	print "Kgram has been set to $k\n\n";
      } elsif ($answer =~ m/n/i) {
	die "You didn't want to define the kgram analysis.\nBye";
      } else {
	($mode,$k,$rand_range,$total_number_of_bases) = giveMeAnotherChance($answer,$mode,$k,$rand_range,$total_number_of_bases);
      }
    }

    if ($rand_range) {
      print "Rand_range: $rand_range\n";
    } else {
      print "You forgot too specify the rand_range you want to use.\nWould you like to specify it now? (y/n)\n";
      $answer = allowManSet();

      if ($answer =~ m/y/i) {

	$rand_range = setManRandRange();
	print "Rand_Range has been set to $rand_range\n\n";

      } elsif ($answer =~ m/n/i) {
	die "You didn't want to define the random range.\nBye";
      } else {
	($mode,$k,$rand_range,$total_number_of_bases) = giveMeAnotherChance($answer,$mode,$k,$rand_range,$total_number_of_bases);
      }
    }

    if ($total_number_of_bases) {
      print "Total number of bases: $total_number_of_bases\n";
    } else {
      print "You forgot too specify the total number of bases you want to use.\nWould you like to specify it now? (y/n)\n";
      $answer = allowManSet();

      if ($answer =~ m/y/i) {

	$total_number_of_bases = setManTNB();
	print "Total number of bases has been set to $total_number_of_bases\n\n";
      } elsif ($answer =~ m/n/i) {
	die "You didn't want to define the total number of bases.\nBye";
      } else {
	($mode,$k,$rand_range,$total_number_of_bases) = giveMeAnotherChance($answer,$mode,$k,$rand_range,$total_number_of_bases);
      }
    }
  }
  return($mode,$k,$rand_range,$total_number_of_bases);
}

sub allowManSet {


  my $answer = <STDIN>;
  chomp($answer);
  return ($answer);
}

sub setManKgram {
  print "Choose kgram for this analysis? (3 4 5 6 7 8 9 10 11)\n";
  my $k = <STDIN>;
  chomp($k);
  return ($k);
}

sub setManRandRange {

  print "Integer to set as the iteration universe of random numbers.(100 to Infinity)\n";
  my $rand_range = <STDIN>;
chomp($rand_range);
  return ($rand_range);

}

sub setManTNB {

  print "Integer to set as the length of the sequence you want to create.(3 to Infinity)\n";
  my $tnb = <STDIN>;
  chomp($tnb);
  return ($tnb);

}

sub setManMode {

  print "Choose a mode to run the analysis under.(real random test)\n";
  my $mode = <STDIN>;
  chomp($mode);
  return ($mode);

}

sub giveMeAnotherChance {
  my ($answer,$mode,$k,$rand_range,$total_number_of_bases) = @_;
  print "You answered: $answer\nNo option is set for this.\nDo you want another chance?(y/n)\n";
  my $answer2 = <STDIN>;
  chomp($answer2);
  if ($answer2 =~ m/y/i) {
    ($mode,$k,$rand_range,$total_number_of_bases) = checkParameters($mode,$k,$rand_range,$total_number_of_bases);
    return($mode,$k,$rand_range,$total_number_of_bases);
  }
  elsif ($answer2 =~ m/n/i) {
    die "You answered - $answer2 -\nGoodbye";
  }
  else {
    die "Stop wasting my time, bye!";
  }
}


sub displayOptions {

  print "\nYou can pass several flags to this script. Here's the list:\n\n";
  print ("\t",'-m => ', 'sets $mode. Three modes are available (random,test or real)',"\n\n");
  print ("\t",'-gram => ', 'sets $k. An integer that specifies the kgram you want to run. Supported (3 4 5 6 7 8 9 10 11)',"\n\n");
  print ("\t",'-rg => ', 'sets $rand_range. An integer that specifies the range of the universe of random numbers to choose from when running under random mode',"\n\n");
  print ("\t",'-tnb => ', 'sets $total_number_of_bases. An integer that specifies the length of the random sequence you want to analyse when running under random mode',"\n\n");
  print ("\t",'-h => ', 'displays this list',"\n\n");
  exit;
}

sub smallTester {

  my @sequence = qw(a a a c t g t a a a c t g);
  return(\@sequence);
}

sub sequenceTest {

  my %sequences = undef;
  $sequences{'seq1'} = '';
  return(\%sequences);

}

sub test {

my $x = 0;
my $string = "Hello";
print "$x foda-se $string\n";

return;
}

