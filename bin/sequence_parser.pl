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
use Switch;
use Carp;

local $| = 1;       #setautoflush, outputs immediately instead of buffering

main();

sub main {

  my $mode = undef;
  my $k =undef;
  my $rand_range = undef;
  my $total_number_of_bases = undef;
  my $sequence = undef;
  my $random_sequence = undef;
  my @genetic_code = qw{a t c g};

  my %opts;
  eval {
    GetOptions( \%opts,
		'm=s' => \$mode,               #Options: 'random','test' or 'real'
		'gram=i' => \$k,               #Options: 3 4 5 6 7 8 9 10 11
		'rg=i' => \$rand_range,            #Options: range of the universe of random numbers to choose form when running under random mode. From (4 to infinity really)
		'tnb=i' => \$total_number_of_bases #Options: to infinity and beyond
	      );
  };
  warn() if $@;

  ($mode,$k,$rand_range,$total_number_of_bases) = checkParameters($mode,$k,$rand_range,$total_number_of_bases);

  if ($mode eq 'real') {
    my $status = 0;
    my $configuration = configuration();
    my ($files_href,$count_files) = getAllFiles($configuration);
    #print Dumper($files_href);
    #exit;
  FILE:
    for my $file(sort keys %{$files_href}) {
      my ($dest_href) = cacheDestination($file,$files_href,$configuration,$k);
      my ($lines) = fileLoader($file,$dest_href);
      my ($kgram_conv) = analyseFile($lines,$file,$dest_href,$k,$mode,$total_number_of_bases);
      #print Dumper($kgram_conv);
      $status += createAnalysisPackage($kgram_conv,$file,$dest_href,$k);
      next FILE;
    }
    if ($status == $count_files) {
      sendMail($mode,$k,$count_files);
    }
  }
  elsif ($mode eq 'test') {

    open (my $fh, "<", "../dna_language_sequences/test/chrTest/seqExp.fa")
      or croak "cannot open < ../dna_language_sequences/test/chrTest/seqExp.fa: $!";
    print "loading /dna_language_sequences/test/chrTest/seqExp.fa\n";
    my @lines = <$fh>;
    print "seqExp.fa is loaded\n";
    close($fh);
    #print "$lines[0]\n";
    chomp($lines[0]);
    #print "$lines[0]\n";
    my $chr_label = $lines[0];
    $chr_label =~ s/^>(\d)\sdna:(chr).*/$2$1/;
    $chr_label .= '_test';
    #print "$chr_label\n";
    #exit;
    my @sequence;
    my @sequences;
    for (my $i =0; $i < scalar @lines; $i++) {
      chomp($lines[$i]);
      $lines[$i] = lc($lines[$i]);
      if ($i > 0) {
	#print "$i\n";
	@sequence = split(//,$lines[$i]);

	#This implementation needs to be uncommented if running the code on base4. If running on base5 we don't need this
	#as it deals with N's and n's 
	#my @sequence2;
	#for (my $j = 0; $j < scalar @sequence; $j++) {
	#  if($sequence[$j] ne 'N' && $sequence[$j] ne 'n') {
	#    push(@sequence2,$sequence[$j]);
	#  }
	#}
	push(@sequences,@sequence);
	#print "$line";
      }
    }
    undef(@lines);
    print "undefined array lines\n";
    #exit;
    $total_number_of_bases = scalar @sequences;
    my ($countA,$countT,$countC,$countG,$countN)  = (0,0,0,0,0);
    for my $base(@sequences) {
      if ($base eq 'a') {
	$countA++;
      }
      elsif ($base eq 't') {
	$countT++;
      }
      elsif ($base eq 'c') {
	$countC++;
      }
      elsif ($base eq 'g') {
	$countG++;
      }
      elsif ($base eq 'n') {
	$countN++;
      }
    }

    my %base_dist;
    $base_dist{'a'} = $countA/$total_number_of_bases;
    $base_dist{'t'} = $countT/$total_number_of_bases;
    $base_dist{'c'} = $countC/$total_number_of_bases;
    $base_dist{'g'} = $countG/$total_number_of_bases;
    $base_dist{'n'} = $countN/$total_number_of_bases;

    print (scalar @sequences,"\n");
    print Dumper(\%base_dist);
    my ($kgram_analysis) = getKgrams(\@sequences,$k,$mode);
    #print Dumper($kgram_analysis);
    my ($kgram_conv) = convertDigitToDNA($kgram_analysis,$k,$mode);
    print Dumper($kgram_conv);
  }
  elsif ($mode eq 'random') {

    my ($configuration) = configurationRandom();
    ($sequence) = randomBaseAssignment($rand_range,\@genetic_code, $total_number_of_bases,$mode);
    #print "@$sequence\n";
    my ($countA,$countT,$countC,$countG)  = (0,0,0,0);
    for my $base(@$sequence) {
      if ($base eq 'a') {
	$countA++;
      }
      elsif ($base eq 't') {
	$countT++;
      }
      elsif ($base eq 'c') {
	$countC++;
      }
      elsif ($base eq 'g') {
	$countG++;
      }
    }

    my %base_dist;
    $base_dist{'a'} = $countA/$total_number_of_bases;
    $base_dist{'t'} = $countT/$total_number_of_bases;
    $base_dist{'c'} = $countC/$total_number_of_bases;
    $base_dist{'g'} = $countG/$total_number_of_bases;

    print (scalar @$sequence,"\n");
    print Dumper(\%base_dist);

    my $chr_label = 'random100million';
    my ($kgram_analysis) = getKgrams($sequence,$k,$mode);
    #print Dumper($kgram_analysis);

    my ($kgram_conv) = convertDigitToDNA($kgram_analysis,$k,$mode);
    print Dumper($kgram_conv);
  }
  # eval {
  #   test();
  # };
  # warn() if $@;

}

sub createAnalysisPackage {

  my ($kgram_conv,$file,$dest_href,$k) = @_;
  my $status = 0;

  mkdir($dest_href->{$file}->{NEW_DIR1});
  mkdir($dest_href->{$file}->{NEW_DIR2});

  my $fh;
  open ($fh, ">", $dest_href->{$file}->{OUTPUT_FILE})  or die "Can't open .$!";
  #print "FILE: $dest_href->{$file}->{OUTPUT}\n";
  #print "Header: $dest_href->{$file}->{HEADER}\n";
  #print "k$k\tcount\n";

  #for my $key(sort keys %{$kgram_conv}) {
  #  print "$key\t$kgram_conv->{$key}\n";
  #}

  print $fh ("$dest_href->{$file}->{HEADER}\n");
  print $fh ("k$k\tcount\n");
  for my $key(sort keys %{$kgram_conv}) {
    print $fh "$key\t$kgram_conv->{$key}\n";
  }

  close($fh);

  $status = 1;

  return($status);
}

sub configuration {

  my %configuration;
  $configuration{fa_sequences_dir} = '/nfs/cancer_translation/Jorge/dna_language_proj/dna_language_sequences/';
  $configuration{analysis_dir} = '/nfs/cancer_translation/Jorge/dna_language_proj/dna_language_analysis/real/';

  return (\%configuration);

}

sub configurationRandom {

  my %configuration;
  $configuration{analysis_dir} = '/nfs/cancer_translation/Jorge/dna_language_proj/dna_language_analysis/random/';

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
    or die "gunzip failed: $GunzipError\n";

  open (my $fh, "<", $output)  or croak "cannot open < $output: $!";
  print "loading $output\n";
  my @lines = <$fh>;
  print "output is loaded\n";
  close($fh);

  fileCleanUp($output);

  return (\@lines);
}

sub fileCleanUp {

  my $file_to_rm = shift;
 SYSTEM_CALL:
  `rm $file_to_rm`;

  return;
}

sub analyseFile{

  my ($lines,$file,$dest_href,$k,$mode,$total_number_of_bases) = @_;

  $lines->[0] =~ s/\n$//;
  $dest_href->{$file}{HEADER} = $lines->[0];

  my @sequence;
  my @sequences;
  for (my $i =1; $i < scalar @$lines; $i++) {
    $lines->[$i] =~ s/\n$//;
    $lines->[$i] = lc($lines->[$i]);
    @sequence = split(//,$lines->[$i]);
    push(@sequences,@sequence);
  }
  $total_number_of_bases = scalar @sequences;
  my ($countA,$countT,$countC,$countG,$countN)  = (0,0,0,0,0);
  for my $base (@sequences) {
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

  print (scalar @sequences,"\n");
  print Dumper(\%base_dist);
  my ($kgram_analysis) = getKgrams(\@sequences,$k,$mode);
  #print Dumper($kgram_analysis);
  my ($kgram_conv) = convertDigitToDNA($kgram_analysis,$k,$mode);
  #print Dumper($kgram_conv);
  return($kgram_conv);
}


sub getKgrams {

  my ($sequence,$k,$mode) = @_;
  #print Dumper($sequence);
  #print "$mode\n";
  my @kgram;
  my %base_hash;
  if ($mode eq 'test' || $mode eq 'real') {
    $base_hash{a} = 0;
    $base_hash{c} = 1;
    $base_hash{g} = 2;
    $base_hash{t} = 3;
    $base_hash{n} = 4;
  }
  elsif ($mode eq 'random') {
    $base_hash{a} = 0;
    $base_hash{c} = 1;
    $base_hash{g} = 2;
    $base_hash{t} = 3;
  }

  my $kminus1 = $k - 1;
  print "KMINUS1: $kminus1\n";
  my $filter = scalar @$sequence - $kminus1;
  my @exps = (0..$kminus1);
  @exps = reverse(@exps);
  for (my $b = 0; $b < $filter; $b++) {
    my @base4;
    my @base5;
    my $word;
    for(my $i = 0; $i < $k; $i++) {
      if ($mode eq 'test' || $mode eq 'real') {
	push(@base5,$base_hash{$sequence->[$b + $i]});
      }
      elsif ($mode eq 'random') {
	push(@base4,$base_hash{$sequence->[$b + $i]});
      }
    }

    my $i_dec;
    if ($mode eq 'test' || $mode eq 'real') {
      $i_dec = base5ToDec(\@base5,\@exps);
      $kgram[$i_dec]++;
    }
    elsif ($mode eq 'random') {
      $i_dec = base4ToDec(\@base4,\@exps);
      $kgram[$i_dec]++;
    }
  }
  return(\@kgram);
}

sub convertDigitToDNA {

  my($kgram_analysis,$k,$mode) = @_;
  #print Dumper($kgram_analysis);
  #exit;
  my %conv_hash;
  my @bases;
  my $base;
  #if ($mode eq 'test' || $mode eq 'real') {
    @bases = qw{a c g t n};
    $base = 5;
  #}
  #elsif ($mode eq 'random') {
  #  @bases = qw{a c g t};
  #  $base = 4;
  #}

  my @k_factors;
  for (my $i=1; $i < $k; $i++) {
    push(@k_factors,$i);
  }
  print Dumper(\@k_factors);

  for (my $i =0; $i < scalar @$kgram_analysis; $i++) {
    my @split_kgram;
    if ($kgram_analysis->[$i] && $kgram_analysis->[$i] > 0) {
      my $kgram = base($i,10,$base);
      @split_kgram = split(//,$kgram);
    }
    #@split_kgram = qw(1 0 1 0);
  #}
    #print 'scalar kgram ',scalar @split_kgram,"\n";
  #print Dumper(\@split_kgram);

    #kfactor_way($k,\@split_kgram,\@k_factors);

    simpler_way($k,\@split_kgram);

    #stupid_way($k,\@split_kgram);

    my $dna_kgram;
    for (my $b = 0; $b < scalar @split_kgram; $b++) {
      #print ($bases[$split_kgram[$b]],"\n");
      #print "$dna_kgram\n";
      $dna_kgram .= $bases[$split_kgram[$b]];
    }
    if ($kgram_analysis->[$i]) {
      $conv_hash{$dna_kgram} = $kgram_analysis->[$i];
    }
  }
  #print Dumper(\%conv_hash);
  return(\%conv_hash);
}

sub simpler_way {

  my($k,$split_kgram) = @_;
  for (scalar @$split_kgram .. $k - 1) {
    unshift(@$split_kgram, 0);
  }
  return;
}

sub kfactor_way {

  my ($k,$split_kgram,$k_factors) = @_;
  #print("Scalar kgram: ",scalar @$split_kgram,"\n");
  #print Dumper($k_factors);
  if (scalar @$split_kgram == 0) {
    my $count =0;
    for (my $i=0; $i < $k ; $i++) {
      unshift(@$split_kgram, 0);
    }
  }
  my $count =0;
  if (scalar @$split_kgram > 0) {
    for my $kf (@$k_factors) {
      for (my $j=1; $j < $kf ;$j++) {
	unshift(@$split_kgram, 0);
      }
    }
  }
  return;
}

sub base4ToDec {
  #Implementation of: htun = (h * n2) + (t * n1) + (u * n0)
  my ($base4,$exps) = @_;
  my $dec_index;
  for (my $j = 0; $j < scalar @$base4; $j++) {
    $dec_index += $base4->[$j] * 4**($exps->[$j]);;
  }
  return($dec_index);
}

sub base5ToDec {
  #Implementation of: htun = (h * n3) + (h * n2) + (t * n1) + (u * n0)
  my ($base5,$exps) = @_;
  my $dec_index;
  for (my $j = 0; $j < scalar @$base5; $j++) {
    $dec_index += $base5->[$j] * 5**($exps->[$j]);;
  }
  return($dec_index);
}

sub base {
  my ($number, $inbase, $outbase) = @_;
  my ($realnum, $output, $i, $digit);
  # Convert the number (which might have letters) to lowercase.
  $number = lc($number);
  # Return undef (or an empty list) if the base is too weird.
  return if $inbase > 36 or $outbase > 36 or
    $inbase < 2 or $outbase < 2;
  # Convert $number from base $inbase to base 10.
  for $digit (reverse split(//, $number)) {
    $digit = ord($digit) - 87 if ord($digit) > 96;
    return if $digit >= $inbase;
    $realnum += $digit * ($inbase ** $i++);
  }
  # Convert the number from base 10 to $outbase.
  # logbase() is defined below.
  for ($i = int(logbase($realnum, $outbase)); $i >= 0; $i--) {
    $digit = int($realnum / ($outbase ** $i));
    $realnum -= $digit * ($outbase ** $i);
    $digit = ord($digit + 49) if ord($digit) > 57;
    if ($digit) {
      $output .= $digit;
    }
  }
  return $output;
}

sub logbase {
  my ($number, $base) = @_;
  return if $number <= 0 or $base <= 0 or $base == 1;
  return log ($number) / log($base);
}

sub randomBaseAssignment {

  my ($range,$genetic_code,$total_number_of_bases) = @_;
  my %base_hash;
  #print "RANGE: $range\n";
  #print Dumper($genetic_code);

  my @random_sequence;

  for (my $j = 0; $j < $total_number_of_bases; $j++) {
    my $random_integer;
    my $base;
    my $base_int;
    my @rands;

    @$genetic_code = shuffle(@$genetic_code);
    for (my $i = 0; $i < 4; $i++) {
      $random_integer = sprintf("%.0f", rand($range));
      push(@rands,$random_integer);
    }
    ($base,$base_int) = randomizedRand($genetic_code,\@rands);
    if (defined $base && $base ne q{}) {
      push(@random_sequence,$base);
    }
  }

  print "@random_sequence\n";
  #exit;
  return(\@random_sequence,\%base_hash);
}

sub randomizedRand {

  my ($bases,$rands) = @_;
  for (my $i = 0; $i < @$bases; $i++) {
    if ($rands->[$i] == 0) {
	$rands->[$i] = 1;
    }
    if ($rands->[$i] % 5 == 0 && $rands->[$i] % 3 == 0) {
      return ($bases->[$i],$rands->[$i]);
    }
    if ($rands->[$i] % 5 == 0 && $rands->[$i] % 2 == 0) {
      return ($bases->[$i],$rands->[$i]);
    }
    if ($rands->[$i] % 5 == 0 && $rands->[$i] % 2 != 0) {
      return ($bases->[$i],$rands->[$i]);
    }
    if ($rands->[$i] % 5 == 0 && $rands->[$i] % 3 != 0) {
      return ($bases->[$i],$rands->[$i]);
    }
    if ($rands->[$i] % 5 != 0 && $rands->[$i] % 3 == 0) {
      return ($bases->[$i],$rands->[$i]);
    }
    if ($rands->[$i] % 5 != 0 && $rands->[$i] % 2 == 0) {
      return ($bases->[$i],$rands->[$i]);
    }
    if ($rands->[$i] % $rands->[$i] == 0 && $rands->[$i] % 5 != 0 && $rands->[$i] % 3 != 0 && $rands->[$i] % 2 != 0) {
      return ($bases->[$i],$rands->[$i]);
    }
    if ($rands->[$i] == 1) {
      return ($bases->[$i],$rands->[$i]);
    }
  }
}

sub sendMail {

  my ($mode,$k,$count_files) = @_;

  my $to_field = 'js21@sanger.ac.uk';

  eval {
    my $msg = MIME::Lite->new(
			      To       => $to_field,
			      Subject  => "The DNA analysis of the $count_files files was  done for kgram: $k",
			      Type     => 'text/plain',
			      Encoding => 'base64',
			      Data     => "Number of file: $count_files/nAnalysis Mode: $mode\nKgram: $k\n"
			     );
    $msg->send;
  };
  warn() if $@;
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


sub stupid_way {

  my ($k,$split_kgram) = @_;

  if ($k == 3) {
    if (scalar @$split_kgram == 0) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    }
    if (scalar @$split_kgram == 1) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 2) {
      unshift(@$split_kgram, 0);
    }
  } elsif ($k == 4) {
    if (scalar @$split_kgram == 0) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    }
    if (scalar @$split_kgram == 1) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 2) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 3) {
      unshift(@$split_kgram, 0);
    }
  } elsif ($k == 5) {
    if (scalar @$split_kgram == 0) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    }
    if (scalar @$split_kgram == 1) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 2) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 3) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 4) {
      unshift(@$split_kgram, 0);
    }
  } elsif ($k == 6) {
    if (scalar @$split_kgram == 0) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    }
    if (scalar @$split_kgram == 1) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 2) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 3) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 4) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 5) {
      unshift(@$split_kgram, 0);
    }
  } elsif ($k == 7) {
    if (scalar @$split_kgram == 0) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    }
    if (scalar @$split_kgram == 1) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 2) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 3) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 4) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 5) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 6) {
      unshift(@$split_kgram, 0);
    }
  } elsif ($k == 8) {
    if (scalar @$split_kgram == 0) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    }
    if (scalar @$split_kgram == 1) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 2) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 3) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 4) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 5) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 6) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 7) {
      unshift(@$split_kgram, 0);
    }
  } elsif ($k == 9) {
    if (scalar @$split_kgram == 0) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    }
    if (scalar @$split_kgram == 1) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 2) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 3) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 4) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 5) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 6) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 7) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 8) {
      unshift(@$split_kgram, 0);
    }
  } elsif ($k == 10) {
    if (scalar @$split_kgram == 0) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    }
    if (scalar @$split_kgram == 1) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 2) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 3) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 4) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 5) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 6) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 7) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 8) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 9) {
      unshift(@$split_kgram, 0);
    }
  } elsif ($k == 11) {
    if (scalar @$split_kgram == 0) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    }
    if (scalar @$split_kgram == 1) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 2) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 3) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 4) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 5) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 6) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 7) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 8) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 9) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    }
    elsif (scalar @$split_kgram == 10) {
      unshift(@$split_kgram, 0);
    }
  }
  elsif ($k == 11) {
    if (scalar @$split_kgram == 1) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 2) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 3) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 4) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 5) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 6) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 7) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 8) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 9) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 10) {
      unshift(@$split_kgram, 0);
      unshift(@$split_kgram, 0);
    } elsif (scalar @$split_kgram == 11) {
      unshift(@$split_kgram, 0);
    }
  }
  return;
}
