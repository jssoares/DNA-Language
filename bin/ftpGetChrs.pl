#!/software/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Storable;
use List::MoreUtils qw(uniq);
use List::Util qw(shuffle);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Net::FTP;
use MIME::Lite;
use Switch;
use Carp;


local $| = 1;       #setautoflush, outputs immediately instead of buffering

main();

sub main {

  my %opts;

  GetOptions( \%opts,
	      'help=s' => \my $help,
	      'h=s' => \my $host,    #Options: 'ensembl' or 'genomes'
	      'b=s' => \my $branch,     #Options: 'plants' or 'fungi' or 'protist' or 'metazoa' or 'bacteria'
	      'm=s' => \my $mode,      #Options: 'ro', 'rw'
	      't=i' => \my $test,   #Options: 1 (run the test sub_routine) or 0 (don't run it)
	    ) or pod2usage(2);

pod2usage(2) if ($help);

  my ($options_status,$ftp_host) = checkOptions($host,$branch,$mode,$test);

  if ($options_status == 1) {
    my ($ftp,$main_dir,$species) = getSpeciesEnsemblFTP($ftp_host);
    if ($test == 1) {
      print "Running test sub_routine\n";
      testFTP($ftp,$main_dir,$mode,$species,$branch);
    }
    elsif ($test == 0) {
      print "Won't run the test sub_routine\n";
    }
  }
  elsif ($options_status == 0) {
    die "Should have died before already. What the hell is going on\n"; 
  }

}

sub getSpeciesEnsemblFTP {

  my ($ftp_host,$branch) = @_;
  my $ftp = Net::FTP->new($ftp_host, Debug => 0)
    or die "Cannot connect to $ftp_host: $@";

  $ftp->login("anonymous",'-anonymous@')
    or die "Cannot login ", $ftp->message;

  my $main_dir;
  my @species;

  unless ($ftp_host =~ m/genomes/) {

    $main_dir = "/pub/current_fasta/";

    $ftp->cwd($main_dir)
	or die "Cannot change working directory ", $ftp->message;

    @species = $ftp->ls();
  }

  else {

    if ($branch) {

      my $branch_dir = "/pub/current/$branch/fasta/";

      $ftp->cwd($branch_dir)
	or die "Cannot change working directory ", $ftp->message;

      unless ($branch eq 'bacteria') {

	$main_dir = $branch_dir;

	@species = $ftp->ls();
      }

      else {

	my @collections = $ftp->ls();

	for my $collection(@collections) {

	  $main_dir = $branch_dir . $collection . '/';

	  $ftp->cwd($main_dir)
	    or die "Cannot change working directory ", $ftp->message;

	  @species = $ftp->ls();
	}
      }
    }
  }

  return($ftp,$main_dir,\@species);

}

sub testFTP {

  my ($ftp,$main_dir,$mode,$species,$branch)  = @_;

  if ($mode ne q{}) {
    if ($mode eq 'rw') {
      "Running under $mode mode. Will download all the files I find\n";
    }
    elsif ($mode eq 'ro') {
      "Running under $mode mode. Will not download any of the files I find\n";
    }
  }

  my %file_counter;

  for my $organism(@$species) {
    unless ($organism =~ m/^ancestral/) {
      my $work_dir = $main_dir . $organism . '/dna';
      print "$work_dir\n";
      $ftp->cwd($work_dir)
	or die "Cannot change working directory ", $ftp->message;
      
      print "SPECIES: $organism\n";

      my @file_list = $ftp->ls();

      my $local_dir;
      if ($mode eq 'rw') {
	print "Creating local directories\n";
	$local_dir = createLocalDirectory($organism);
      }

      my $count = 0;
      unless ($branch eq 'bacteria') {
	for (my $i = 0; $i < scalar @file_list; $i++) {
	  if ($file_list[$i] =~ m/\.dna\.chromosome\.\d{1,}/ || $file_list[$i] =~ m/\.dna\.chromosome\.[XxYyIiVv|MT]/) {

	    my $local_file_path;

	    if ($mode eq 'rw') {
	      $local_file_path = $local_dir . '/' . $file_list[$i];
	    }

	    print "getting $file_list[$i]\n";
	    if ($mode eq 'rw') {
	      my $status = $ftp->get($file_list[$i],$local_file_path);
	      if ($status ne $local_file_path) {
		print "failed to get File: $file_list[$i]\n";
	      }
	    }
	    $count++;
	  }
	}
	$file_counter{$organism} = $count;
      }
      else {
	for (my $i = 0; $i < scalar @file_list; $i++) {
	  if ($file_list[$i] =~ m/\.dna\.chromosome\.Chromosome/) {

	    my $local_file_path;
	    if ($mode eq 'rw') {
	      $local_file_path = $local_dir . '/' . $file_list[$i];
	    }

	    print "getting $file_list[$i]\n";
	    if ($mode eq 'rw') {
	      my $status = $ftp->get($file_list[$i],$local_file_path);
	      if ($status ne $local_file_path) {
		print "failed to get File: $file_list[$i]\n";
	      }
	    }
	    $count++;
	  }
	}
	$file_counter{$organism} = $count;
      }
    }
  }

  for my $organism(keys %file_counter) {
    if ($file_counter{$organism} == 0) {
      my $work_dir = $main_dir . $organism . '/dna/';
      print "$work_dir\n";
      $ftp->cwd($work_dir)
	or die "Cannot change working directory ", $ftp->message;

      my $local_dir;
      if ($mode eq 'rw') {
	print "Creating local directories\n";
	$local_dir = createLocalDirectory($organism);
      }
  
      print "SPECIES: $organism\n";
      my @file_list = $ftp->ls();
      for (my $i = 0; $i < scalar @file_list; $i++) {

	if ($file_list[$i] !~ m/\.dna\.chromosome\.\d{1,}/ && $file_list[$i] !~ m/\.dna\.chromosome\.[XxYyIiVv|MT]/ && $file_list[$i] !~ m/\.dna\.chromosome\.Chromosome/ && $file_list[$i] =~ m/\.dna\.toplevel/) {

	  my $local_file_path;

	  if ($mode eq 'rw') {

	    $local_file_path = $local_dir . '/' . $file_list[$i];
	  }

	  print "getting $file_list[$i]\n";

	  if ($mode eq 'rw') {

	    my $status = $ftp->get($file_list[$i],$local_file_path);

	    if ($status ne $local_file_path) {

	      print "failed to get File: $file_list[$i]\n";
	    } 
	  }
	}
      }
    }
  }
}

sub createLocalDirectory {

  my $organism = shift;

  my $local_dir = "../new_dna_language_sequences/$organism";
  if (! -e $local_dir) {
    mkdir $local_dir
      or die "Cannot create dir $local_dir: $!";
  }
  return ($local_dir);
}

sub checkOptions {

  my ($host,$branch,$mode,$test) = @_;

  my $genome_branches = allowedBranches();
  my $modes = allowedModes();
  my $options_status = 0;

  my $ftp_host;

  if ($host) {
    if ($branch) {
      if ($mode) {
	if ($test) {
	  if ($test == 1 || $test == 0) {
	    if ($modes->{$mode}) {
	      if ($host eq 'ensembl') {
		$ftp_host = "ftp.ensembl.org";
		$options_status = 1;
	      } 
	      elsif ($host eq 'genomes') {
		if ($genome_branches->{$branch}) {
		  $ftp_host = "ftp.ensemblgenomes.org";
		  $options_status = 1;
		}
		else {
		  die "Unsupported branch $branch. Bye\n";
		}
	      }
	      else {
		die "Unknown host: $host specified. Bye\n";
	      }
	    }
	  }
	  else {
	    die "Unknown test flag: $test specified. Bye\n";
	  }
	}
	else {
	  die "Test option not specified. Bye\n";
	}
      }
      else {
	die "Mode option not specified. Bye\n";
      }
    }
    else {
      die "Branch option not specified. Bye\n";
    }
  }
  else {
    die "Host option not specified. Bye\n";
  }

  return($options_status,$ftp_host);
}


sub allowedBranches {

  my %genome_branches;
  $genome_branches{plants} = 1;
  $genome_branches{fungi} = 1;
  $genome_branches{protist} = 1;
  $genome_branches{metazoa} = 1;
  $genome_branches{bacteria} = 1;

  return(\%genome_branches);
}

sub allowedModes {

  my %modes;
  $modes{rw} = 1;
  $modes{ro} = 1;
  return(\%modes);
}


__END__

=head1 NAME

ftpGetChrs.pl.pl

=head1 SYNOPSIS

perl ftpGetChrs.pl -h [ 'ensembl' || 'genomes' ] -b [ 'plants' || 'fungi' || 'protist' || 'metazoa' || 'bacteria' ] -m [ 'ro' || 'rw' ] -t [ 1 || 0 ]

=head1 AUTHOR

Jorge Soares , E<lt>js21@sanger.ac.ukE<gt>

__END__
=head1 COPYRIGHT AND LICENSE

Totally free to use wherever or whenever you want. This is my personal work

=cut
