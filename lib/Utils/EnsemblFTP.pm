package EnsemblFTP;

use strict;
use warnings;
use Data::Dumper;
use Net::FTP;
use Carp;


our $VERSION = "1.00";
our ($AUTOLOAD, %ok_field);

# Get/set methods to be handled by AUTOLOAD
for my $attr (
	      qw(
		  ftp_obj
		  main_dir
		  species
	       )
	     ) { $ok_field{$attr}++; }

sub new {

  my($class, %args) = @_;

  #The EnsemblFTP object needs a set of arguments at its initialisation.
  #These are essential for setting the attributes bellow. (Hence the verificaton and default values)

  my $self = bless({}, $class);
  my $host = exists $args{host} ? $args{host} : 'no_host';            #the ftp host
  my $branch = exists $args{branch} ? $args{branch} : 'no_branch';    #the branch. Like plants, fungi, protists, metazoa or bacteria. Flag only needed if connecting to genomes
  my $rmode = exists $args{rmode} ? $args{rmode} : 'no_mode';            #the mode. Like Read/Write (rw) or Read Only (ro)
  my $test = exists $args{test} ? $args{test} : '2';                  #the test flag. To run the code in test mode. Needs to be set in the caller.

 TODO:
  #I need to have a look at the implementation of the setting of the test variable.
  #Maybe doesn't really belong here.

  my $ftp_host;
  #ftp_host needs to be either 'ensembl' or 'genomes'
  if ($args{host} eq 'ensembl') {
    $ftp_host = "ftp.ensembl.org";
  } 
  elsif ($args{host} eq 'genomes') {
    $ftp_host = "ftp.ensemblgenomes.org";
  }
  else {
    croak "No known host\n";
  }

  #Setting objects initial attributes
  $self->ftp_host($ftp_host);
  $self->branch($branch);
  $self->rmode($rmode);
  $self->test($test);

  return $self;

}

sub getSpeciesEnsemblFTP {

  #Object method to set ftp_obj, main_dir and species.
  #This method navigates through two ensembl ftp sites.
  #It sets up information about the fasta file locations
  #for the species specified in the command line options.
  my $self = shift;

  #ftp_host should be either 'ensembl' or 'genomes'
  my $ftp = Net::FTP->new($self->ftp_host, Debug => 0)
    or croak "Cannot connect to " . $self->ftp_host . ": $@";

  #As specified in the ensembl ftp site docs
  $ftp->login("anonymous",'-anonymous@')
    or croak "Cannot login ", $ftp->message;

  my $main_dir;
  my @species;

  #Path navigation through the two FTP sites is different.
  unless ($self->ftp_host =~ m/genomes/) {

    #Moving to the current version of the fasta files in the ensembl ftp site
    $main_dir = "/pub/current_fasta/";

    $ftp->cwd($main_dir)
	or croak "Cannot change working directory ", $ftp->message;

    @species = $ftp->ls();      #Gathering the list of available species
  }

  else {
    #Checking once again if -branch is set, but it should never be caught at this point.
    #Should be upstream.
    if ($self->branch) {
      my $branch_dir = "/pub/current/" . $self->branch . "/fasta/";
      #Moving to the current version of the fasta files in the genomes ftp site
      $ftp->cwd($branch_dir)
	or croak "Cannot change working directory ", $ftp->message;

      #Bacteria are then sub-divided into collections
      unless ($self->branch eq 'bacteria') {
	$main_dir = $branch_dir;
	@species = $ftp->ls();      #Gathering the list of available species
      }
      else {
	my @collections = $ftp->ls();
	for my $collection(@collections) {
	  $main_dir = $branch_dir . $collection . '/';
	  $ftp->cwd($main_dir)
	    or croak "Cannot change working directory ", $ftp->message;

	  @species = $ftp->ls();      #Gathering the list of available species
	}
      }
    }
    else {
      croak "No branch information. I came from getSpeciesEnsemblFTP"
    }
  }

  #Setting attributes
  $self->ftp_obj($ftp);
  $self->main_dir($main_dir);
  $self->species(\@species);
  return $self;

  #Legacy non OO code
  #return($ftp,$main_dir,\@species);
}

sub testFTP {

  my $self = shift;

  #This object method checks once again if the flags test and mode have been set.
  #It's probably redundant, but will only remove it if it impacts on performance.
  #We shall see.
  unless (!defined $self->test) {
    unless (!defined $self->rmode) {
      if ($self->rmode ne q{}) {

	if ($self->rmode eq 'rw') {
	  "Running under " . $self->rmode . " rmode. Will download all the files I find\n";
	} elsif ($self->rmode eq 'ro') {
	  "Running under " . $self->rmode . " rmode. Will not download any of the files I find\n";
	}

      }
      my %file_counter;

      #The species attribute that was set by getEnsemblSpecies
      for my $organism (@{ $self->species }) {

	#I think archaea lives in ancestral, but not sure. Need to check this
	unless ($organism =~ m/^ancestral/) {

	  #Moving to the organism's dir in the ftp site
	  my $work_dir = $self->main_dir . $organism . '/dna';
	  print "$work_dir\n";
	  $self->ftp_obj->cwd($work_dir)
	    or croak "Cannot change working directory ", $self->ftp_obj->message;

	  print "SPECIES: $organism\n";

	  #Gathering a list of files present in the dir
	  my @file_list = $self->ftp_obj->ls();

	  #At this point, if run under rw mode,
	  #a directory structure is created for each organism, on the local machine
	  my $local_dir;
	  if ($self->rmode eq 'rw') {
	    print "Creating local directories\n";
	    $local_dir = createLocalDirectory($organism);
	  }

	  #FTP navigation for the bacteria branch is different
	  my $count = 0;
	  if (defined $self->branch) {
	    unless ($self->branch eq 'bacteria') {

	      for (my $i = 0; $i < scalar @file_list; $i++) {
		#dirname is terminated by the chromosome number and by X | Y | V | MT (mitochondrial)
		#I don't know what V stands for anymore.
		if ($file_list[$i] =~ m/\.dna\.chromosome\.\d{1,}/ || $file_list[$i] =~ m/\.dna\.chromosome\.[XxYyIiVv|MT]/) {

		  my $local_file_path;

		  #Establishing the file path to write locally
		  if ($self->rmode eq 'rw') {
		    $local_file_path = $local_dir . '/' . $file_list[$i];
		  }

		  print "FTP get  $file_list[$i]\n";
		  if ($self->rmode eq 'rw') {
		    #Writing the file locally. And I think I'm handling the exception well, as well.
		    #But I need to check if the $status is at any point equal to the the $local_file_path

		    my $status = $self->ftp_obj->get($file_list[$i],$local_file_path);
		    if ($status ne $local_file_path) {
		      print "failed to get File: $file_list[$i]\n";
		    }

		  }
		  $count++;
		}
	      }
	      #Setting the file counter for each organism
	      $file_counter{$organism} = $count;

	    } else {
	      #Dealing with the bacteria FTP navigation and file writing

	      for (my $i = 0; $i < scalar @file_list; $i++) {
		if ($file_list[$i] =~ m/\.dna\.chromosome\.Chromosome/) {

		  my $local_file_path;
		  if ($self->rmode eq 'rw') {
		    $local_file_path = $local_dir . '/' . $file_list[$i];
		  }

		  print "FTP get $file_list[$i]\n";
		  if ($self->rmode eq 'rw') {
		    my $status = $self->ftp_obj->get($file_list[$i],$local_file_path);
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
      } #Matches       for my $organism (@{ $self->species }) {


      #Navigating through the FTP site for directories of  organisms that have a slightly different layout
      #I don't know why I do this here. Probably found out later?
      #Anyway will keep it here but if code needs optimisation this might be a good place to start
      for my $organism (keys %file_counter) {
	if ($file_counter{$organism} == 0) {
	  my $work_dir = $self->main_dir . $organism . '/dna/';
	  print "$work_dir\n";
	  $self->ftp_obj->cwd($work_dir)
	    or croak "Cannot change working directory ", $self->ftp_obj->message;

	  my $local_dir;
	  if ($self->rmode eq 'rw') {
	    print "Creating local directories\n";
	    $local_dir = createLocalDirectory($organism);
	  }

	  print "SPECIES: $organism\n";
	  my @file_list = $self->ftp_obj->ls();
	  for (my $i = 0; $i < scalar @file_list; $i++) {

	    #Going deeper in the FTP site into .toplevel
	    if ($file_list[$i] !~ m/\.dna\.chromosome\.\d{1,}/ && $file_list[$i] !~ m/\.dna\.chromosome\.[XxYyIiVv|MT]/ && $file_list[$i] !~ m/\.dna\.chromosome\.Chromosome/ && $file_list[$i] =~ m/\.dna\.toplevel/) {

	      my $local_file_path;

	      if ($self->rmode eq 'rw') {

		$local_file_path = $local_dir . '/' . $file_list[$i];
	      }

	      print "getting $file_list[$i]\n";

	      if ($self->rmode eq 'rw') {

		my $status = $self->ftp_obj->get($file_list[$i],$local_file_path);

		if ($status ne $local_file_path) {

		  print "failed to get File: $file_list[$i]\n";
		}
	      }
	    }
	  }
	}
      }
    } else {
      unless ($self->test == 0) {
	croak "-m [rmode] flag required if -t [test] flag \=\= " . $self->test . "\n";
      }
    }
  }
  else {
    croak "Can't run a test without defining -t flag\n";
  }
}

sub createLocalDirectory {

  my $organism = shift;

  my $local_dir = "../new_dna_language_sequences/$organism";
  if (! -e $local_dir) {
    mkdir $local_dir
      or croak "Cannot create dir $local_dir: $!";
  }
  return ($local_dir);
}

sub ftp_host {
  my $self = shift;
  if( @_ ) {
    my $ftp_host = shift;
    $self->{ftp_host} = $ftp_host;
  }
  return $self->{ftp_host};
}

sub branch {
  my $self = shift;
  if( @_ ) {
    my $branch = shift;
    $self->{branch} = $branch;
  }
  return $self->{branch};
}

sub rmode {
  my $self = shift;
  if( @_ ) {
    my $rmode = shift;
    $self->{rmode} = $rmode;
  }
  return $self->{rmode};
}

sub test {
  my $self = shift;
  if( @_ ) {
    my $test = shift;
    $self->{test} = $test;
  }
  return $self->{test};
}


sub AUTOLOAD {

  my $self = shift;
  my $attr = $AUTOLOAD;
  $attr =~ s/.*:://;

  return unless $attr =~ /[^A-Z]/;  #skip DESTROY and all-cap methods

  croak "invalid attribute method: ->$attr()"
    unless $ok_field{$attr};
  $self->{$attr} = shift if @_;

  return $self->{$attr};
}

1;
