package DNA::Language::CommandLine::LanguageAnalysis;

# ABSTRACT: 

=head1 SYNOPSIS


=cut

use Moose;
use Getopt::Long qw(GetOptionsFromArray);
use DNA::Language;

has 'args'        => ( is => 'ro', isa => 'ArrayRef', required => 1 );
has 'script_name' => ( is => 'ro', isa => 'Str',      required => 1 );
has 'help'        => ( is => 'rw', isa => 'Bool',     default  => 0 );

has 'fasta_file' => ( is => 'rw', isa => 'Str' );
has 'rand_range' => ( is => 'rw', isa => 'Str' );
has 'total_number_of_bases' => ( is => 'rw', isa => 'Int' );
has 'mode' => ( is => 'rw', isa => 'Str', lazy => 1, default => q(real) );
has 'k' => ( is => 'rw', isa => 'Int', lazy => 1, default => 3 );
has 'host' => ( is => 'rw', isa => 'Str', lazy => 1, default => q() );
has 'branch' => ( is => 'rw', isa => 'Str', lazy => 1, default => q() );
has 'rmode' => ( is => 'rw', isa => 'Str', lazy => 1, default => q(ro) );
has 'test' => ( is => 'rw', isa => 'Bool', lazy => 1, default => 0 );

sub BUILD {

    my ($self) = @_;

    my ( $fasta_file, $mode, $k, $rand_range, $total_number_of_bases, $host, $branch,
        $rmode, $test, $help );

    GetOptionsFromArray(
			$self->args,
			'f|fasta:s' => \$fasta_file,
        'rg|rand_range:s' => \$rand_range,
        'n|tnb=s'         => \$total_number_of_bases,
        'm|mode:s'        => \$mode,
        'k|kgramm:s'      => \$k,
        'ho|host=s'       => \$host,
        'b|branch=s'      => \$branch,
        'r|rmode=s'       => \$rmode,
        't|test'          => \$test,
        'h|help'          => \$help,
    );

    $self->mode($mode)             if ( defined($mode) );
    $self->k($k)                   if ( defined($k) );
    $self->rand_range($rand_range) if ( defined($rand_range) );
    $self->total_number_of_bases($total_number_of_bases)
      if ( defined($total_number_of_bases) );
    $self->host($host)     if ( defined($host) );
    $self->branch($branch) if ( defined($branch) );
    $self->rmode($rmode)   if ( defined($rmode) );
    $self->test($test)     if ( defined($test) );

}

sub run {

  my($self) = @_;
  
  ( $self->rand_range && $self->total_number_of_bases && $self->fasta_file ) or die <<USAGE;

Usage:
  -f|fasta                   <Options: fasta file to analyse >
  -rg|rand_range             <Options: range of the universe of random numbers to choose form when running under random mode. From (4 to infinity really)>
  -tnb|total_number_of_bases <Options: to infinity and beyond>
  -m|mode                    <Options: 'random','test' or 'real'>
  -k|kgram                   <Options: 3 4 5 6 7 8 9 10 11>
  -h|host                    <Options: 'ensembl' or 'genomes'>
  -b|branch                  <Options: 'plants' or 'fungi' or 'protist' or 'metazoa' or 'bacteria'>
  -r|rmode,                  <Options: 'ro', 'rw'>
  -t|test                    <Options: 1 (run the test sub_routine) or 0 (don't run it)>
  -?|help                    <Print this message>

USAGE

  my $lanalysis = DNA::Language->new(
				     fasta_file => $self->fasta_file,
				     mode => $self->mode,
				     k => $self->k,
				     rand_range => $self->rand_range,
				     tnb => $self->total_number_of_bases,
				     host => $self->host,
				     branch => $self->branch,
				     rmode =>$self->rmode,
				     test => $self->test
				    );
  $lanalysis->run;
  
}



no Moose;
__PACKAGE__->meta->make_immutable;
1;
