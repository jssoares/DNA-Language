package DNA::Language;

use Moose;
use DNA::Language::Utils::FastaParser;
use DNA::Language::Utils::SimpleBaseCounter;
use DNA::Language::Utils::DNAcryptic;

has 'fasta_file' => ( is => 'rw', isa => 'Str', required => 1 );
has 'rand_range' => ( is => 'rw', isa => 'Str', required => 1 );
has 'tnb' => ( is => 'rw', isa => 'Int', required => 1 );
has 'mode' => ( is => 'rw', isa => 'Str', required => 1 );
has 'k' => ( is => 'rw', isa => 'Int', , required => 1 );
has 'host' => ( is => 'rw', isa => 'Str', required => 1 );
has 'branch' => ( is => 'rw', isa => 'Str', required => 1 );
has 'rmode' => ( is => 'rw', isa => 'Str', required => 1 );
has 'test' => ( is => 'rw', isa => 'Bool', required => 1 );

has 'genetic_code' => ( is => 'rw', isa => 'ArrayRef', lazy => 1, builder => '_build_genetic_code' );
has 'bio_seq' => ( is => 'rw', isa => 'Bio::SeqIO', lazy => 1, builder => '_build_bio_seq' );
has 'basic_stats' => ( is => 'rw', isa => 'HashRef', lazy => 1, builder => '_build_basic_stats' );
has 'seq_array' => ( is => 'rw', isa => 'ArrayRef', lazy => 1, builder => '_build_seq_array' );
has 'dna_cryptic' => ( is => 'rw', isa => 'Maybe[DNA::Language::Utils::DNAcryptic]', lazy => 1, builder => '_build_dna_cryptic' );

sub _build_genetic_code {

  my ($self) = @_;
  my @genetic_code = qw{a t c g};  
  return(\@genetic_code);
}

sub _build_bio_seq {

  my ($self) = @_;
  return(DNA::Language::Utils::FastaParser->new(fasta_file => $self->fasta_file)->bio_seq);
}

sub _build_basic_stats {

  my ($self) = @_;
  my %basic_stats;
  while (my $seq_obj = $self->bio_seq->next_seq){
    my $monomers_freq = DNA::Language::Utils::SimpleBaseCounter->new(seq => $seq_obj)->monomers_freq;
    $basic_stats{$seq_obj->display_id} = $monomers_freq;
  }
  return(\%basic_stats);
}

sub _build_seq_array {

  my ($self) = @_;
  my @seq_array;
  while (my $seq_obj = $self->bio_seq->next_seq){
    my $seq_string = $seq_obj->seq;
    @seq_array = split(//,$seq_string);
  }
  #print Dumper(\@seq_array);
  return(\@seq_array);
}

sub _build_dna_cryptic {
  my ($self) = @_;
  my $dna_cryptic = DNA::Language::Utils::DNAcryptic->new(
							  k => $self->k,
							  mode => $self->mode,
							  sequence => $self->seq_array
							 );
  $dna_cryptic->kgrams;
  $dna_cryptic->digit_to_dna;
  return($dna_cryptic);
}

sub run {

  my($self) = @_;
  $self->genetic_code;
  $self->bio_seq;
  #$self->basic_stats;
  $self->seq_array;
  $self->dna_cryptic;


}




no Moose;
__PACKAGE__->meta->make_immutable;
1;
