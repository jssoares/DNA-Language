package DNA::Language::Utils::DNAcryptic;

use Moose;
use Carp;

has k => ( is => 'rw', isa => 'Int', lazy => 1,  default => q(3) );
has mode => ( is => 'rw', isa => 'Str', lazy => 1, default => q(no_mode) );
has sequence => ( is => 'rw', isa => 'ArrayRef', lazy=> 1, default => q(no_sequence) );
has kgrams => ( is => 'rw', isa => 'ArrayRef', lazy => 1, builder => '_build_kgrams');
has digit_to_dna => ( is => 'rw', isa => 'HashRef', lazy => 1, builder => '_build_digit_to_dna');


sub _build_kgrams {

  my ($self) = @_;
  my $k = $self->k;

  print "HELLO\n";
  
  my %base_hash;
  if ($self->mode eq 'test' || $self->mode eq 'real') {
    $base_hash{a} = 0;
    $base_hash{c} = 1;
    $base_hash{g} = 2;
    $base_hash{t} = 3;
    $base_hash{n} = 4;
  }
  elsif ($self->mode eq 'random') {
    $base_hash{a} = 0;
    $base_hash{c} = 1;
    $base_hash{g} = 2;
    $base_hash{t} = 3;
  }

  my $kminus1 = $k - 1;
  print "KMINUS1: $kminus1\n";
  my $filter;

  eval {
    $filter = scalar @{ $self->sequence } - $kminus1;
  };
  confess if $@;

  my @exps = (0..$kminus1);
  @exps = reverse(@exps);

  my @kgram;

  #foreach letter minus 1 in the sequence
  for (my $b = 0; $b < $filter; $b++) {

    my @base4;
    my @base5;

    #print "$b,";
    #Split words, convert them into base 5 and load them into @baseX array
    for(my $i = 0; $i < $k; $i++) {
      #print ($self->sequence->[$b + $i]);
      if ($self->mode eq 'test' || $self->mode eq 'real') {
	push(@base5,$base_hash{$self->sequence->[$b + $i]});
      }
    }

    #Convert base 5 words into base 10
    #Maybe I have to create several @kgram arrays to stop:
    #ERROR - Modification of non-creatable array value attempted, subscript...

    eval {
      my $i_dec;
      if ($self->mode eq 'test' || $self->mode eq 'real') {
	my $dec_index;
	for (my $j = 0; $j < scalar @base5; $j++) {
	  $dec_index += $base5[$j] * 5**($exps[$j]);;
	}

	if ($dec_index < 999999999) {
	$i_dec = _base5ToDec(\@base5,\@exps);
	#print "\t$dec_index\n";
	#increment the count on the index $dec_index of the @kgram
	$kgram[$i_dec]++;

	}
      }
    };
    warn if $@;
  }

  return(\@kgram);
}


sub _build_digit_to_dna {

  my ($self) = @_;
  my $k = $self->k;

  my (%digit_to_dna, @bases, $base);

  if ($self->mode eq 'test' || $self->mode eq 'real') {
    @bases = qw{a c g t n};
    $base = 5;
  }
  elsif ($self->mode eq 'random') {
    @bases = qw{a c g t};
    $base = 4;
  }

  my @split_kgram;
  if (defined $self->kgrams) {
    for (my $i = 0; $i < scalar @{ $self->kgrams }; $i++) {
      if ($self->kgrams->[$i] && $self->kgrams->[$i] > 0) {
	my $kgram = _base($i,10,$base);
	@split_kgram = split(//,$kgram);
	_simpler_way($k,\@split_kgram);
	#_stupid_way($k,\@split_kgram);
	my $dna_kgram;

	for (my $b = 0; $b < scalar @split_kgram; $b++) {
	  $dna_kgram .= $bases[$split_kgram[$b]];
	}
	if ($self->kgrams->[$i]) {
	  $digit_to_dna{$dna_kgram} = $self->kgrams->[$i];
	}
      }
    }
    return(\%digit_to_dna);
  }
}

sub _simpler_way {

  my($k,$split_kgram) = @_;
  for (scalar @$split_kgram .. $k - 1) {
    unshift(@$split_kgram, 0);
  }
  return;
}

sub _base4ToDec {
  #Implementation of: htun = (h * n2) + (t * n1) + (u * n0)
  my ($base4,$exps) = @_;
  my $dec_index;
  for (my $j = 0; $j < scalar @$base4; $j++) {
    $dec_index += $base4->[$j] * 4**($exps->[$j]);;
  }
  return($dec_index);
}

sub _base5ToDec {
  #Implementation of: htun = (h * n3) + (h * n2) + (t * n1) + (u * n0)
  my ($base5,$exps) = @_;
  my $dec_index;
  for (my $j = 0; $j < scalar @$base5; $j++) {
    $dec_index += $base5->[$j] * 5**($exps->[$j]);;
  }
  return($dec_index);
}

sub _base {
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
  eval {
    for ($i = int(_logbase($realnum, $outbase)); $i >= 0; $i--) {
      $digit = int($realnum / ($outbase ** $i));
      $realnum -= $digit * ($outbase ** $i);
      $digit = ord($digit + 49) if ord($digit) > 57;
      if (defined $digit && $digit ne q{}) {
	$output .= $digit;
      }
    }
  };
  confess if $@;

  return $output;
}

sub _logbase {
  my ($number, $base) = @_;
  return if $number <= 0 or $base <= 0 or $base == 1;
  return log ($number) / log($base);
}

#Only here for legacy reasons and for me to try to understand what I was doing wiht this.
sub _stupid_way {

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

no Moose;
__PACKAGE__->meta->make_immutable;
1;
