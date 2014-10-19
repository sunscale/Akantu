#!/usr/bin/perl
use strict;
use feature "say";

sub mean {
    my (@table) = @_;
    my $sum = 0;
    for (@table) { $sum += $_; }
    return $sum/($#table + 1);
}

sub variance {
    my ($mean, @table) = @_;
    my $sum = 0;
    for (@table) { $sum += ($_ - $mean)**2; }
    return $sum/($#table + 1);
}

my $fh;
my @files = ();

my $filename = shift @ARGV;
my $tol = shift @ARGV;

open($fh, $filename) or die "Cannot open $filename";

my @file = <$fh>;
my @mesure = ();
my $nb_line = 0;
READ: for (@file) {
  $nb_line++;
  next READ if $nb_line == 1;
  my @line = split /,/;
  $line[$#line] =~ /([+-]?[0-9]+\.?[0-9]*([eE][+-]?[0-9]+)?)/;
  push @mesure, $1;
}

close($fh);

my $mean = mean(@mesure);
my $std_devi = sqrt(variance($mean, @mesure));

if (($std_devi/$mean <= $tol)) {
  say "Pass: Mean $mean (std devi $std_devi) -> ".($std_devi/$mean)." <= $tol";
  exit 0;
} else {
  say "FAIL: Mean $mean (std devi $std_devi) -> ".($std_devi/$mean)." > $tol";
  exit 1;
}
