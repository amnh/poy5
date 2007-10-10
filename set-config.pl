#!/usr/bin/perl -w

# This file should be called as:

# ./set-config.pl [--parallel] [--ncurses] [--graphics] [--verbose]

# By default, each of the three options (parallel, ncurses, graphics) will be
# turned off.  Specify these options to explicitly turn them on.  Specifying
# --verbose will display the contents of the generated config file.


my $parallel = 0;
my $curses = 0;
my $graphics = 0;
my $verbose = 0;

foreach (@ARGV) {
  $parallel = 1 if m/^-*parallel$/;
  $curses = 1   if m/^-*n?curses$/;
  $graphics = 1 if m/^-*graphics$/;
  $verbose = 1  if m/^-*v(erbose)?$/;
}

open(CONFIG, "<config");
my @config = <CONFIG>;
close(CONFIG);

@config = grep !/USEPARALLEL|USENCURSES/, @config;

unshift
  @config,
  ("export USEPARALLEL = " . ($parallel ? "true\n" : "false\n")),
  ("export USENCURSES = "  . ($curses   ? "true\n" : "false\n"));
#   ("export USEGRAPHICS = " . ($graphics ? "true\n" : "false\n"));

if ($verbose) {
  foreach (@config) {
    print $_;
  }
}

open(CONFIG, ">config");
foreach (@config) {
  print CONFIG $_;
}
close(CONFIG);
