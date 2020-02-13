#!/usr/bin/env perl
# EPN, Sat Jan  5 07:20:01 2019
use warnings;
use strict;
use Getopt::Long;

my $usage;
$usage  = "perl list-get-tax-efetch.pl [OPTIONS] <file with list of sequences>\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t-v:          be verbose, output commands as they are run\n";
$usage .= "\t\t--unknownok: output UNKNOWN-TAXONOMY (instead of failing) if taxonomy of a seq is unknown because it is not in nuccore (wgs, maybe?)\n";

my $opt_unknownok = undef;
my $be_verbose = 0;
&GetOptions(
  "v"         => \$be_verbose,
  "unknownok" => \$opt_unknownok, 
    );

my $do_unknownok = 0;
if(defined $opt_unknownok) { $do_unknownok = 1; }

if(scalar(@ARGV) != 1) { die $usage; }
my ($list_file) = (@ARGV);

if(! -s $list_file) { die "ERROR fasta file $list_file does not exist or is empty"; }


# read list file
my @seqname_list_str_A = ();
my $seqname_list_str = "";
my $list_size = 0;
my $block_size = 100;
my $max_width = 0;
open(LIST, $list_file) || die "ERROR unable to open $list_file for reading";
my $line;
while($line = <LIST>) { 
  chomp $line;
  my $orig_seqname = $line;
  if(length($orig_seqname) > $max_width) { $max_width = length($orig_seqname); }
  my $acc = undef;
  if($orig_seqname =~ /^gi\|\d+\|\S+\|(\S+)\.\d+\|.*/) { 
    $acc = $1;
  }
  elsif($orig_seqname =~ /^(\S+)\.\d+/) { 
    $acc = $1;
  }
  else { 
    die "ERROR unable to parse sequence name $orig_seqname";
  }
  if($seqname_list_str ne "") { $seqname_list_str .= ","; }
  $seqname_list_str .= $acc;
  $list_size++;
  if($list_size == $block_size) { 
    push(@seqname_list_str_A, $seqname_list_str);
    $seqname_list_str = "";
    $list_size = 0;
  }
}
close(LIST);
if($list_size > 0) { 
    push(@seqname_list_str_A, $seqname_list_str);
}  

# run efetch to get tax strings
my $tax_file = $list_file . ".tax";

my $nlists = scalar(@seqname_list_str_A);

my $eutils_cmd = "";
if(-e $tax_file) { unlink $tax_file; }
foreach $seqname_list_str (@seqname_list_str_A) { 
  $eutils_cmd = "esearch -db nuccore -query \"$seqname_list_str\" | efetch -format gpc | xtract -insd INSDSeq_taxonomy >> $tax_file";
  run_command($eutils_cmd, $be_verbose);
}

# parse efetch output
open(TAX, $tax_file) || die "ERROR unable to open $tax_file for reading";
my %tax_H = ();
while($line = <TAX>) { 
  chomp $line;
  if($line =~ m/\w/) { 
    my @el_A = split("\t", $line);
    if(scalar(@el_A) != 2) { die "ERROR failed to read exactly 2 tab delimited tokens in efetch tax output file line: $line"; }
    my ($accver, $tax_str) = ($el_A[0], $el_A[1]);
    my $acc = $accver;
    $acc =~ s/\.\d+$//;
    $tax_H{$acc} = $tax_str;
  }
}
close(TAX);

# go through the list file and output tax strings in the headers
open(LIST, $list_file) || die "ERROR unable to open $list_file for reading";
while($line = <LIST>) { 
  chomp $line;
  my $orig_seqname = $line;
  my $acc = undef;
  if($orig_seqname =~ /^gi\|\d+\|\S+\|(\S+)\.\d+\|.*/) { 
    $acc = $1;
  }
  elsif($orig_seqname =~ /^(\S+)\.\d+/) { 
    $acc = $1;
  }
  else { 
    die "ERROR unable to parse sequence name $orig_seqname";
  }
  my $cur_tax = undef;
  if(exists $tax_H{$acc}) {
    $cur_tax = $tax_H{$acc};
  }
  else { # maybe from the WGS database? Not in nuccore at least... can get around this by using srcchk, but then you don't get the full tax strings
    if($do_unknownok) { 
      $cur_tax = "UNKNOWN-TAXONOMY";
    }
    else { 
      die "ERROR, no tax info for sequence $orig_seqname ($acc), use --unknownok to allow this";
    }
  }
  printf("%-*s %s\n", $max_width, $orig_seqname, $cur_tax);
}
close(LIST);

exit 0;

#################################################################
# Subroutine:  run_command()
# Incept:      EPN, Mon Dec 19 10:43:45 2016
#
# Purpose:     Runs a command using system() and exits in error 
#              if the command fails.
##
# Arguments:
#   $cmd:         command to run, with a "system" command;
#   $be_verbose:  '1' to output command to stdout, '0' not to
#
# Returns:    amount of time the command took, in seconds
#
# Dies:       if $cmd fails
#################################################################
sub run_command {
  my $sub_name = "run_command()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd, $be_verbose) = @_;

  if($be_verbose) { 
    print ("Running cmd: $cmd\n"); 
  }

  system($cmd);

  if($? != 0) { 
    die "ERROR in $sub_name, the following command failed:\n$cmd\n";
  }

  return;
}
