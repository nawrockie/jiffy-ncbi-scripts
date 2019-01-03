#!/usr/bin/env perl
# EPN, Tue Jan  1 06:32:21 2019
use warnings;
use strict;
use Getopt::Long;

my $usage;
$usage  = "perl fasta-add-tax-as-desc-efetch.pl [OPTIONS] <fasta file (cannot be read from stdin)>\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t--prepend:   add taxonomy string at beginning of existing seq description [default: replace current desc with tax string]\n";
$usage .= "\t\t--append:    add taxonomy string at end       of existing seq description [default: replace current desc with tax string]\n";
$usage .= "\t\t--unknownok: output UNKNOWN-TAXONOMY (instead of failing) if taxonomy of a seq is unknown because it is not in nuccore (wgs, maybe?)\n";

my $opt_prepend   = undef;
my $opt_append    = undef;
my $opt_unknownok = undef;
&GetOptions( "prepend"   => \$opt_prepend, 
             "append"    => \$opt_append,
             "unknownok" => \$opt_unknownok);

my $do_prepend = 0;
my $do_append = 0;
my $do_unknownok = 0;
if(defined $opt_prepend)   { $do_prepend = 1; }
if(defined $opt_append)    { $do_append = 1; }
if(defined $opt_unknownok) { $do_unknownok = 1; }

if($do_prepend && $do_append) { die "ERROR choose one of --prepend or --append"; }

if(scalar(@ARGV) != 1) { die $usage; }
my ($fasta_file) = (@ARGV);

if(! -s $fasta_file) { die "ERROR fasta file $fasta_file does not exist or is empty"; }

my $be_verbose = 0;

# first get all sequence names from the fasta file using esl-seqstat
my $seqstat_file = $fasta_file . ".a.seqstat";
my $seqstat_cmd = "esl-seqstat -a $fasta_file > $seqstat_file";
run_command($seqstat_cmd, $be_verbose);

# parse seqstat file to get list of sequence names
my $seqname_list_str = "";
open(SEQSTAT, $seqstat_file) || die "ERROR unable to open $seqstat_file for reading";
my $line;
while($line = <SEQSTAT>) { 
  if($line =~ /^\=\s+(\S+)/) { 
    my $orig_seqname = $1;
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
  }
}
close(SEQSTAT);

# run efetch to get tax strings
my $tax_file = $fasta_file . ".tax";
my $eutils_cmd = "esearch -db nuccore -query \"$seqname_list_str\" | efetch -format gpc | xtract -insd INSDSeq_taxonomy > $tax_file";
run_command($eutils_cmd, $be_verbose);

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

# go through the fasta file and output tax strings in the headers
open(FASTA, $fasta_file) || die "ERROR unable to open $fasta_file for reading";
while($line = <FASTA>) { 
  if($line =~ m/^\>/) { 
    chomp $line;
    if($line =~ m/^\>(\S+)\s*(.*)$/) { 
      my ($orig_seqname, $orig_seqdesc) = ($1, $2);
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
      print ">" . $orig_seqname . " ";
      if($do_append) { print $orig_seqdesc . ";;"; }
      print $cur_tax; 
      if($do_prepend) { print ";;" . $orig_seqdesc; }
      print "\n";
    }
    else { 
      die "ERROR unable to parse fasta header line $line";
    }
  }
  else { # not a header line, regurgitate
    print $line;
  }
}
close(FASTA);

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
