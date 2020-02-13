#!/usr/bin/env perl
# EPN, Tue Jan  1 06:32:21 2019
use warnings;
use strict;
use Getopt::Long;

my $usage;
$usage  = "perl fasta-add-tax-as-desc-efetch.pl [OPTIONS] <fasta file (cannot be read from stdin)>\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t--name <s>:  add taxonomy string tokens <d1> to <d2> to the sequence name (<s>=<d1>-<d2>, or 'all' for all) [default: replace current desc with tax string]\n";
$usage .= "\t\t--prepend:   add taxonomy string at beginning of existing seq description [default: replace current desc with tax string]\n";
$usage .= "\t\t--append:    add taxonomy string at end       of existing seq description [default: replace current desc with tax string]\n";
$usage .= "\t\t--unknownok: output UNKNOWN-TAXONOMY (instead of failing) if taxonomy of a seq is unknown because it is not in nuccore (wgs, maybe?)\n";

my $opt_name      = undef;
my $opt_prepend   = undef;
my $opt_append    = undef;
my $opt_unknownok = undef;
&GetOptions( "name=s"    => \$opt_name,
             "prepend"   => \$opt_prepend, 
             "append"    => \$opt_append,
             "unknownok" => \$opt_unknownok);

my $do_name    = 0;
my $do_prepend = 0;
my $do_append = 0;
my $do_unknownok = 0;
if(defined $opt_name)      { $do_name      = 1; }
if(defined $opt_prepend)   { $do_prepend   = 1; }
if(defined $opt_append)    { $do_append    = 1; }
if(defined $opt_unknownok) { $do_unknownok = 1; }

if($do_name    && $do_append)  { die "ERROR choose one of --name or --append"; }
if($do_name    && $do_prepend) { die "ERROR choose one of --name or --prepend"; }
if($do_prepend && $do_append)  { die "ERROR choose one of --prepend or --append"; }

my $name_d1 = undef;
my $name_d2 = undef;
my $name_all = 0;
if(defined $opt_name) { 
  if($opt_name eq "all") { 
    $name_all = 1;
  }
  else { 
    if($opt_name =~ /^(\d+)\-(\d+)$/) {
      ($name_d1, $name_d2) = ($1, $2);
      if($name_d1 >= $name_d2) { 
        die "ERROR with --name <s>, <s> must equal \"<d1>-<d2>\" with <d1> <= <d2>";
      }
    }
    else { 
      die "ERROR with --name <s>, <s> must equal \"all\" or \"<d1>-<d2>\" with <d1> <= <d2>";
    }
  }
}

if(scalar(@ARGV) != 1) { die $usage; }
my ($fasta_file) = (@ARGV);

if(! -s $fasta_file) { die "ERROR fasta file $fasta_file does not exist or is empty"; }

my $be_verbose = 0;

# first get all sequence names from the fasta file using esl-seqstat
my $seqstat_file = $fasta_file . ".a.seqstat";
my $seqstat_cmd = "esl-seqstat -a $fasta_file > $seqstat_file";
run_command($seqstat_cmd, $be_verbose);

# read list file
my @seqname_list_str_A = ();
my $seqname_list_str = "";
my $list_size = 0;
my $block_size = 100;
my $max_width = 0;
open(SEQSTAT, $seqstat_file) || die "ERROR unable to open $seqstat_file for reading";
my $line;
while($line = <SEQSTAT>) { 
  chomp $line;
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
    $list_size++;
    if($list_size == $block_size) { 
      push(@seqname_list_str_A, $seqname_list_str);
      $seqname_list_str = "";
      $list_size = 0;
    }
  }
}
close(SEQSTAT);
if($list_size > 0) { 
    push(@seqname_list_str_A, $seqname_list_str);
}  

# run efetch to get tax strings
my $tax_file = $fasta_file . ".tax";

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
      if($do_name) { 
        # create new name from tax string 
        my $new_name = $orig_seqname . ";"; 
        my $tax_str = "";
        if($name_all) { 
          $tax_str = $cur_tax;
        }
        else { 
          my @tax_A = split("; ", $cur_tax); 
          my $tax_start = $name_d1;
          my $tax_stop  = $name_d2;
          if($tax_stop > (scalar(@tax_A))) { $tax_stop = scalar(@tax_A); }
          for(my $t = ($tax_start-1); $t <= ($tax_stop-1); $t++) { 
            $new_name .= $tax_A[$t] . ";";
          }
        }
        print ">" . $new_name . "\n";
      }
      else { 
        # modify desc with tax string
        print ">" . $orig_seqname . " ";
        if($do_append) { print $orig_seqdesc . ";;"; }
        print $cur_tax; 
        if($do_prepend) { print ";;" . $orig_seqdesc; }
        print "\n";
      }
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
