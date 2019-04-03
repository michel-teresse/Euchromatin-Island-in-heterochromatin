# This script, "search.pl" is meant to identify euchromatic island within heterochromatin.
# it is based on the chromatin state of each base. The script detect the typical eucrhromatic
# states that are surrounded by a certain length (determined by the user) of state 9 and 8.
#
# Copyright (c) 2017 Michel TERESE
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

# The 2 following varaibles have to be implemented to the command line
my $MIN_LEN_STATE_9 = 0;
my $MIN_LEN_STATE_8 = 0;

my $FILE_CHROM_STATES; # implemented to the command line
my $verbose;  # 1 if -v is used

my @chr_states_coord;
my @chr_states_string = ('', '', '', '', '');

#================================================================================
# Read the file $FILE_CHROM_STATES in tables @chr_states_coord and @chr_states_string
#================================================================================
sub read_states {
  open my $f, '<', $FILE_CHROM_STATES or die "$FILE_CHROM_STATES $!\n";
  while (<$f>) {
    my ($chr, $start, $end, $state) = split;
    push( @{$chr_states_coord[$chr -1]}, [$start, $end] );
    $chr_states_string[$chr -1] .= $state;
  }
  close $f;

#   for my $chr ( 0 .. $#chr_states_coord ) {
#     my $states = $chr_states_coord[$chr];
#     for my $coord (@{$states}) {
#       my ($start, $end) = @{$coord};
#       print( int($chr) +1, ", $start, $end\n");
#     }
#   }
}


#================================================================================
# Return the coordinates (start, end) of the line $line "zero based" for the chormosome $chr
#================================================================================
sub get_coord {
  my ($chr, $line) = @_;
 
  my $states = $chr_states_coord[$chr -1];
  my $coord = @{$states}[$line];
  my ($start, $end) = @{$coord};
  return ($start, $end);
}


#================================================================================
# Check the state 9 lenght is greater or equal to $MIN_LEN_STATE_9
#================================================================================
sub check_state_9 {
  my ($start, $end) = @_;

  my $len = $end - $start +1;
  if ( $len >= $MIN_LEN_STATE_9 ) {
    print "9 : $start, $end. Len= $len OK\n" if $verbose;
    return 1;
  }
  else {
    print "9 : $start, $end. Len= $len :-(\n" if $verbose;
    return 0;
  }
}


#================================================================================
# Check the state 8 lenght is greater or equal to $MIN_LEN_STATE_8
#================================================================================
sub check_state_8 {
  my ($start, $end) = @_;

  my $len = $end - $start +1;
  if ( $len >= $MIN_LEN_STATE_8 ) {
    print "8 : $start, $end. Len= $len OK\n" if $verbose;
    return 1;
  }
  else {
    print "8 : $start, $end. Len= $len :-(\n" if $verbose;
    return 0;
  }
}


#================================================================================
# Look for the patern in the chromatin states list $states_list for the chromosom $chr
#================================================================================
sub search_patern {
  my ($chr, $state_list) = @_;

  while ( $state_list =~ /(98)((?!9).+?)(89)/gc ) {
    my $borne_inf_line_start = $-[1];
    my $borne_inf_line_stop = $+[1] -1;

    my $borne_ilot_line_start = $-[2];
    my $borne_ilot_line_stop = $+[2] -1;

    my $borne_sup_line_start = $-[3];
    my $borne_sup_line_stop = $+[3] -1;

    # Place the start of the next search in $state_list one base before
    # i.e. on the "9" on the right. The next island left border can be the right border 
    # of the previous island.It is mandatory to include the last right border found 
    # in the next search.
    pos($state_list) = $borne_sup_line_stop;
    
#     print "borne inf line start= $borne_inf_line_start end= $borne_inf_line_stop\n"; 
#     print "borne ilot line start= $borne_ilot_line_start end= $borne_ilot_line_stop\n"; 
#     print "borne sup line start= $borne_sup_line_start end= $borne_sup_line_stop\n"; 

    my ($start, $end, $dummy);
    ($start, $dummy) = get_coord( $chr, $borne_ilot_line_start );
    ($dummy, $end)   = get_coord( $chr, $borne_ilot_line_stop );
    my $ilot_len = $end - $start +1;

    if ( check_state_9( get_coord( $chr, $borne_inf_line_start ) ) &&
         check_state_8( get_coord( $chr, $borne_inf_line_stop ) ) &&
         check_state_8( get_coord( $chr, $borne_sup_line_start ) ) &&
         check_state_9( get_coord( $chr, $borne_sup_line_stop ) )
       ) {

      print "$chr\t$start\t$end\n";
      print "len= $ilot_len\n" if $verbose;
    }
    else {
      print "rejected: len= $ilot_len\n" if $verbose;
    }

    print "------\n\n" if $verbose;
  }
}


#================================================================================
# Check parameters
#================================================================================
sub check_param {
	GetOptions ( 'verbose|v'  => \$verbose)
		or die("Syntax: $0 [-v]\n");
		
	if (@ARGV != 3 || ! -f $ARGV[0]) {
	    die "syntax: $0 [-v] chrom_state_header_less MIN_LEN_STATE_8 MIN_LEN_STATE_9\n";
	}

	($FILE_CHROM_STATES, $MIN_LEN_STATE_8, $MIN_LEN_STATE_9) = @ARGV;
}


#================================================================================
# main
#================================================================================
check_param();
read_states();

for my $chr (1..5) {
  search_patern( $chr, $chr_states_string[$chr -1] );
}

