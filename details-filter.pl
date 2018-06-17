#!/usr/local/bin/perl
## simplify and standardize ClogP details table by removing ASCII 'boxes', rounding to two digits after decimal
#  input can be details for a single compound or multiple compounds
#  handles VAX (BDRIVE/verbose), UNIX (clogp/verbose & extended), & Windows (Bio-Loom/verbose) input
#  after processing, output from all platforms is directly comparable via 'diff'

use strict 'vars';
use vars qw(@table $name $smi);

use Math::Round qw(/nearest/);

while (<>) {
    exit if eof();

    # skip electronic detail sub-table (VAX only)
	if (/React:/) {
		do {
			$_ = <>;
            exit if eof();
		} until /^>>>/;
    }

    exit if eof();

	# skip leading ASCII 'boxing' characters
	next unless /^ ?\|/ || /^>>>/ || /^C:/;

	# tidy up: left/right margin box chars, extraneous VAX escape sequenc𝑒
	s/^\s*\|//; s/\|$//; 
 	s/\000/ /;

    # standardize variant numerical notations, and round to 2 digits after decimal
    unless (/(Name|SMILES|C):/i) {
	    s/  0\./   ./; s/ -0\./  -./;
        my($n) = /(......)$/; 
        my $r = nearest(.01, $n);
        $r = sprintf "%5.2f", $r;
        $r =~ s/([ \-])0/ $1/;
        s/......$/$r/;
    }

	# extract SMILES and name (VAX input)
	if ( /^>>>/) {
		<>; $_ = <>;
		if (/NAME:/) {
			s/^\s*\|//; s/\|$//; 
			s/\s+$//; 
			($name) = /NAME: (.*)/; 
		} else {
			$name = 'noname';
		}
		if (!/SMILES:/) { do { $_= <> } until /SMILES:/ }
 		($smi) = /SMILES: (\S+)/;
		do {
			$_ = <>;
 			if (/SMILES:/) { /SMILES: (\S+)/; $smi .=  $1 }
		} until /Contribution/;

	# extract name, and look up SMILES (Bio-Loom input)
	} elsif ( /^C:/ ) {
		$name = substr($_,21);
        chomp $name;
        $name =~ s/ # .*//;
		do { $_ = <> } until /Contribution/;
        $smi = `grep ' $name\$' /bb/db/MlogP/master.smi | head -1`;
        $smi =~ s/ .*//;
        chomp $smi;

	# extract SMILES and name (UNIX input)
	} elsif ( /Name:/i || /SMILES:/ ) {
		s/  +/ /g;

        # build up name (possibly multi-line)
		if (/Name:/i) {
 			($name) = /NAME: (.*)|$/i; 
			$name =~ s/ $//;
 			do {
				$_ = <>;
				s/  +/ /g;

                # munge text to extract just name, working around extraneous material
 				if (/NAME:/i) { 
                    /NAME: (.*)|$/i; 
                    $name .= $1;  
                    $name =~ s/\s*\|$//;
                    $name =~ s/(\w)(-?\d+\.\d+)/$1 $2/;
                    $name =~ s/#2([\-\d])/#2 $1/;
                    $name =~ s/([)*])([\d.]+)$/$1 $2/;
                    $name =~ s/([^ ])\*/$1 */;
                    $name =~ s/\*([^ ])/* $1/;
                    $name =~ s/\* \*/**/;
                }
			} until /SMILES:/;
		}
 		($smi) = /SMILES: (\S+)/;

        # build up SMILES (possibly multi-line)
		do { $_ = <>; if (/SMILES:/) { /SMILES: (\S+)/; $smi .=  $1 } } until /Contribution/;

	# we have reached the end of the table: print everything, reset '@table'
	} elsif (/RESULT/) {
        
        # blank out any fragment DB and ClogP version info
		s/DB=\d\d/     /;
		s/CLOGP=\d.[\d ]{2}/          /;
		s/CLOGP/     /;
		s/ \d\.\d\d \|/      |/;
		s#\| \d/\d\d[ +]\|#|      |#;

		print "Name: $name\n" if $name;
		print "SMILES: $smi\n";
		foreach my $i (@table) { print "$i\n" }
		undef @table;
		print "$_\n";

	# all other lines added to table 
	} else {
		chop; push(@table,$_);

	} 
}
