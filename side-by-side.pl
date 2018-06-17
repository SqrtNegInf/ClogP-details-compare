#!/usr/local/bin/perl
## compare ClogP details tables, side-by-side
#  David H.  2018-06-14
#  bug: another fragment gets doubled due to this calc. fragment: *[o+]* [aa]

use warnings;
use strict 'vars';

use CGI qw/:standard :cgi/;
use CGI::Carp qw(fatalsToBrowser);
use List::AllUtils qw(any extract_by zip_by);

my $DEBUG = 0;

# external tools/files we'll need
$ENV{BB_ROOT} = '/biobyte/ClogP/dev';
my $clogp   = "$ENV{BB_ROOT}/bin/clogp";
my $usmiles = "$ENV{BB_ROOT}/bin/usmiles";
my $dfilter = '/biobyte/bin/details-filter.pl';
my $master  = '/bb/db/MlogP/master.smi';
my $pcre     = '/usr/local/bin/pcregrep';

# display form
print 
    header,
    start_html('ClogP details side-by-side'),
	start_form(-method=>'GET'),
    h3('ClogP details side-by-side'),
	"<b>Compound 1: </b>",textfield(-name=>'c1', -size=>100), p,
	"<b>Compound 2: </b>",textfield(-name=>'c2', -size=>100), p,
	submit(-name=>'Compare'), ' ',
	defaults('Reset'),
	end_form, p, hr, p;

# do not proceed unless have two compound candidates, and both resolve to valid SMILES
my $c1  = CGI::param('c1');
my $c2  = CGI::param('c2');
bail() unless $c1 && $c2;
my($smi1,$nam1) = getsmi($c1);
my($smi2,$nam2) = getsmi($c2);
bail() unless $smi1 && $smi2;

# get details tables, discard 'SMILES:' line
my @table1 = `$clogp -2 '$smi1' | $dfilter`;
my @table2 = `$clogp -2 '$smi2' | $dfilter`;
shift @table1;
shift @table2;

# display results
print '<pre><br>';
side_by_side($smi1,$nam1,\@table1,$smi2,$nam2,\@table2);
print '</pre>';
print end_html;

#####
# arrange details tables side-by-side, pairing lines as feasible, and annotating +/- changes
sub side_by_side {
my($smi1,$nam1,$table1,$smi2,$nam2,$table2) = @_;

# pre-flight special sorts for 'Fragment' and 'Proximity' corrections
my @table1 = sort_sections(@$table1);
my @table2 = sort_sections(@$table2);

# collate raw details tables by correction type
my @raw;
for my $s (extract_sections(@table1,@table2)) {

    # need at least one item to work with
#print "s1: $s<br>";
    $s =~ s/([^.])\*/$1./g;
    $s =~ s/[()]/./g;
#print "s2: $s<br>";
    my @s1 = grep { /^$s/ } @table1;
    my @s2 = grep { /^$s/ } @table2;
    next unless @s1 or @s2;

    # pad shorter table
    my $l1 = +@s1;
    my $l2 = +@s2;
    my $del= abs ($l1-$l2);
    my $pad = '          |      |                                        |          |      ';
    if ($l1>$l2) { for (1..$del) { push @s2, $pad } }
    else         { for (1..$del) { push @s1, $pad } }

    # abbreviate and build raw side-by-side; section labels only at left margin
    while (my $a = abbrev(shift @s1), my $b = abbrev(shift @s2)) { 
        my $label = $a =~ /^ / ? substr($b,0,10) : substr($a,0,10);
        push @raw, [split //, sprintf qq{%10s|%65s ## %s}, $label, substr($a,11), substr($b,11)];
    }
}   
    
# narrow tables by removing blank columns (invert, filter, re-invert)
my(@narrow,@cooked);
for my $r (zip_by { [ @_ ] } @raw) {
    next unless any { $_ ne ' ' } @$r;
    push @narrow, [@$r];
}
for my $r (zip_by { [ @_ ] } @narrow) {
    push @cooked, join '', @$r;
}

# sanity check
if ($DEBUG) {
my($lcnt,$rcnt);
for my $line (@cooked) {
    my($l,$r) = split /##/, $line;
    $lcnt++ if ' ' ne substr($l, -1, 1);
    $rcnt++ if ' ' ne substr($r, -1, 1);
}
printf qq{\n<p><font color="red">ERROR: $smi1 $nam1 vs $smi2 $nam2</font>\n} 
    if (1+$lcnt) ne scalar(@table1) or (1+$rcnt) ne scalar(@table2);
}

# format & display: header
my($l,$r) = split /##/, $cooked[0];
printf "$nam1 %s $nam2<br>", ' ' x (1+length($l)-length($nam1));
printf "%s   %s<p>", ('_' x length($l)), ('_' x length($r)); 

# format & display: details table
for my $l (@cooked) {
    my($dn,$del,$dir,$sign,$col);
    my($a,$an,$b,$bn) = $l =~ m/.*?\|(.*?\|\s*([^|]+))##(.*?\|\s*([^|]+))$/;

    # annotate direction/degree of correction difference
    $an = 0 if $an eq ' ' || $an eq '#';
    $bn = 0 if $bn eq ' ' || $bn eq '#';
    $dn  = abs($bn - $an);
    $del = ' ';
    $dir = ' ';
    if ($dn > 0) {
        $sign = $bn > $an ? '+' : '-';
        $del  = '.';
        $del  = $sign if $dn > 0.11;
        $del .= $sign if $dn > 0.50;
        $del .= $sign x int(-.25+$dn);
        $dir = $sign eq '-' ? '>' : '<';
    }

    # de-emphasize identical/similar lines
    $col = nws($a,$b) ? 'lightgray' : 'black';
    $col = 'gray' if $dn > 0 && $dn <= 0.10;

    # modify line with annotations, legibility tweaks
    $l =~ s!##!<font color="blue"><b> $dir </b></font>!;
    $l =~ s!(\|\d+) !<b>$1</b> !g   if $dn > 0;
    $l =~ s!, (\d+) !, <b>$1</b> !g if $dn > 0;

    # at last, output!
    print qq{<font color="$col">$l</font> $del<br>};
}

# 2D images below table for reference, if available
print '<p><hr><p>';
print '        '; printf qq{<img src="%s">}, get_gif_ref($smi1);
print '        '; printf qq{<img src="%s">}, get_gif_ref($smi2);

}

#####
# abbreviate descriptions, but maintain field width
sub abbrev {
return unless defined $_[0];
my($l) = @_;
chomp $l;
my @f = split /\|/, $l;

# 'class' field
$f[0] =~ s/^Fragment/Frag/;
$f[0] =~ s/^ExFragment/ExFrag/;
$f[0] =~ s/^Proximity/Proxim/;
$f[0] =~ s/^Electronic/Elect/;
$f[0] =~ s/^Benzylbond/Benzyl/;
$f[0] =~ s/^Fragbranch/Fragbr/;
$f[0] =~ s/^HaloInFrag/HaloFr/;
$f[0] =~ s/^Sigma-I/SigmaI/;
$f[0] =~ s/^Glycoside/Glyco/;
$f[0] =~ s/^OrthoEther/OrthoE/;
$f[0] =~ s/^Zwitterion/Zwitt/;
$f[0] =~ s/^Dpl shield/Dpl sh/;
$f[0] .= ' ' x (10-length($f[0]));

# 'description' field
$f[2] =~ s/extended hetero-aromatic/ext. hetero-ar/;
$f[2] =~ s/benzyl bond to/benzyl bond,/;
$f[2] =~ s/simple aromatics?/simple arom/;
$f[2] =~ s/isolating carbons/iso-C/;
$f[2] =~ s/iso-C's/iso-C/;
$f[2] =~ s/potential interactions?;/potential,/;
$f[2] =~ s/multihetero-5-ring/multihet-5-ring/;
$f[2] =~ s/measured fragment value/measured/;
$f[2] =~ s/from measured value/from measured/;
$f[2] =~ s/phenyl-fragment/phenyl-frag/;
$f[2] =~ s/non-halogen, polar group branch(es)?/non-halo, polar gr. branch/;
#$f[2] =~ s/polar group branch/polar gr. branch/;
$f[2] =~ s/ and 0 alicyclic//;
$f[2] =~ s/, 0 triple//;
$f[2] =~ s/0 double, //;
$f[2] =~ s/and 0 cluster //;
$f[2] =~ s/ and 0 cluster branch(es)?//;
$f[2] =~ s/0 chain and//;
$f[2] =~ s/chain and/chain,/;
$f[2] =~ s/Ether in a/Ether,/;
$f[2] =~ s/not applicable in modeling//;
$f[2] =~ s/in nature//;
$f[2] =~ s/\(net\)//;

$f[2] =~ s/\s+/ /g;
$f[2] .= ' ' x (40-length($f[2]));

# 'comment' field
$f[3] =~ s/WithinRing/Within/;
$f[3] =~ s/FusedRings/Fused/;
$f[3] =~ s/JoindeRing/Joined/;
$f[3] =~ s/Calculated/Calc/;
$f[3] =~ s/Pyridinium/Pyridin/;
$f[3] =~ s/Ext.Fusion/ExFusion/;
$f[3] =~ s/Keto.enol/Ket-Enol/;
$f[3] .= ' ' x (10-length($f[3]));

return join '|', @f;
}

#####
# sort sections to align them for upcoming comparison
sub sort_sections {
my(@t) = @_;

# arrange 'Fragment' entries alphabetically, based on 'description' field
my @fr = sort { lc(substr($a,19,38)) cmp lc(substr($b,19,38)) } extract_by { /^Fragment/ } @t;

# arrange 'Proximity' entries by custom sort
my @pr = sort { prcmp($a) cmp prcmp($b) } extract_by { /^Proximity/ } @t;

# combine sorted section with any remaining unsorted
return @fr, @pr, @t;
}

#####
# re-arrange fields to compare by: 1) proximity type 2) correction value 3) description
sub prcmp { 
my($l) = @_; 
chomp $l;
my @f = split /\|/, $l;
$l = join '|', $f[0], $f[1], $f[3], $f[4], $f[2];
$l =~ s/\(.*?\)/ /; 
$l =~ s/#\d+/ /g; 
$l =~ s/\s+/ /g; 
return $l;
} 

#####
# test string equality, modulo white-space, and non-significant numerics (fragment #, ring #, integer counts, etc)
sub nws { 
my($a,$b) = @_;

# H-bond & Ortho
$a =~ s/Ring ?\d+//;
$b =~ s/Ring ?\d+//;
$a =~ s/\d+ potential, \S+//;

# Electronic
$a =~ s/\d+ potential, \S+//;
$b =~ s/\d+ potential, \S+//;

# Fusion
$a =~ s/\d+ extended//;
$b =~ s/\d+ extended//;

#Fragments
$a =~ s/\(#\d+\)//;
$b =~ s/\(#\d+\)//;
$a =~ s/ #(\s\d|\d\d) //; 
$b =~ s/ #(\s\d|\d\d) //; 
$a =~ s/Fragments #\d+ & #\d+//; 
$b =~ s/Fragments #\d+ & #\d+//; 
$a =~ s!1 pair .#\d+/\d+.!!; 
$b =~ s!1 pair .#\d+/\d+.!!; 
$a =~ s!pairs .#\d+/\d+,\d+/\d+.!!;
$b =~ s!pairs .#\d+/\d+,\d+/\d+.!!;

# Ortho
$a =~ s!ortho .#[b0-9]+/[b0-9]+.!!;
$b =~ s!ortho .#[b0-9]+/[b0-9]+.!!;

# strip white-space and test
$a =~ s/\s//g; 
$b =~ s/\s//g; 
return $a eq $b; 
}

#####
# load table section definition patterns (NB: ordering is significant, controls output order)
sub extract_sections {
my(@t) = @_;

return 
idio('Fragment  . ...  .','\].*','\[',@t),
"Carbon.*aromatic",
"Carbon.*aliphatic",
"Fusion.*extended aromatic.*Joined",
"Fusion.*extended aromatic.*Fused",
"Fusion.*extended hetero.*Joined",
"Fusion.*extended hetero.*Fused",
"ExFragment.Hydrog",
"ExFragment.Branch",
"ExFragment.Bonds",
"ExFragment.Mbonds",
"ExFragment.CMbond",
"Benzylbond.Hetero",
"Benzylbond.Simple",
"Fragbranch",
"OrthoEther. ",
"OrthoEther.Vinyl",
"OrthoEther.Aromat",
"Dpl shield. ",
"Dpl shield.XCX",
"Dpl shield.XCCX",
"Proximity .XCCY",
"Proximity .XCCCY",
"Proximity .PCCY",
"Proximity .Y1-X1",
"Proximity .Y1-X2",
"Proximity .Y1-X3",
"Proximity .Y2-X1",
"Proximity .Y2-X2",
"Proximity .Y2-X3",
"Proximity .Y3-X1",
"Proximity .Y3-X2",
"Proximity .Y3-X3",
"Proximity .YCY",
"Proximity .YCCY.*1 pair",
"Proximity .YCCY.*2 pair",
"Proximity .YCCY.*3 pair",
"Proximity .YCCCY",
"Proximity .YCC=C",
"Electronic.*Within",
"Electronic.*Joined",
"Electronic.*Fused",
"Electronic.*Ext-Fusion",
"Glycoside . ",
"Glycoside .Global",
"Glycoside .Ring",
"H-bond",
"Ortho     .Benz",
"Ortho     .Ring",
"Zwitterion",
"Sigma-I",
"Steroid",
"HaloInFrag",
"Screen    .CyProp",
"Screen    .K-BiPh",
idio('Screen    .SMARTS.','','',@t),
"Screen    . ",
"RESULT";

}

# special handling for section with non-template 'Description' fields (fragments, SMARTS)
sub idio {
my($prefix,$pat1,$pat2,@table) = @_;
my(%S,@sk);
for my $s (sort { substr($a,19,38) cmp substr($b,19,38) } extract_by { /^$prefix/ } @table) {
    my $s = substr($s,19,38);
    $s =~ s/[+*()]/./g;
    $s =~ s!$pat1!.!g if $pat1;
    $s =~ s!$pat2!.!g if $pat2;
    $s =~ s/  .*//;
    $S{"$prefix $s"}++;
}
for my $k (sort keys %S) { push @sk, $k }
return @sk;
}

#####
# get canonical SMILES from input (SMILES/name)
sub getsmi {
my($id) = @_;
my($nam,$match);

my $smi = `$usmiles '$id'`; chomp $smi;
if ($smi) {
    $match = `grep '^$id ' $master`;
} else {
    $id =~ tr/a-z/A-Z/;
    $id =~ s/[^A-Z0-9]//g;
    $match = `grep ' $id\$' $master`;
}
chomp $match;
($smi,$nam) = split / /, $match if $match;

return $smi, $nam;
}

#####
# return path to image for SMILES
sub get_gif_ref {
my($smi) = @_;

my $depict_dir = 'depictions/thor_depictions';
my $serial = '/biobyte/htdocs/bb-server/databases/thor.serialno';
my $qsmi = quotemeta $smi;

my $result = `$pcre '^$qsmi ' $serial`;
my $dd;
if ($result) {
    chomp $result;
     my($d) = $result =~ / (\d+)/;
     my($d3,$d2,$d1) = $d =~ /^(((.).).)/;
    if    ( $d <  1000 ) { $dd = '00'; }
    elsif ( $d < 10000 ) { $dd = sprintf "%02d", $d1; }
    else                 { $dd = $d < 100000 ? $d2 : $d3; }
    $result = "/$depict_dir/$dd/$d.GIF";

# try harder
} else {
    $depict_dir = 'depictions/qsar_depictions';
    $serial = '/biobyte/htdocs/bb-server/databases/bio.serialno';
    $result = `$pcre '^$qsmi ' $serial`;
    if ($result) {
    chomp $result;
    my($d) = $result =~ / (\d+)/;
    my($d3,$d2,$d1) = $d =~ /^(((.).).)/;
    if    ( $d <  1000 ) { $dd = '00'; }
    elsif ( $d < 10000 ) { $dd = sprintf "%02d", $d1; }
    else                 { $dd = $d < 100000 ? $d2 : $d3; }
    $result = "/$depict_dir/$dd/$d.GIF";
    }
}

if (! $result) {
    my $fullpath = '/bb/db/master/ZLN-dict/delta/babel-bb';
    my $webpath = '/depictions/babel-bb';
    $serial = $fullpath . '/babel.serial';
    my $query = qq{$pcre '^$qsmi ' $serial};
    $result = `$pcre '^$qsmi ' $serial`;
    chomp $result;
    my($file) = $result =~ / (.*)/;
    $result = "$webpath/$file.GIF";
}

$result = "depictions/black_370_277.GIF" if $result =~ m#/\.GIF#; # placehoder if no image found
return $result
}

#####
# warn user about insufficient input
sub bail {
    print '<h3><b>Need two valid SMILES or compound names!</b>' if defined $ENV{QUERY_STRING} && $ENV{QUERY_STRING} =~ /compare/i;
    exit;
}
