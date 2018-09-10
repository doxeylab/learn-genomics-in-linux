#!/usr/bin/perl
#glimmer_2_gbk.perl - created Thu Aug 29 13:28:50 BST 2013

use strict;
use warnings;

print "FEATURES             Location/Qualifiers\n";

for(<>)
{
    chomp;
    next if /^>/;
    my($orf, $start, $end, $frame, $score) = split(/\s+/, $_);

    if($end >= $start)
    {
	print "     CDS             $start..$end\n";
    }
    else
    {
	print "     CDS             complement($start..$end)\n";
    }
	
    print     "                     /note=\"Predicted by glimmer3, score $score\"\n";
    print     "                     /gene=\"$orf\"\n";
}
