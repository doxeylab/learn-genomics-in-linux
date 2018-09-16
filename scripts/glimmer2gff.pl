#!/usr/bin/perl -w

my $usage = "Usage: $0 <sequence name> [<glimmer coord file(s)>]\n";

die $usage unless @ARGV;
my $seqname = shift;

while (<>) {
    chomp;
    s/^\s*//;
    my ($id, $start, $end, $rest) = split /\s+/, $_, 4;
    my $strand;
    ($start, $end, $strand) = $end > $start ? ($start, $end, "+") : ($end, $start, "-");
    print join ("\t", $seqname, "glimmer", "gene", $start, $end, ".", $strand, ".", $rest), "\n";
}
