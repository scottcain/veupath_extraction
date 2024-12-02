#!/usr/bin/perl
use strict;
use warnings;

#use Bio::GFF3::LowLevel qw/ gff3_parse_feature /
use URI::Escape;

my $FILE = $ARGV[0];
open FILE, "<$FILE" or die "couldnt open $FILE for reading: $!";

#get jbrowse track config info from the GFF header
my %confhash;
while (<FILE>) {
    chomp;

    if (/^#/) {
        if ( /\s([^\s]+)\s=\s(.*)/) {
            $confhash{$1} = $2;
	}
    } 
    last;
}

#create a (hopefully) cleaner name
my $sorted;
if (defined $confhash{'label'}) {
    $sorted = $confhash{'label'};
}
elsif (defined $confhash{'key'}) {
    $sorted = $confhash{'key'};
}
else {
    $sorted = $FILE;
}
$sorted =~ s/ /_/g;
$sorted = escape($sorted);
$sorted =. 'sorted.gff3';

open CONF, ">$sorted.json" or die "couldnt open $sorted.json for writing: $!";

#deal with track specific json
if (defined $confhash{'label'}) { }
if (defined $confhash{'key'}) { }
if (defined $confhash{'category'}) { }
if (defined $confhash{'metadata.dataset'}) { }
if (defined $confhash{'metadata.description'}) { }
if (defined $confhash{'metadata.subcategory'}) { }
if (defined $confhash{'subParts'}) { }
if (defined $confhash{'metadata.attribution'}) { }

#print out jbrowse track json
print CONF <<END;
"tracks" : [
    {
        
    }
]

END


#bgzip/tabix the GFF file
system("jbrowse sort-gff $FILE > $sorted") == 0 or die $!;
system("bgzip $sorted");
system("tabix $sorted.gz");



