#!perl
use strict;
use warnings;
use local::lib;
#use Data::Dumper;

use JSON;
use Bio::GFF3::LowLevel qw / gff3_format_feature  /;

my $ASSEMBLY = 'pfal3D7';
my $OUT = $ASSEMBLY.'_gene%3Aannotation2.gff';
open OUT, ">$OUT" or die "couldn't open $OUT for writing: $!";

my $contiglist_blob;
{
    local $/ = undef;
    my $file = "pf_seqlist.json";
    open FF, "<$file" or die "couldn't open $file: $!";
    $contiglist_blob = <FF>;
    close FF;
}

my $contig_json = JSON->new->decode($contiglist_blob);

for my $contig_info (@{$contig_json}) {
    my $contig_name = $$contig_info{'name'};
   
    my $fetch_url = "https://plasmodb.org/a/service/jbrowse/features/$contig_name?feature=gene%3Aannotation2&start=0&end=2038340";

    system("curl", '-o', $contig_name."_gene%3Aannotation2.json",$fetch_url);

    my $blob;
    {
    local $/ = undef;
    my $file = $contig_name."_gene%3Aannotation2.json";
    open FF, "<$file" or die "couldn't open $file: $!";
    $blob = <FF>;
    close FF;
    unlink $file;
    }

    my $json = JSON->new->decode($blob);

    for my $feature (@{$$json{'features'}}) {
        &parse_line(undef, undef, $feature);
    }
}
close OUT;
exit 0;

sub parse_line {
    my $parent_id = shift;
    my $contig    = shift;
    my $f         = shift;

    my %gff3_feature;
    if ($contig) {
        $gff3_feature{'seq_id'} = $contig;
    }
    else {
        $gff3_feature{'seq_id'} = $$f{'contig'};
    }
    $gff3_feature{'source'} = $$f{'source'};
    $gff3_feature{'type'}   = $$f{'type'};
    $gff3_feature{'start'}  = $$f{'startm'}; # the base coord start
    $gff3_feature{'end'}    = $$f{'end'};
    $gff3_feature{'phase'}  = $$f{'phase'};
    if ($$f{'strand'} == 1) {
        $gff3_feature{'strand'} = '+';
    }
    elsif ($$f{'strand'} == -1) {
        $gff3_feature{'strand'} = '-';
    }
    else {
        $gff3_feature{'strand'} = $$f{'strand'};
    }
    $gff3_feature{'score'}  = $$f{'score'};

    my %attributes;
    $attributes{'Parent'} = [ $parent_id ] if $parent_id;
    $attributes{'ID'}     = [ $$f{'feature_id'} ]; #attributes have to be lists
    $attributes{'Name'}   = [ $$f{'name'} ];
    $attributes{'Note'}   = [ $$f{'note'} ];

    for my $key (keys %{$f}) {
         next if $key =~ /subfeatures|feature_id|note|name|contig|source|type|start|end|phase|strand|score/;
         $attributes{$key} = [ $$f{$key} ];
    }

    $gff3_feature{'attributes'} = \%attributes;

    print OUT gff3_format_feature(\%gff3_feature);

    # print JSON->new->pretty->encode(\@{$f{'subfeatures'}});

    for my $sub (@{$$f{'subfeatures'}}) {
        &parse_line(@{$attributes{'ID'}}[0], $gff3_feature{'seq_id'}, $sub);
    }
}
