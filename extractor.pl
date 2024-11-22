#!perl
use strict;
use warnings;
use local::lib;

use JSON;
use Bio::GFF3::LowLevel qw / gff3_format_feature  /;


my $blob;
{
    local $/ = undef;
    my $file = "chr_1_pf_gene_annotation2.json";
    open FF, "<$file" or die "couldn't open $file: $!";
    $blob = <FF>;
    close FF;
}

my $json = JSON->new->decode($blob);

for my $feature (@{$$json{'features'}}) {
    &parse_line(undef, undef, $feature);
}

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

    print gff3_format_feature(\%gff3_feature);

    # print JSON->new->pretty->encode(\@{$f{'subfeatures'}});

    for my $sub (@{$$f{'subfeatures'}}) {
        &parse_line(@{$attributes{'ID'}}[0], $gff3_feature{'seq_id'}, $sub);
    }
}
