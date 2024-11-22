#!perl
use strict;
use warnings;
use local::lib;
#use Data::Dumper;

use JSON;
use Bio::GFF3::LowLevel qw / gff3_format_feature  /;

my $ASSEMBLY = 'pfal3D7';

#this is the standard tracks.conf; ignores for the moment other track configs
#"/a/service/jbrowse/dnaseq/pfal3D7" BAMs/BWs
#"/a/service/jbrowse/rnaseq/pfal3D7" BAMs/BWs
#"/a/service/jbrowse/chipseq/pfal3D7" BWs
#"/a/service/jbrowse/rnaseqJunctions/pfal3D7"
#"/a/service/jbrowse/organismSpecific/pfal3D7"
my $TRACK_CONF = "tracks.conf";
my %vuepath_track_info;
my %lh;
my $track_key;
open TRACKCONF, "<$TRACK_CONF" or die "couldn't open $TRACK_CONF for reading: $!";
while (<TRACKCONF>) {
    next if /^#/;
    chomp;

    if (/^\[(.*)\]/) {
        $track_key = $1;
    }
    elsif (/^([^=]+)=(.*)/) {
        $vuepath_track_info{$track_key}{$1} = $2;
    }

}
close TRACKCONF;


my $contiglist_blob;
{
    local $/ = undef;
    my $file = "pf_seqlist.json";
    open FF, "<$file" or die "couldn't open $file: $!";
    $contiglist_blob = <FF>;
    close FF;
}
my $contig_json = JSON->new->decode($contiglist_blob);

for my $tr_key (keys %vuepath_track_info) {
    next if $tr_key eq 'tracks.refseq';
    next if exists $vuepath_track_info{$tr_key}{'query.edName'};
    next if exists $vuepath_track_info{$tr_key}{'query.edname'};
    next if $vuepath_track_info{$tr_key}{'query.feature'} eq 'ReferenceSequence';

    my $track_name = $vuepath_track_info{$tr_key}{'key'};
    my $query_feature = $vuepath_track_info{$tr_key}{'query.feature'};
    next unless $query_feature;
    $query_feature =~ s/\./%3A/g;

    my $OUT = $ASSEMBLY.'_'.$track_name.'.gff';
    open OUT, ">$OUT" or die "couldn't open $OUT for writing: $!"; 

    print OUT "##gff-version 3\n";
    print OUT "# VEuPathDB track config info:\n";
    for my $config_key (keys %{$vuepath_track_info{$tr_key}}) {
        print OUT "# $config_key = $vuepath_track_info{$tr_key}{$config_key}\n";
    }

    for my $contig_info (@{$contig_json}) {
        my $contig_name = $$contig_info{'name'};
	my $contig_end  = $$contig_info{'end'};
   
        my $fetch_url = "https://plasmodb.org/a/service/jbrowse/features/$contig_name?feature=$query_feature&start=0&end=$contig_end";
        warn $fetch_url;

        system("curl", '-o', $contig_name."_$query_feature",$fetch_url);

        my $blob;
        {
        local $/ = undef;
        my $file = $contig_name."_$query_feature";
        open FF, "<$file" or die "couldn't open $file: $!";
        $blob = <FF>;
        close FF;
        unlink $file;
        }

        my $json = JSON->new->decode($blob);

        for my $feature (@{$$json{'features'}}) {
            &parse_line(undef, $contig_name, $feature);
        }
    }
    close OUT;
    sleep 10;
    die;
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

    print OUT gff3_format_feature(\%gff3_feature);

    # print JSON->new->pretty->encode(\@{$f{'subfeatures'}});

    for my $sub (@{$$f{'subfeatures'}}) {
        &parse_line(@{$attributes{'ID'}}[0], $gff3_feature{'seq_id'}, $sub);
    }
}
