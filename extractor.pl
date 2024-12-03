#!perl
use strict;
use warnings;
use local::lib;
use Data::Dumper;

use JSON;
use Bio::GFF3::LowLevel qw / gff3_format_feature  /;
use URI::Escape;

my $ASSEMBLY = $ARGV[0];
my $DIVISION = $ARGV[1];
my $SEQINFO              = "https://$DIVISION/a/service/jbrowse/seq/$ASSEMBLY";
my $TRACKINFO            = "https://$DIVISION/a/jbrowse/tracks/$ASSEMBLY/tracks.conf";
my $ASSEMBLYSPECIFICCONT = "https://$DIVISION/a/service/jbrowse/organismSpecific/$ASSEMBLY"; #json file
my $RNASEQJUNCTIONS      = "https://$DIVISION/a/service/jbrowse/rnaseqJunctions/$ASSEMBLY"; #json file

##
# check for items to skip
##
my %skip_done;
while (<DATA>) {
    chomp;
    $skip_done{$_}++;
}

#this is the standard tracks.conf; ignores for the moment other track configs
#"/a/service/jbrowse/dnaseq/pfal3D7" BAMs/BWs
#"/a/service/jbrowse/rnaseq/pfal3D7" BAMs/BWs
#"/a/service/jbrowse/chipseq/pfal3D7" BWs

###
#fetch tracks.conf
###
my $TRACK_CONF = "tracks.conf";
system("curl --retry 5 -o $TRACK_CONF $TRACKINFO") == 0 or die;
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

###
#fetch other json config files
###
my $JUNC_FILE = $ASSEMBLY."_junc.json";
system("curl --retry 5 -o $JUNC_FILE $RNASEQJUNCTIONS") == 0 or die;
my $junc_blob;
{
    local $/ = undef;
    my $file = $JUNC_FILE;
    open FF, "<$file" or die "couldn't open $file: $!";
    $junc_blob = <FF>;
    close FF;
}
my $junc_tracks = JSON->new->decode($junc_blob);

my $SPEC_SPEC = $ASSEMBLY."_species_specific.json";
system("curl --retry 5 -o $SPEC_SPEC $ASSEMBLYSPECIFICCONT") == 0 or die;
my $ss_blob;
{
    local $/ = undef;
    my $file = $SPEC_SPEC;
    open FF, "<$file" or die "couldn't open $file: $!";
    $ss_blob = <FF>;
    close FF;
}
my $ss_tracks = JSON->new->decode($ss_blob);

###
#   Map the jbrowse json configs to the internal config hash
###
#print Dumper($$ss_tracks{'tracks'});
#print @{$$junc_tracks{'tracks'}};
#print "\n----------\n";
#print @{$$ss_tracks{'tracks'}};
push @{$$ss_tracks{'tracks'}}, @{$$junc_tracks{'tracks'}};
for my $hashref (@{$$ss_tracks{'tracks'}}) {
	#next unless defined $$hashref{'key'};
    my $track_key = $$hashref{'key'};
    next unless $$hashref{'storeClass'} =~ /SeqFeature\/REST/i;
    next if defined $$hashref{'subtracks'};

    for my $key (keys %{$hashref}) {
	    
            next if $key eq 'key';
	    if ($key eq 'query') {
                for my $qkey (keys %{$$hashref{'query'}}) {
                    $vuepath_track_info{$track_key}{"query.$qkey"} = $$hashref{'query'}{$qkey};
		    #print "$track_key\t$qkey\n";
	        }
                next;
	    }
            if ($key eq 'metadata') {
                for my $mkey (keys %{$$hashref{'metadata'}}) {
                    $vuepath_track_info{$track_key}{"metadata.$mkey"} = $$hashref{'metadata'}{$mkey};
                }
  	        next;
	    }

   	    $vuepath_track_info{$track_key}{$key} = $$hashref{$key};
    }
}

#print JSON->new->pretty->encode(\%vuepath_track_info);
#die;
###
#fetch contig/refseqs file
###
my $ASSEMBLY_CONTIG_FILE = $ASSEMBLY."_contigs.json";
system("curl --retry 5 -o $ASSEMBLY_CONTIG_FILE $SEQINFO") == 0 or die;

my $contiglist_blob;
{
    local $/ = undef;
    my $file = $ASSEMBLY_CONTIG_FILE;
    open FF, "<$file" or die "couldn't open $file: $!";
    $contiglist_blob = <FF>;
    close FF;
}
my $contig_json = JSON->new->decode($contiglist_blob);

###
# loop through all of the tracks
###

my $log = $ASSEMBLY."_tracks.log";
open LOG, ">$log" or die;
for my $tr_key (keys %vuepath_track_info) {
	#next unless (defined $vuepath_track_info{$tr_key}{'query.panId'});
    next if (defined $vuepath_track_info{$tr_key}{'query.feature'} &&
            $vuepath_track_info{$tr_key}{'query.feature'} eq 'ReferenceSequence');
    next unless $vuepath_track_info{$tr_key}{'storeClass'} =~ /SeqFeature\/REST/i;
    #   next unless defined $vuepath_track_info{$tr_key}{'query.edName'};

    #reasonably sure this isn't need anymore
    #  if (defined $vuepath_track_info{$tr_key}{'query.edname'}) {
    #    $vuepath_track_info{$tr_key}{'query.edName'} = $vuepath_track_info{$tr_key}{'query.edname'};
    #}

    my @GETarray;
    for my $config_key (keys %{$vuepath_track_info{$tr_key}}) {
        if ( $config_key =~ /query\.(.*)/ ) {
            my $subkey = $1;
            push @GETarray, "$subkey=".uri_escape($vuepath_track_info{$tr_key}{$config_key});
        }
    }
    my $GETstr     = join("&", @GETarray) if @GETarray >0 ;
    my $GETnamestr = join("_", @GETarray) if @GETarray >0 ;
    $GETnamestr   =~ s/\=/_/g;

    my $track_name    = uri_escape($vuepath_track_info{$tr_key}{'key'});
    $track_name = uri_escape($vuepath_track_info{$tr_key}{'label'}) unless $track_name;
    $track_name =~ s/\ /_/g;
    #my $query_feature = uri_escape($vuepath_track_info{$tr_key}{'query.feature'});
    #my $edname_str    = '';
    #if (defined $vuepath_track_info{$tr_key}{'query.edName'}) {
    #    $edname_str = '&edName='.uri_escape($vuepath_track_info{$tr_key}{'query.edName'});
    #} else {
    #    die Dumper($vuepath_track_info{$tr_key});
    #}
    #next unless $query_feature;
    #    warn $edname_str;

    my  $OUT = $ASSEMBLY.'_'.$track_name;
        $OUT .= "_$GETnamestr" if $GETnamestr;
	$OUT .= ".gff";
    next if $skip_done{$OUT}; # this one was completed on an earlier run
    if (!$track_name && !$GETnamestr) {
        print LOG "skipping because it will likely fail: $OUT\n";
	print LOG Dumper($vuepath_track_info{$tr_key}) . "\n";
	next;
    }

    open OUT, ">$OUT" or die "couldn't open $OUT for writing: $!"; 

    print OUT "##gff-version 3\n";
    print OUT "# VEuPathDB track config info:\n";
    for my $config_key (keys %{$vuepath_track_info{$tr_key}}) {
        print OUT "# $config_key = $vuepath_track_info{$tr_key}{$config_key}\n";
    }

    for my $contig_info (@{$contig_json}) {
        my $contig_name = $$contig_info{'name'};
	my $contig_end  = $$contig_info{'end'};
  
	my $json_file = $contig_name.'_'.$OUT.'.json';
        my $fetch_url = "https://$DIVISION/a/service/jbrowse/features/$contig_name?start=0&end=$contig_end";
	   $fetch_url .= '&'.$GETstr if $GETstr;
	my $curl = "curl --retry 5 -o $json_file \"$fetch_url\"";
        warn $curl;

        system( $curl ) == 0 or die $!;

        my $blob;
        {
        local $/ = undef;
        my $file = $json_file;
        open FF, "<$file" or die "couldn't open $file: $!";
        $blob = <FF>;
        close FF;
        unlink $file;
        }

	if ($blob =~ /^Cannot find/ ) {
            print LOG "got a 'cannot find' for $fetch_url\n";
	    next;
	}
        my $json = JSON->new->decode($blob) or print LOG "might die here: $fetch_url\n";
	next unless $json;

        for my $feature (@{$$json{'features'}}) {
            &parse_line(undef, $contig_name, $feature);
        }
	sleep 10; # 3 seconds between curls
    }
    close OUT;
    system("bzip2 $OUT");

    print LOG "$OUT\n";
    sleep 60; # 10 seconds between tracks
    #die;
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
    $gff3_feature{'type'}   = $$f{'type'} eq 'gff' ? defined $$f{'soterm'} ? $$f{'soterm'} : 'region' : $$f{'type'};
    $gff3_feature{'start'}  = $$f{'startm'}; # the base coord start
    $gff3_feature{'end'}    = $$f{'end'};
    $gff3_feature{'phase'}  = $$f{'phase'};
    if (!defined $$f{'strand'}) {
        $gff3_feature{'strand'} = '.';
    }
    elsif ($$f{'strand'} eq '1' or $$f{'strand'} eq '+1') {
        $gff3_feature{'strand'} = '+';
    }
    elsif ($$f{'strand'} eq '-1') {
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

# if we need to re run for an assembly, put the names of the completed GFF files here so they'll be skipped
__DATA__
