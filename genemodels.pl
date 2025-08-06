#!perl
use strict;
use warnings;
use local::lib;
use Data::Dumper;

use JSON;
use Bio::GFF3::LowLevel qw / gff3_format_feature  /;
use URI::Escape;

my $DIVISION = $ARGV[0];
my $STARTING_ORGSPEC = $ARGV[1];

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


#read provided organismspecific.json to get a list of dataset ids

my $starting_blob;
{
    local $/ = undef;
    my $file =  $STARTING_ORGSPEC;
    open FF, "<$file" or die "couldnt open $file: $!";
    $starting_blob = <FF>;
    close FF;
}
my $starting_json = JSON->new->decode($starting_blob);
my @datasetids = keys %{$$starting_json{datasets}};

#print $datasetids[0]."\n";
# loop through the datasetids and fetch the gene models and apollo additions.

for my $datasetid (@datasetids) {

    next if $skip_done{$datasetid};

    my $ASSEMBLY = $datasetid;
    my $SEQINFO              = "https://$DIVISION/a/service/jbrowse/seq/$ASSEMBLY";
    #    my $TRACKINFO            = "https://$DIVISION/a/jbrowse/tracks/$ASSEMBLY/tracks.conf";
    my $ASSEMBLYSPECIFICCONT = "https://$DIVISION/a/service/jbrowse/organismSpecific/$ASSEMBLY"; #json file
    #    my $RNASEQJUNCTIONS      = "https://$DIVISION/a/service/jbrowse/rnaseqJunctions/$ASSEMBLY"; #json file
    my $TRACKLIST            = "https://$DIVISION/a/service/jbrowse/tracks/$ASSEMBLY/trackList.json";

###
#fetch contig/refseqs file
###
    my $skip_contigs = 0;
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
    my $fetch_tracklist = ($contiglist_blob =~ /^HTTP 404 Not Found/) ? 1 : 0;

    my $contig_json;
    if (!$fetch_tracklist) { #contents probably json
        $contig_json = JSON->new->decode($contiglist_blob);
    }
    else { #contents not json, get the fai instead
        system("curl --retry 5 -o ".$ASSEMBLY."_tracklist.json $TRACKLIST") == 0 or die;

        my $tracklist_blob;
        {
            local $/ = undef;
            my $file = $ASSEMBLY."_tracklist.json";
            open FF, "<$file" or die "couldn't open $file: $!";
            $tracklist_blob = <FF>;
            close FF;
        }
        my $tracklist_json = JSON->new->decode($tracklist_blob);

        system("curl --retry 5 -o $ASSEMBLY.fai \"https://$DIVISION$$tracklist_json{refSeqs}\"") ;

        open FAI, "<$ASSEMBLY.fai" or die "couldn't open $ASSEMBLY.fai: $!";
        my $count = 0;
        while (<FAI>) {
            my @la = split("\t", $_);
            my %temp_hash;
            $temp_hash{'name'} = $la[0];
            $temp_hash{'end'}  = $la[1];
            push @{$contig_json}, \%temp_hash;
            $count++;
        }
        close FAI;
        $skip_contigs = 1 if $count >300;
    }

###
#fetch tracks.conf
###
  #    my $TRACK_CONF = "tracks.conf";
  #  system("curl --retry 5 -o $TRACK_CONF $TRACKINFO") == 0 or die;
  my %vuepath_track_info;
  #  my %lh;
  #  my $track_key;

  #  open TRACKCONF, "<$TRACK_CONF" or die "couldn't open $TRACK_CONF for reading: $!";
  #  while (<TRACKCONF>) {
  #      next if /^#/;
  #      chomp;
  #
  #      if (/^\[(.*)\]/) {
  #          $track_key = $1;
  #      }
  #      elsif (/^([^=]+)=(.*)/) {
  #         $vuepath_track_info{$track_key}{$1} = $2;
  #      }
  #
  #  }
  #  close TRACKCONF;

###
#fetch other json config files
###
  #  my $JUNC_FILE = $ASSEMBLY."_junc.json";
  #  system("curl --retry 5 -o $JUNC_FILE $RNASEQJUNCTIONS") == 0 or die;
  #  my $junc_blob;
  #  {
  #      local $/ = undef;
  #      my $file = $JUNC_FILE;
  #      open FF, "<$file" or die "couldn't open $file: $!";
  #      $junc_blob = <FF>;
  #      close FF;?log
  #
  #  }
  #  my $junc_tracks = JSON->new->decode($junc_blob);
  #
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
    my $log = $ASSEMBLY."_tracks.log";
    open LOG, ">$log" or die;
    if ($ss_blob =~ /^Cannot open/) {
        print LOG "failed to get species specific json\n" ;       
        next;
    }

    my $ss_tracks = JSON->new->decode($ss_blob) or warn "$SPEC_SPEC failed" ;

###
#   Map the jbrowse json configs to the internal config hash
###
#print Dumper($$ss_tracks{'tracks'});
#print @{$$junc_tracks{'tracks'}};
#print "\n----------\n";
#print @{$$ss_tracks{'tracks'}};
    
    #    push @{$$ss_tracks{'tracks'}}, @{$$junc_tracks{'tracks'}};
    for my $hashref (@{$$ss_tracks{'tracks'}}) {
	#next unless defined $$hashref{'key'};
        my $track_key = $$hashref{'key'};
        next unless ($track_key =~ /Annotated Transcripts/ 
                  or $track_key =~ /Community annotations from Apollo/);
        #      next unless $$hashref{'storeClass'} =~ /SeqFeature\/REST/i;
        #      next if defined $$hashref{'subtracks'};

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

    for my $tr_key (keys %vuepath_track_info) {
	#next unless (defined $vuepath_track_info{$tr_key}{'query.panId'});
        next if (defined $vuepath_track_info{$tr_key}{'query.feature'} &&
                $vuepath_track_info{$tr_key}{'query.feature'} eq 'ReferenceSequence');
              #next unless $vuepath_track_info{$tr_key}{'storeClass'} =~ /SeqFeature\/REST/i;
    #   next unless defined $vuepath_track_info{$tr_key}{'query.edName'};

    #reasonably sure this isn't need anymore
    #  if (defined $vuepath_track_info{$tr_key}{'query.edname'}) {
    #    $vuepath_track_info{$tr_key}{'query.edName'} = $vuepath_track_info{$tr_key}{'query.edname'};
    #}

        if ($vuepath_track_info{$tr_key}{'storeClass'} =~ /GFF3Tabix/) {
           my $filename = $ASSEMBLY . "_" . $tr_key . ".gz";
           my $url = "https://$DIVISION" . $vuepath_track_info{$tr_key}{'urlTemplate'} ;
           my $curl = "curl --retry 5 -o \"$filename\" \"$url\"";
           warn $curl;
           system($curl) == 0 or die "failed to fetch $url";
           next;
        }

        next if $skip_contigs;
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
        next if ($track_name =~ /Tandem%20Repeats/);
        next if ($track_name =~ /Gene%20Density/);
        next if ($track_name =~ /Low%20Complexity%20Regions/);
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
          	$OUT = substr($OUT, 0, 200) if length($OUT) > 200;
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
	          sleep 15; # 3 seconds between curls
        }
        close OUT;
        system("bzip2 $OUT");

        print LOG "$OUT\n";
        sleep 60; # 10 seconds between tracks
        #die;
    }

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
    $gff3_feature{'start'}  = $$f{'startm'} < 1 ? 1 : $$f{'startm'}; # the base coord start
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
