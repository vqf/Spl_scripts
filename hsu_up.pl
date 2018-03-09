#!/usr/bin/perl -w
use strict;
use Storable;

my $result = {};
my $debug = {};

my $WN = 100000;
my $ln = 0;
warn("Reading stdin...\n");
while (<>){
  next if (/^#/);
  chomp;
  my @f = split(/\t/);
  my $data = pop @f;
  my @d = split(/;/, $data);
  my $tmp = {};
  foreach my $d (@d){
    my ($k, $v) = $d =~ /(\S+)\s+\"([^\"]+)\"/;
    $tmp->{$k} = $v;
  }
  my $chr = $f[0];
  my $f = min($f[3], $f[4]);
  my $t = max($f[3], $f[4]);
  my $type = $f[2];
  my $gname = $tmp->{'gene_name'};
  next unless ($gname);
  if ($type eq 'gene'){
    my $strand = 1;
    if ($f[6] eq '-'){
      $strand = -1;
    }
    $result->{'genes'}{$gname} = $tmp;
    $result->{'genes'}{$gname}{'chr'} = $chr;
    $result->{'genes'}{$gname}{'from'} = $f;
    $result->{'genes'}{$gname}{'to'} = $t;
    $result->{'genes'}{$gname}{'strand'} = $strand;
  }
  elsif($type eq 'transcript'){
    my $tid = $tmp->{'transcript_id'};
    $result->{'genes'}{$gname}{'transcripts'}{$tid} = $tmp;
    $result->{'genes'}{$gname}{'transcripts'}{$tid}{'chr'} = $chr;
    $result->{'genes'}{$gname}{'transcripts'}{$tid}{'from'} = $f;
    $result->{'genes'}{$gname}{'transcripts'}{$tid}{'to'} = $t;
  }
  elsif($type eq 'exon'){
    my $tid = $tmp->{'transcript_id'};
    my $en = $tmp->{'exon_number'};
    $result->{'genes'}{$gname}{'transcripts'}{$tid}{'exons'}{$en} = $tmp;
    $result->{'genes'}{$gname}{'transcripts'}{$tid}{'exons'}{$en}{'chr'} = $chr;
    $result->{'genes'}{$gname}{'transcripts'}{$tid}{'exons'}{$en}{'from'} = $f;
    $result->{'genes'}{$gname}{'transcripts'}{$tid}{'exons'}{$en}{'to'} = $t;
  }
  $debug->{$f[1]} = 1 unless exists($debug->{$f[1]});
  $ln++;
  if ($ln > $WN){
    $ln = 0;
    warn("\t$WN more lines read\n");
    #last;
  }
}
warn("Considering splicing sites...\n");
foreach my $gene (keys %{$result->{'genes'}}){
  my $common_spl = {};
  my $strand = $result->{'genes'}{$gene}{'strand'};
  my $chr = $result->{'genes'}{$gene}{'chr'};
  my $nt = scalar(keys %{$result->{'genes'}{$gene}{'transcripts'}});
  foreach my $transcript (keys %{$result->{'genes'}{$gene}{'transcripts'}}){
    next unless(scalar (keys %{$result->{'genes'}{$gene}{'transcripts'}{$transcript}{'exons'}}) > 1);
    my $coor = [];
    foreach my $exon (keys %{$result->{'genes'}{$gene}{'transcripts'}{$transcript}{'exons'}}){
      push @$coor, $result->{'genes'}{$gene}{'transcripts'}{$transcript}{'exons'}{$exon}{'from'};
      push @$coor, $result->{'genes'}{$gene}{'transcripts'}{$transcript}{'exons'}{$exon}{'to'};
    }
    my @sc = sort{$a <=> $b} @$coor;
    shift @sc; pop @sc;
    while (scalar @sc){
      my $f = shift @sc; my $t = shift @sc;
      my $k = "$chr\:$f\-$t";
      my $p5 = ($strand == -1) ? $t : $f;
      my $p3 = ($strand == -1) ? $f : $t;
      $common_spl->{$k}{'n'}++;
      $common_spl->{$k}{'5p'} = $p5;
      $common_spl->{$k}{'3p'} = $p3;
    }
  }
  foreach my $s (keys %$common_spl){
    if ($common_spl->{$s}{'n'} == $nt){
      if (exists($result->{'splices'}{$s})){
        warn("Discarding site $s. It belongs to $result->{'splices'}{$s}{'gene'}{'gene_name'} and $gene\n");
        $result->{'splices'}{$s}{'rep'} = 1;
      }
      else{
        $result->{'splices'}{$s}{'gene'} = $result->{'genes'}{$gene};
        my $p5 = $common_spl->{$s}{'5p'};
        my $p3 = $common_spl->{$s}{'3p'};
        $result->{'5p'}{$p5}{'gene'} = $gene;
        $result->{'5p'}{$p5}{'3p'}{$p3} = 1;
        $result->{'3p'}{$p3}{'gene'} = $gene;
        $result->{'3p'}{$p3}{'5p'}{$p5} = 1;
      }
    }
  }
}
Storable::store($result, 'hsu.hsh');

sub debug{
  my $hash = shift;
  my $i = 0;
  my $level = 0;
  my $MAX = 30;
  my $uid = "\&__$i\__";
  my $xml = '<root>'.$uid.'</root>';
  my $id_list = {$uid => $hash};
  while (scalar keys %$id_list){
    my $new_id_list = {};
    $level++;
    foreach my $id (keys %$id_list){
      my $temp_xml = '';
      my $href = $id_list->{$id};
      if (ref($href) eq 'ARRAY'){
        my $counter = 0;
        foreach my $val (@$href){
          $i++;
          $uid = "\&__$i\__";
          $new_id_list->{$uid} = $val;
          $temp_xml .= "\<c_$level\_$counter\>$uid\<\/c_$level\_$counter\>";
          $counter++;
          last if ($counter > $MAX);
        }
      }
      elsif (ref($href) eq 'HASH' || ref($href) eq __PACKAGE__){
        my $counter = 0;
        foreach my $key (keys %$href){
          $i++;
          $uid = "\&__$i\__";
          $new_id_list->{$uid} = $href->{$key};
          my $safe = '';
          if (substr($key,0,1) =~ /[^a-zA-Z]/){
            $safe = 'a';
          }
          $temp_xml .= "\<$safe$key\>$uid\<\/$safe$key\>";
          $counter++;
          last if ($counter > $MAX);
        }
      }
      else{
        $href =~ s/[<>]//g;
        $href = '<![CDATA['.$href.']]>';
        $temp_xml .= $href;
      }
      $temp_xml = '_empty_' unless ($temp_xml);
      die ("$id\t$temp_xml\n") unless ($xml =~ /$id/);
      $xml =~ s/$id/$temp_xml/;
    }
    $id_list = $new_id_list;
  }
  return $xml;
}

sub min{
  my $result = shift;
  foreach my $v (@_){
    $result = $v if ($v < $result);
  }
  return $result;
}

sub max{
  my $result = shift;
  foreach my $v (@_){
    $result = $v if ($v > $result);
  }
  return $result;
}