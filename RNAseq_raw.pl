#!/usr/bin/perl -w

# Gets information about RNAseq coverages and splicing sites and stores that
# information into a hash. This version tries to unify RNAseq.pl, fill_model.pl
# and fill_model_hack.pl

use strict;
use Bio::DB::Sam;
use Storable;


#my $scount = 'samtools view -c ';
my $shead  = 'samtools view -H ';

my $bam = shift;
my $outfile = shift || '';
my @temp = split /[\/\\]/, $0;
my $pname = pop @temp;
die("Use: perl $pname bam_file\n")
  unless ($bam && -e $bam);
  

my $sam = Bio::DB::Bam->open($bam);


warn ("Retrieving file $bam...\n");

#print debug($h);

warn("Counting reads...\n");
my $fs = "$bam\.flagstat";
my $N = 0;
if (!-e $fs){
  `samtools flagstat $bam > $fs`;
}
open (IN, $fs) or die ('yo');
  while (<IN>){
    /(\d+)\s.+?mapped/;
    $N = $1 if ($1);
    last if ($N);
  }
close IN;
die("Failed $bam flagstat\n") unless ($N);

my $head = $sam->header();#`$shead $bam`;  
my $chrs = $head->target_name();
my $is_bchr = $chrs->[0] =~ /chr/i ? 1 : 0;

my $rfile = get_name($bam);

my $dcounter = 0;
my $r = {'splice' => {}};
my $cchr = '';
while (my $align = $sam->read1) {
  my $chr     = $chrs->[$align->tid];
  my $start     = $align->pos+1;
  my $end       = $align->calend;
  my $cigar     = $align->cigar_str;
  if (!$cchr || $chr ne $cchr){
    $cchr = $chr;
    warn("File $bam, chromosome $chr...\n");
  }
  my $fchr = uc($chr);
  $fchr =~ s/chr//i;
  $fchr = 'MT' if ($chr eq 'M');
  if (get_splice($r->{'splice'}, $align, $chr)){
    $dcounter++;
  }
}

$outfile = "$rfile\_raw\.hsh" unless ($outfile);
Storable::store($r, $outfile);


sub get_splice{
  my $splicing = shift;
  my $read = shift;
  my $chr = shift;
  my $cigar = $read->cigar_str();
  return 0 unless ($cigar =~ /N/);
  my $start = $read->pos() + 1;
  my @cig = $cigar =~ /(\d+\w)/g;
  my $current = $start;
  foreach my $code (@cig){
    my ($ncode, $lcode) = $code =~ /(\d+)(\w)/;
    warn("Unknown cigar code $lcode\n") unless ($lcode =~ /[MIDNSHP]/);
    if ($lcode eq 'N'){
      my $end = $current + $ncode;
      my $key = "$chr\:$current\_$end";
      $splicing->{$key}{'n'}++;
    }
    if ($lcode eq 'D' || $lcode eq 'M' ||
        $lcode eq 'S' || $lcode eq 'N'){
      $current += $ncode;
    }
  }
  return 1;
}

sub get_name{
  my $n = shift;
  my @y = split(/[\\\/]/, $n);
  my $r = pop @y;
  my @y2 = split(/\./, $r);
  my $r2 = shift @y2;
  return $r2;
}

sub shallow_copy{
  my $h1 = shift;
  my $h2 = shift;
  foreach my $key (keys %$h1){
    $h2->{$key} = $h1->{$key};
  }
  return 0;
}

sub min{
  my $result = shift;
  while(my $v = shift){
    $result = $v if ($v < $result);
  }
  return $result;
}

sub max{
  my $result = shift;
  while(my $v = shift){
    $result = $v if ($v > $result);
  }
  return $result;
}



sub debug{
  my $hash = shift;
  my $MAX = 10;
  my $i = 0;
  my $level = 0;
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
        foreach my $key (keys %$href){
          $i++;
          $uid = "\&__$i\__";
          $new_id_list->{$uid} = $href->{$key};
          my $safe = '';
          if (substr($key,0,1) =~ /[^a-zA-Z]/){
            $safe = 'a';
          }
          $temp_xml .= "\<$safe$key\>$uid\<\/$safe$key\>";
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
