#!/usr/bin/perl -w
use strict;

use Storable;
use Bio::DB::Sam;
use Cwd;

my $folder = getcwd();

my $table = {'files' => []};



opendir(DIR, $folder);
while (my $f = readdir(DIR)){
  if ($f =~ /(.+?)_raw\.hsh$/){
    my $nroot = $1;
    my $bfile = "$nroot\.bam";
    if (-e $bfile){
      warn("Reading $f...\n");
      my $t = Storable::retrieve($f);
      $table->{'cases'}{$nroot}{'hash'} = $t;
      $table->{'cases'}{$nroot}{'bam'} = Bio::DB::Bam->open($bfile);
    }
    else{
      warn("Could not find $bfile. It will not be considered\n");
    }
  }
}

closedir DIR;
warn("Preparing list of positions\n");
my $total_sites = {};
foreach my $case (keys %{$table->{'cases'}}){
  warn("Adding $case...\n");
  foreach my $splice (keys(%{$table->{'cases'}{$case}{'hash'}{'splice'}})){
    $total_sites->{$splice} = 1;
  }
}
print scalar(keys %$total_sites) . " splicing sites computed\n";

my @k = keys %$total_sites;
warn("Sorting keys\n");
my $coords = keysorter(\@k);
my $result = {'withspl' => {}, 'without' => {}};

warn("Computing numbers...\n");
my $rtable = {};
foreach my $case (keys %{$table->{'cases'}}){
  warn("\tCase $case...\n");
  my $ret_cov = {};
  my $bam = $table->{'cases'}{$case}{'bam'};
  my $retained = 0;
  my $cpos = 0;
  my $cchr = '';
  my $ccoords = [];
  my $head = $bam->header();#`$shead $bam`;  
  my $chrs = $head->target_name();
  my $cache = [];
  while (my $align = $bam->read1()){
    my $chr     = $chrs->[$align->tid];
    my $start     = $align->pos+1;
    my $end       = $align->calend;
    my $cigar     = $align->cigar_str;
    if (!$cchr || $chr ne $cchr){
      $cchr = $chr;
      $cpos = 0;
      if (!exists($coords->{$chr})){
        warn("No splicing sites in chr $chr\n");
      }
      @$ccoords = @{$coords->{$chr}};
      $retained = 0;
      warn("\t\tChromosome $chr...\n");
    }
    next if ($end < $cpos);
    next unless ($cigar =~ /^\d+M$/);
    while ($start > $cpos && scalar(@$ccoords)){
      $rtable->{$case}{$chr}{$cpos} = $retained;
      print "$chr\t$cpos\t$retained\n";
      $cpos = shift @$ccoords;
      if (scalar(@$cache)){
        if ($cache->[scalar(@$cache) - 1] < $cpos){
          $cache = [];        
        }
        else{
          my $cc = $cache->[0];
           while ($cc < $cpos){
            shift @$cache;
            $cc = $cache->[0];
          }
        }
      }
      $retained = scalar(@$cache) || 0;
    }
    if ($start <= $cpos && $end >= $cpos){
      $retained++;
      push @$cache, $end;
    }
    
  }
}
warn("Storing results\n");
Storable::store($rtable, "retained.hsh");

sub keysorter{
  my $in = shift;
  my $result = {};
  my $tmp = {};
  foreach my $c (@$in){
    my ($chr, $f, $t) = split(/[:_]/, $c);
    die("$c\t$chr\t$f\t$t\n") unless ($f);
    $tmp->{$chr}{$f} = 1;
    $tmp->{$chr}{$t - 1} = 1;
  }
  foreach my $chr (keys %$tmp){
    my @s = sort {$a <=> $b} keys %{$tmp->{$chr}};
    $result->{$chr} = \@s;
  }
  return $result;
}