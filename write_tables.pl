#!/usr/bin/perl -w
use strict;

use Storable;
use Cwd;

my $folder = getcwd();

my $table = {'files' => []};
my $retained = {};
warn("Reading splicing assignment...\n");
my $hsu = Storable::retrieve('/mnt/nas1/vqf/julia/cufflinks/hsu.hsh');

opendir(DIR, $folder);
while (my $f = readdir(DIR)){
  if ($f =~ /(.+?)_raw\.hsh$/){
    my $nroot = $1;
    warn("Reading $f...\n");
    my $t = Storable::retrieve($f);
    $table->{'cases'}{$nroot}{'hash'} = $t;
  }
  elsif ($f =~ /\.hsh$/){ #Retained
    $retained = Storable::retrieve($f);
  }
}
closedir DIR;

my @cases = sort keys %{$table->{'cases'}};
my $header = "Splice\tGene\t" . join("\t", @cases);

warn("Preparing list of positions\n");
my $total_sites = {};
foreach my $case (keys %{$table->{'cases'}}){
  warn("Adding $case...\n");
  foreach my $splice (keys(%{$table->{'cases'}{$case}{'hash'}{'splice'}})){
    my ($chr, $f, $t) = split(/[:_]/, $splice);
    $f--;
    my $nspl = "$chr:$f\-$t";
    #$f++;
    #$t--;
    $total_sites->{$splice} = {
      'nspl' => $nspl,
      'chr' => $chr,
      'from' => $f,
      'to' => $t
    };
  }
}

my $tablespl = "$header\n";
my $tableret = "$header\n";
my $tablerat = "$header\ttotal_reads\n";
my $table3p = "$header\ttotal_reads\n";

foreach my $spl (sort {($total_sites->{$a}{'chr'} cmp $total_sites->{$b}{'chr'}) ||
                        $total_sites->{$a}{'from'} <=> $total_sites->{$b}{'from'}} keys %$total_sites){
  my $lspl = '';
  my $lret = '';
  my $lrat = '';
  my $l3p = '';
  my $ingene = '-';
  my $nspl = $total_sites->{$spl}{'nspl'};
  my $f = $total_sites->{$spl}{'from'};
  my $t = $total_sites->{$spl}{'to'};
  if (exists($hsu->{'splices'}{$nspl}{'gene'}) && !exists($hsu->{'splices'}{$nspl}{'rep'})){
    $ingene = $hsu->{'splices'}{$nspl}{'gene'}{'gene_name'};
  }
  #elsif ((exists($hsu->{'5p'}{$f}) && $hsu->{'splices'}{'5p'}{$f}) || (exists($hsu->{'5p'}{$t}) && $hsu->{'splices'}{'5p'}{$t})){ # New 3p site
  #  $ingene = $hsu->{'splices'}{'5p'}{$f} if (exists($hsu->{'splices'}{'5p'}{$f}) && $hsu->{'splices'}{'5p'}{$f});
  #  $ingene = $hsu->{'splices'}{'5p'}{$t} if (exists($hsu->{'splices'}{'5p'}{$t}) && $hsu->{'splices'}{'5p'}{$t});
  #  #die("$f\t$t\t$hsu->{'5p'}{$f}\t-\t$hsu->{'5p'}{$t}\n") unless ($ingene);
  #  foreach my $c (@cases){
  #    my $s = $table->{'cases'}{$c}{'hash'}{'splice'}{$spl}{'n'} || 0;
  #    $l3p .= "$s\t";
  #  }
  #  $table3p .= "$spl\t$ingene\t$l3p\n"; 
  #}
  #else{
    my $total_reads = 0;
    foreach my $c (@cases){
      my $chr = $total_sites->{$spl}{'chr'};
      my $from = $total_sites->{$spl}{'from'};
      my $to = $total_sites->{$spl}{'to'};
      my $s = $table->{'cases'}{$c}{'hash'}{'splice'}{$spl}{'n'} || 0;
      $lspl .= "$s\t";
      my $f = 0; my $t = 0;
      $f = $retained->{$c}{$chr}{$from} if (exists($retained->{$c}{$chr}{$from}));
      $t = $retained->{$c}{$chr}{$to} if (exists($retained->{$c}{$chr}{$to}));
      my $r = $f + $t;
      $lret .= "$r\t";
      my $d = (log($s + 1) - log($r + 1))/log(2);
      $lrat .= "$d\t";
      $total_reads += $r + $s;
    }
    $tablespl .= "$spl\t$ingene\t$lspl\n";
    $tableret .= "$spl\t$ingene\t$lret\n";
    $tablerat .= "$spl\t$ingene\t$lrat\t$total_reads\n";
  #}
}

open (OUT, ">splicing.txt");
print OUT $tablespl;
close OUT;

open (OUT, ">retained.txt");
print OUT $tableret;
close OUT;

open (OUT, ">ratio.txt");
print OUT $tablerat;
close OUT;

open (OUT, ">alt3p.txt");
print OUT $table3p;
close OUT;


sub debug{
  my $hash = shift;
  my $i = 0;
  my $level = 0;
  my $MAX = 10;
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