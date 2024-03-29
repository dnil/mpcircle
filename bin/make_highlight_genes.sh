#!/bin/bash

# usage: make_highlight_genes.sh [prefix]

# Needs a copy of the kgXref and knownGene tables from UCSC. Download eg as:
# curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/kgXref.txt.gz
# to the refs subdir or set the path in KGXREF  and KNOWNGENE accordingly.
# Or link from the refs subdir. 

# Daniel Nilsson, 2015

#KGXREF="/Users/daniel/src/FindTranslocations/scipts/kgXref.txt"
#KNOWNGENE="/Users/daniel/src/FindTranslocations/scipts/knownGene.txt"
KGXREF="refs/kgXref.txt"
KNOWNGENE="refs/knownGene.txt"

PREFIX=$1

if [ -z "$PREFIX" ]
then
  # default is all tabs in the variations directory
    PREFIX=""
fi

for sample_inter_tab in data/variations/$PREFIX*inter*tab; do
    samplebase=`basename $sample_inter_tab`;

    highlight_file=data/highlight/${samplebase%%_inter_chr_events_annotated_vsDB.tab}.highlight_genes.txt

    # empty highlight output
    echo -n > $highlight_file

   for gene in `perl -ne '@r=split(/\t+/,$_); if($r[0] < 3) {if (@r[13]=~m/gene_([^.]+.\d+)/) {print $1,"\n";} if (@r[14]=~m/gene_([^.]+.\d+)/) {print $1,"\n"; }}' $sample_inter_tab`;
    do
	grep -w $gene $KNOWNGENE | perl -ne '@r=split(/\t/,$_); @r[1]=~/chr([\dxymXYM]+)/; my $chrnr=$1; my $start=@r[3]; my $stop=@r[4]; if($chrnr ne "") { print "hs$chrnr $start $stop "; } else { exit 1 }' && grep -w $gene $KGXREF |perl -ne '@r=split(/\t/,$_); print @r[4], "\n";'
    done |sort|uniq >> $highlight_file

    sample_intra_tab=data/variations/${samplebase%%_inter_chr_events_annotated_vsDB.tab}_intra_chr_events_annotated_vsDB.tab
    for gene in `perl -ne '@r=split(/\t/,$_); if($r[0] < 3) {if (@r[13]=~m/gene_([^.]+.\d+)/) {print $1,"\n";} if (@r[14]=~m/gene_([^.]+.\d+)/) {print $1,"\n"; }}' $sample_intra_tab`;
    do
	grep -w $gene $KNOWNGENE | perl -ne '@r=split(/\t/,$_); @r[1]=~/chr([\dxymXYM]+)/; my $chrnr=$1; my $start=@r[3]; my $stop=@r[4]; if($chrnr ne "") { print "hs$chrnr $start $stop "; } else { exit 1 }' && grep -w $gene $KGXREF |perl -ne '@r=split(/\t/,$_); print @r[4], "\n";'
    done |sort|uniq >> $highlight_file
done

# highlight files without gene names.. 
for file in data/highlight/$PREFIX*.highlight_genes.txt ; do 
    cut -f1-3 -d\  $file > ${file%%.highlight_genes.txt}.highlight.txt
done
