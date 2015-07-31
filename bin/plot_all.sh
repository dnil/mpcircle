#!/bin/bash
#
# Daniel Nilsson, 2015
#
# usage: plot_all.sh [prefix] [template_circos_config]

PREFIX=$1
TEMPLATE_CONF=$2 

if [ -z "$PREFIX" ]
then
  # default is to attempt plots for all files in the respecive driving data dirs (CNVnator output, variation tabs, and ROI-files)
    PREFIX=""
fi

if [ -z "$TEMPLATE_CONF" ] 
then
    TEMPLATE_CONF=etc/MP_template_w_genes_zoom_shrinkable.conf
fi

cd data/CNVnator

# convert cnvnator output to regular .seg files (as from ArrayCGH or such).
for file in ${PREFIX}*cnvnator.out;
do
    sample=${file%%.cnvnator.out};
    perl -e 'print "Sample\tChr\tStart\tStop\tCall\tPval\tRDnorm\n";' > $sample.cnvnator.seg;
    perl -ne 'chomp; @r=split(/\t/); ($chr,$start,$end)=$r[1]=~/([\w\d]+):(\d+)-(\d+)/; print join("\t","'$sample'",$chr, $start,$end, $r[0], $r[4], $r[3]-1),"\n"' $file >> $sample.cnvnator.seg ;
done

# join them for easy viz in IGV

cat ${PREFIX}*cnvnator.seg > ${PREFIX}all_cnvnator_joint.seg

head ${PREFIX}all_cnvnator_joint.seg > ${PREFIX}all_cnvnator_joint.fixed.seg 
grep -v ^Sample ${PREFIX}all_cnvnator_joint.seg >> ${PREFIX}all_cnvnator_joint.fixed.seg 

# Convert seg files to tracks that will display easiliy in e.g. circos.
for file in ${PREFIX}*cnvnator.seg;
do
    awk '{print $2, $3, $4, $7}' $file |sed -e's/chr/hs/g;' > ${file%%.seg}.track;
done

cd ../..

#for sample in data/variations/*inter*tab ; do samplebase=`basename $sample`; awk 'BEGIN {IFT='\t'} ($1==1) { print $2,$3,$4,$5,$6,$7 }' $sample |grep -v \_gl | grep -v chrUn |sed 's/chr/hs/g;' > data/FT_links/${samplebase%%.tab}.links ; done

# convert inter & intra FindTranslocation-tab files to links
for sample in data/variations/${PREFIX}*inter*tab ;
do 
    samplebase=`basename $sample`;
    (awk 'BEGIN {IFT='\t'} ($1<3) { print $2,$3,$4,$5,$6,$7 }' $sample |grep -v \_gl | grep -v chrUn |sed 's/chr/hs/g;'; samplebase=`basename $sample`; awk 'BEGIN {IFT='\t'} ($1<3) { print $2,$3,$4,$5,$6,$7 }' data/variations/${samplebase%%_inter_chr_events_annotated_vsDB.tab}_intra_chr_events_annotated_vsDB.tab |grep -v \_gl | grep -v chrUn |sed 's/chr/hs/g;';) > data/FT_links/${samplebase%%.tab}.links ; 
done

# Add gene highlights
./bin/make_highlight_genes.sh ${PREFIX}

# prepare config files based on an MP template config, CNVnator tracks and links.
for sample in data/ROI/${PREFIX}*.bed ; do 
    samplebase=`basename $sample`;
    sampleid=${samplebase%%.bed};

    # zoom plots to ROI-chrs
    zoomchrs=`cut -f 1 data/ROI/${sampleid}.bed |sort|uniq |perl -e '$out="hs"; while(<>) { chomp; if ($_ ne "") { if ($out ne "hs") {$out.=";hs".$_ ;} else { $out.=$_; } } } if ($out ne "hs") { print $out; }'`;

    sed -e 's/_LINKS_/'data\\/FT_links\\/${sampleid}_inter_chr_events_annotated_vsDB.links'/; s/_CNVNATORTRACK_/'data\\/CNVnator\\/$sampleid.cnvnator.track'/; s/_SAMPLEOUT_/'plots\\/$sampleid'/; s/_GENESLABEL_/'data\\/highlight\\/${sampleid}.highlight_genes.txt'/; s/_GENESNOLABEL_/'data\\/highlight\\/$sampleid.highlight.txt'/; s/_ZOOMCHRS_/'$zoomchrs'/;' $TEMPLATE_CONF  > plots/$sampleid.conf ;

    # run CIRCOS on the generated config files..
    circos -conf plots/$sampleid.conf -debug_group summary,timer > plots/$sampleid.run;
done

