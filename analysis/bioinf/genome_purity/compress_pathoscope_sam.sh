#!/usr/bin/sh
#compressing pathoscope sam files

for i in SRR1555301 SRR1555302 SRR1555303 SRR1555304 SRR1555305 SRR1555306 SRR1555307 SRR1555308 SRR1555309 SRR1555310;
do
	 samtools view -b -S $i/$i\_1_tr.sam -o $i/$i\_1_tr.bam	
	 rm $i/$i\_1_tr.sam
done


for i in SRR1393710 SRR1393711 SRR1393713 SRR1393714 SRR1393716 SRR1393718 SRR1393719 SRR1393721;
do
	 samtools view -b -S $i/$i\_tr.sam -o $i/$i\_tr.bam	
	 rm $i/$i\_tr.sam
done