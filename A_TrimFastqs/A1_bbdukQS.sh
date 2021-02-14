#! /bin/bash

for f in *.fastq.gz; do
	~/bbmap/bbduk.sh in=$f out=${f}_trimmed_clean  ref=ref/polyA.fa.gz,ref/truseq_rna.fa.gz k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 int=f -Xmx27g
	mv ${f}_trimmed_clean tc/$f
done	