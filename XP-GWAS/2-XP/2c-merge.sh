#!/bin/bash

cp binomial.ac1 merged.binomial.ac.txt
for i in `seq 2 45`; do
	echo $i
	cat binomial.ac$i >> merged.binomial.ac.txt
done

mkdir 1-split_binomial
mv binomial.ac* 1-split_binomial

