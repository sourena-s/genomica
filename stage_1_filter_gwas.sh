#!/bin/bash

rho_path=/projects/0/einf2700/sourena/genica-2025/gwas-databases/panukbb/sumstats/EUR

cat $rho_path/*.psy_rg.rho|head -n1

cat $rho_path/*.psy_rg.rho|awk 'function abs(x) {  if (x == "NA" || x == "" || x ~ /[^0-9eE+.-]/ || x > 1.1 || x < -1.1) x = 0 ;return (x<0 ? -x: x) }  
{thresh=0.5
	if ( abs($3) > thresh || abs($4)> thresh || abs($6) > thresh || abs($8) > thresh || abs($9) > thresh || abs($10) > thresh || abs($11) > thresh || abs($12) > thresh || abs($13) > thresh ) print $0;} '
