#!/bin/bash
# Need two mates of paired end reads mapped seperately, then merge two Chimeric.out.junction files and select duplicate reads (means both two mates are chimera, and this is not output in the joined mapping), merge this with the joined mapping Chimeric.out.junction file

# Two mate Chimeric.out.junction files
export MATE1=$1
export MATE2=$2
export JOINED=$3
# Python script path
export PYPATH='/home/JCheng/tools'

#cat $1 $2 | awk '{if (x[$10]) { x_count[$10]++; print $0; if (x_count[$10] == 1) { print x[$10] } } x[$10] = $0}' | cat - $3 > $3.fixed

cat $1 $2 | ${PYPATH}/fixreadname.py - | awk '{if (x[$10]) { x_count[$10]++; print $0; if (x_count[$10] == 1) { print x[$10] } } x[$10] = $0}' | cat - $3 | ${PYPATH}/fixreadname.py - > $3.fixed
