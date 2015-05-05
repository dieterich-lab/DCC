export FILES=$1
export READFIX=$2

cd ${READFIX}

while read names
do


# First fix the strand of mate2
#cd ${READFIX}/${names}/STAR/mate2
~/tools/fixmate2.awk ${READFIX}/${names}/STAR/mate2/${names}_2Chimeric.out.junction > ${READFIX}/${names}/STAR/mate2/${names}_2Chimeric.out.junction.fixed
#cd ${READFIX}/${names}/STAR/
# Then fix the joined Chimeric.out.junction

~/tools/fixfunction.sh ${READFIX}/${names}/STAR/mate1/${names}_1Chimeric.out.junction ${READFIX}/${names}/STAR/mate2/${names}_2Chimeric.out.junction.fixed ${READFIX}/${names}/STAR/${names}Chimeric.out.junction

cd ${READFIX}
done < $1
