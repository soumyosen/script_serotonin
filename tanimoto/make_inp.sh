
file="tcB.sh"
suffix="BA BB BC BD BE BF BG BH BI BJ"

echo "DIR=/sc/arion/projects/H_filizm02a/soumyo/serotonin_receptor/structure" >> $file
echo "DIR1=/sc/arion/projects/H_filizm02a/ZINC20_3d/fingerprints/20210406" >> $file
echo " " >> $file

for value in $suffix
do
echo "\"\$SCHRODINGER/utilities/canvasFPMatrix\" -ifp2 \$DIR1/${value}_out.fp  -ifp \$DIR/pimethixene.fp  -metric tanimoto -HOST Chimera-Long-lic03-2020-42:32 -TMPDIR /sc/arion/scratch/sens05/.tmp_schrodinger/ -JOB tc$value -ocsv \$DIR/tanimoto/results/tc_$value.csv" >> $file
done



