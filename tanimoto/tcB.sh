DIR=/sc/arion/projects/H_filizm02a/soumyo/serotonin_receptor/structure
DIR1=/sc/arion/projects/H_filizm02a/ZINC20_3d/fingerprints/20210406
 
"$SCHRODINGER/utilities/canvasFPMatrix" -ifp2 $DIR1/BA_out.fp  -ifp $DIR/pimethixene.fp  -metric tanimoto -HOST Chimera-Long-lic03-2020-42:32 -TMPDIR /sc/arion/scratch/sens05/.tmp_schrodinger/ -JOB tcBA -ocsv $DIR/tanimoto/results/tc_BA.csv
"$SCHRODINGER/utilities/canvasFPMatrix" -ifp2 $DIR1/BB_out.fp  -ifp $DIR/pimethixene.fp  -metric tanimoto -HOST Chimera-Long-lic03-2020-42:32 -TMPDIR /sc/arion/scratch/sens05/.tmp_schrodinger/ -JOB tcBB -ocsv $DIR/tanimoto/results/tc_BB.csv
"$SCHRODINGER/utilities/canvasFPMatrix" -ifp2 $DIR1/BC_out.fp  -ifp $DIR/pimethixene.fp  -metric tanimoto -HOST Chimera-Long-lic03-2020-42:32 -TMPDIR /sc/arion/scratch/sens05/.tmp_schrodinger/ -JOB tcBC -ocsv $DIR/tanimoto/results/tc_BC.csv
"$SCHRODINGER/utilities/canvasFPMatrix" -ifp2 $DIR1/BD_out.fp  -ifp $DIR/pimethixene.fp  -metric tanimoto -HOST Chimera-Long-lic03-2020-42:32 -TMPDIR /sc/arion/scratch/sens05/.tmp_schrodinger/ -JOB tcBD -ocsv $DIR/tanimoto/results/tc_BD.csv
"$SCHRODINGER/utilities/canvasFPMatrix" -ifp2 $DIR1/BE_out.fp  -ifp $DIR/pimethixene.fp  -metric tanimoto -HOST Chimera-Long-lic03-2020-42:32 -TMPDIR /sc/arion/scratch/sens05/.tmp_schrodinger/ -JOB tcBE -ocsv $DIR/tanimoto/results/tc_BE.csv
"$SCHRODINGER/utilities/canvasFPMatrix" -ifp2 $DIR1/BF_out.fp  -ifp $DIR/pimethixene.fp  -metric tanimoto -HOST Chimera-Long-lic03-2020-42:32 -TMPDIR /sc/arion/scratch/sens05/.tmp_schrodinger/ -JOB tcBF -ocsv $DIR/tanimoto/results/tc_BF.csv
"$SCHRODINGER/utilities/canvasFPMatrix" -ifp2 $DIR1/BG_out.fp  -ifp $DIR/pimethixene.fp  -metric tanimoto -HOST Chimera-Long-lic03-2020-42:32 -TMPDIR /sc/arion/scratch/sens05/.tmp_schrodinger/ -JOB tcBG -ocsv $DIR/tanimoto/results/tc_BG.csv
"$SCHRODINGER/utilities/canvasFPMatrix" -ifp2 $DIR1/BH_out.fp  -ifp $DIR/pimethixene.fp  -metric tanimoto -HOST Chimera-Long-lic03-2020-42:32 -TMPDIR /sc/arion/scratch/sens05/.tmp_schrodinger/ -JOB tcBH -ocsv $DIR/tanimoto/results/tc_BH.csv
"$SCHRODINGER/utilities/canvasFPMatrix" -ifp2 $DIR1/BI_out.fp  -ifp $DIR/pimethixene.fp  -metric tanimoto -HOST Chimera-Long-lic03-2020-42:32 -TMPDIR /sc/arion/scratch/sens05/.tmp_schrodinger/ -JOB tcBI -ocsv $DIR/tanimoto/results/tc_BI.csv
"$SCHRODINGER/utilities/canvasFPMatrix" -ifp2 $DIR1/BJ_out.fp  -ifp $DIR/pimethixene.fp  -metric tanimoto -HOST Chimera-Long-lic03-2020-42:32 -TMPDIR /sc/arion/scratch/sens05/.tmp_schrodinger/ -JOB tcBJ -ocsv $DIR/tanimoto/results/tc_BJ.csv
