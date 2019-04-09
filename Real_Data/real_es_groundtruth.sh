 ##################################
 ## CROSSVAL GROUND TRUTH ES
 ##################################

#PATH1: directory where your SPMs, design files and subject masks are located

for I in {1..10}
do

    #Copy all files to target file where flameo results of this iteration will be stored
    mkdir PATH1/groundtruth_$I
    cp PATH1/varcope* PATH1/groundtruth_$I
    cp PATH1/cope* PATH1/groundtruth_$I
    cp PATH1/subj_1_OVERALLMASK.nii PATH1/groundtruth_$I
    cp PATH1/design.* PATH1/groundtruth_$I

    #Remove one file for the leave-one-out cross validation, step I
    cd ~PATH1/groundtruth_$I
    rm varcope_$I\.nii
    rm cope_$I\.nii

    #Merge all needed files in 4D
    fslmerge -t allcopes.nii cope*
    fslmerge -t allvarcopes.nii varcope*

    #Do group analysis
    flameo --cope=allcopes --vc=allvarcopes --mask=subj_1_OVERALLMASK --ld=stats --dm=design.mat --cs=design.grp --tc=design.con --runmode=fe

    #Copy relevant files
    cp PATH1/groundtruth_$I/stats/cope1.nii.gz PATH1/truth_ES_$I.nii.gz

    #Remove all other files
    rm -r PATH1/groundtruth_$I

done
