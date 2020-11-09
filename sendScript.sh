#!/bin/bash

# Submission script for GridEngine (GE). Each job will 
# be executed via the jobScript.sh
# This jobScript supports up to 7 parameters. Edit 
# the user specific part of the script according to 
# your program.
#
# Input to the script is a filelist with 1 file per line.
# For each file a job is started. With the parameter 
# nFilesPerJob a comma separated filelist will be 
# generated and handed to the job script. This feature
# is usefull when running many small jobs. Each
# job has its own logfile. All needed directories for the
# logfiles will be created if non existing.
#                                                                       

# IMPORTANT: the hera/prometheus cluster jobs will only
# see the /hera file system. All needed scripts, programs
# and parameters have to be located on /hera or the job
# will crash. This script syncs your working dir to the submission
# dir on /hera . Make sure your scripts use the submission dir!
# Software should be taken from /cvmfs/hades.gsi.de/install/
#
# job log files will be named like inputfiles. If nFilesPerJob > 1
# the log files will contain the partnumber.
#
######################################################################
#   CONFIGURATION           

user=$(whoami)
currentdir=$(pwd | xargs -i basename {}) 
currentDir=$(pwd)
target=${2}     #C or PE
momentum=${3}   #656,690,748,800
generation=${1} #gen0b, gen1 ...
submmissionbase=/lustre/nyx/hades/user/${user}/sub/exp_${4}

#submissiondir=${submmissionbase}/${currentdir}/program
#    outputdir=/lustre/nyx/hades/user/${user}/elscat/${generation}/${target}/${momentum}   # outputdir for files AND logFiles
#pathoutputlog=/lustre/nyx/hades/user/${user}/elscat/${generation}/logs                                 # protocol from batch farm for each file
submissiondir=${submmissionbase}/${currentdir}

#submissiondir=${submmissionbase}/${currentdir}/new_dieleAna
    outputdir=/lustre/nyx/hades/user/${user}/out_dilep/${generation}/${target}/${momentum}   # outputdir for files AND logFiles
#outputdir=/lustre/nyx/hades/user/${user}/TofCalibration/elscat_test/files
pathoutputlog=/lustre/nyx/hades/user/${user}/out_dilep/${generation}/logs                                 # protocol from batch farm for each file
#pathoutputlog=/lustre/nyx/hades/user/${user}/TofCalibration/elscat_test/logs

#nFilesPerJob=10
nFilesPerJob=10

# number of files to be analyzed by 1 job (default==1)
    jobscript=${submmissionbase}/${currentdir}/jobScript.sh
# exec script (full path, call without dot, set it executable!)
     filename=${target}${generation}${momentum}_${4}
   #  filename=${target}${generation}${momentum}_back_max1_fitted7_nomomcut
  # filename=${target}${generation}${momentum}_back_max1_fitted6_momcut
   	  # filename=${target}${generation}${momentum}_back_max1_fitted4_momcut

  #   filename=${target}${generation}${momentum}_back                            # filename of log file if nFilesPerJob > 1 (partnumber will be appended)
# filename of log file if nFilesPerJob > 1 (partnumber will be appended)
#filename=TofCal_test

par1=/cvmfs/hades.gsi.de/install/5.34.34/hydra2-4.9h/defall.sh     # optional par1 : environment script
par2=${submissiondir}/build/dieleAna                               # optional par2 : executable
par3=""                                                        	   # optional par3 : input file list
par4=${outputdir}                                                  # optional par4 : outputfile (part number will be appended (_num.root))
par5=1000000000                                                    # optional par5 : number of events
par6=""                                                            # optional par6
par7=""                                                           # optional par7
resources="--mem=2000 --time=0-4:00:00"                            # runtime < 10h, mem < 2GB
              
#filelist=${currentDir}/filelists/${generation}/${target}mom${momentum}_${generation}.txt # file list in local dir! not in submissiondir!!!
jobarrayFile="jobarray_${generation}_${target}_${momentum}.dat"
        
# jobarrayFile="jobarrayFile_${1}_${2}_${3}.dat"       
        
filelist=${currentDir}/filelists/${generation}/${target}_${momentum}.list # file list in local dir! not in submissiondir!!!
#filelist=${currentDir}/filelists/${generation}/${target}mom${momentum}_${generation}.txt # file list in local dir! not in submissiondir!!!
#filelist=/lustre/nyx/hades/user/prodrig/newdileptons/exp/completeSTAT/filelists/gen1/PEmom690_gen1_aug14.txt
#filelist=/lustre/nyx/hades/user/prodrig/newdileptons/exp/completeSTAT/new_dieleAna/PEmom690_gen1_aug14.txt
######################################################################

nFiles=$( cat $filelist | wc -l)

#---------------------------------------------------------------------
# create sub dirs
if [ ! -d $submmissionbase ]
then
    echo "===> CREATE SUBMISSIONBASEDIR : $submmissionbase"
    mkdir -p $submmissionbase
else
    echo "===> USE SUBMISSIONBASEDIR : $submmissionbase"
fi

#---------------------------------------------------------------------
# output dirs

if [ ! -d $outputdir ]
then
   echo "===> CREATE OUTPUTDIR : $outputdir"
   mkdir -p $outputdir
else
   echo "===> USE OUTPUTDIR : $outputdir"
fi

if [ ! -d $pathoutputlog ]
then
   echo "===> CREATE LOGDIR : $pathoutputlog"
   mkdir -p $pathoutputlog
else
   echo "===> USE LOGDIR : $pathoutputlog"
fi
#---------------------------------------------------------------------


ctF=0          # counter for file number
ctJ=0          # counter for job number
partNumber=0   # counter for part number

#---------------------------------------------------------------------
# read the files list into an job array
if [ -f $jobarrayFile ]
then
  rm -f $jobarrayFile
fi

echo "===> CREATING JOB ARRAY FILE"

declare -a jobarray
ct1=0
for file in $(cat $filelist)
do
   jobarray[$ct1]=$file
   ((ct1+=1))
done
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# loop over the job array and submit parts with
# nFilesPerJob to SLURM

while ((ctF<$nFiles))
do

     #---------------------------------------------------------------------
     # build comma separated file list
     # per job
     if [ $nFilesPerJob -gt 1 ]
     then
        infileList=${jobarray[${ctF}]}
        ((ctF+=1))
        for (( ctList=1;ctList<$nFilesPerJob; ctList++ ))
        do   	
            if [ $ctF -lt ${nFiles} ]
            then
               infileList="${infileList},${jobarray[${ctF}]}"
               ((ctF+=1))
            fi
        done
     else 
        infileList=${jobarray[${ctF}]}
        ((ctF+=1))
     fi

     ######################################################################
     #  SEND NEW JOB (USER SPECIFIC)
     
     par3=${infileList}

           #defall.sh prog filelist outdir  nev
     echo "${par1} ${par2} ${par3} ${par4} ${par5} ${par6} ${par7}" >>  $jobarrayFile
     

     ######################################################################
     
done
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# sync the local modified stuff 
# to the submission dir
echo "===> SYNC CURENTDIR TO SUBMISSIONDIR : rsync  -vHa $currentDir ${submmissionbase}"
rsync  -vHa $currentDir ${submmissionbase}/

syncStat=$?

if [ ! $syncStat -eq 0 ]
then
     echo "===> ERROR : SYNCHRONIZATION ENCOUNTERED PROBLEMS"
else

  echo "-------------------------------------------------"


  nFiles=$( cat $jobarrayFile | wc -l)
  ctsend=0
  block=500
  while ((${ctsend} * ${block} < ${nFiles}))
  do
     ((start=${ctsend}*${block}+1))
     ((stop= ${start}+${block}-1))
     ((rest=${nFiles}-${start}))
     if [ $rest -le $block ]
     then
        ((stop=$start+$rest))
     fi

     command="--array=${start}-${stop} ${resources} -D ${submissiondir}  --output=${pathoutputlog}/slurm-%A_%a.out ${jobscript} ${submissiondir}/${jobarrayFile} ${pathoutputlog} ${filename}"
     echo $command
     sbatch $command

     ((ctsend+=1))
  done

  echo "${nFiles} jobs for ${generation}, ${target} with ${momentum} MeV/c are submitted"
fi
