#!/bin/tcsh
echo Test Script


# Number of files to be used by each job
set files_per_job = 5

set filenum = 0
# set jobnum = 0

rm -f fileinp_$filenum
rm -f histograms_MC_Data$filenum.C
rm -f jobs_submit

set count = 0
foreach line ("`cat h2tauFall13.txt`")
 @ count = $count + 1
 
	echo "$line" >> fileinp_$filenum
#	echo $count
	if ($count % $files_per_job == 0) then
	 rm -f test_data$filenum.py
         foreach line2 ("`cat template.src`")
	 echo "$line2"
            if ( "$line2" == " #PALLE") then
             cat fileinp_$filenum >> jobs2_$filenum
            else
             echo "$line2" >> jobs2_$filenum
            endif 
         end 
         rm -f fileinp_$filenum
         sed s/treeXX/tree$filenum/ jobs2_$filenum >> test_data$filenum.py
         rm -f jobs2_$filenum
         echo "qsub -f -q" jobs_$filenum  " >> job_output_"$filenum>> jobs_submit
         
	 @ filenum = $filenum + 1
	  rm -f fileinp_$filenum  
	endif
end
if ( -r fileinp_$filenum) then
	 rm -f test_data$filenum.py
         foreach line2 ("`cat template.src`")        
           if ( "$line2" == " #PALLE ") then
            cat fileinp_$filenum >> jobs2_$filenum
           else
            echo "$line2" >> jobs2_$filenum
           endif 
         end
         sed s/treeXX/tree$filenum/ jobs2_$filenum >> test_data$filenum.py
         rm -f jobs2_$filenum
         rm -f fileinp_$filenum
         echo "qsub -f -q" jobs_$filenum  " >> job_output_"$filenum>> jobs_submit
endif

echo " Jobs produced"
