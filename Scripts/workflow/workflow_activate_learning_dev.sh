#!bin/bash

work_dir=${PWD}
prefix_name=LiF_iter01  # prefix name for this workflow, used for the scf calculations
min_dist=1.4    # minimum distance between two atoms
box_limit=13    # box limit for the simulation box
max_fp_num=50  # maximum number of single point calculations
sample_method=pynep  # sampling method 'uniform' 'random' 'pynep'
pynep_sample_dist=0.01  # distance for pynep sampling

source ${GPUMDkit_path}/Scripts/workflow/submit_template.sh

#print some info
echo "********************************************" 
echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S")  
echo "WORK_DIR =" ${work_dir} 
echo "********************************************" 

# Check if the required files exist

if [ -f nep.txt ] && [ -f nep.in ] && [ -f train.xyz ] && [ -f run.in ] && [ -f INCAR ] && [ -f POTCAR ] && [ -f KPOINTS ]; then
    if [ $(find . -maxdepth 1 -name "*.xyz" | wc -l) -eq 2 ]; then
        sample_xyz_file=$(ls *.xyz | grep -v "train.xyz")
        sample_struct_num=$(grep -c Lat ${sample_xyz_file})
        echo "Found the exyz file: $sample_xyz_file"
        echo "There are ${sample_struct_num} structs in the ${sample_xyz_file}."
    else
        echo "Error: There should be exactly one exyz file (except for train.xyz) in the current directory."
        exit 1
    fi
    echo "All required files exist."
    echo "Starting the workflow:"
else
    echo "Please put nep.in nep.txt run.in INCAR POTCAR KPOINTS in the current directory."
fi

cd ${work_dir}
mkdir 00.modev 
cp $sample_xyz_file ${work_dir}/00.modev
cd ${work_dir}/00.modev
(echo 3; echo 302) | gpumdkit.sh >> /dev/null
cp ${work_dir}/{nep.in,nep.txt,run.in} ${work_dir}/00.modev/md
echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S") "Starting 00.modev step ..." 
submit_gpumd_array modev ${sample_struct_num}
sbatch submit.slurm
echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S") "${sample_struct_num} tasks had been submitted."

# Wait for all tasks to finish
while true; do
    finished_tasks_md=$(grep "Finished running GPUMD." sample_*/log | wc -l)
    error_tasks_md=$(grep "Error" sample_*/log | wc -l)

    if [ "$error_tasks_md" -ne 0 ]; then
        echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S") "Error: MD simulation encountered an error. Exiting..." 
        grep "Error" sample_*/log 
        exit 1
    fi
    if [ $finished_tasks_md -eq $sample_struct_num ]; then
        break
    fi
    sleep 30
done

echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S") "All modev tasks have finished. Starting analysis ..." 

mkdir ${work_dir}/01.select
cp ${work_dir}/train.xyz ${work_dir}/01.select
cat sample_*/dump.xyz >> ${work_dir}/01.select/modev_sampled_structs.xyz

echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S") "Analysis the min_dist in modev_sampled_structs.xyz" 
actual_min_dist=$(python ${GPUMDkit_path}/Scripts/sample_structures/get_min_dist.py ${work_dir}/01.select/modev_sampled_structs.xyz | awk '{print $4}')

if [ $(awk 'BEGIN {print ('$actual_min_dist' < '$min_dist')}') -eq 1 ]; then
    echo "The actual minimum distance ($actual_min_dist) between two atoms is less than the specified value ($min_dist)."
    echo "Filtering the structs based on the min_dist you specified."
    python ${GPUMDkit_path}/Scripts/sample_structures/filter_structs_by_distance.py ${work_dir}/01.select/modev_sampled_structs.xyz ${min_dist}
else
    echo "The actual minimum distance ($actual_min_dist) between two atoms is greater than the specified value ($min_dist)."
    echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S") "Analysis the box in filtered_modev_sampled_structs.xyz" 
    python ${GPUMDkit_path}/Scripts/analysis/filter_exyz_by_box.py ${work_dir}/01.select/filtered_modev_sampled_structs.xyz ${box_limit}
    echo " The box limit is $box_limit. filtered structs are saved in filtered_by_box.xyz"
fi

# Check the value of sample_method
case $sample_method in
    "uniform")
        echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S") "Performing uniform sampling..." 
        (echo 2; echo 201; echo "filtered_by_box.xyz ${sample_method} ${max_fp_num}") | gpumdkit.sh >> /dev/null
        mv sampled_structures.xyz selected.xyz
        selected_struct_num=$(grep -c Lat selected.xyz)
        ;;
    "random")
        echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S") "Performing random sampling..." 
        (echo 2; echo 201; echo "filtered_by_box.xyz ${sample_method} ${max_fp_num}") | gpumdkit.sh >> /dev/null
        mv sampled_structures.xyz selected.xyz
        selected_struct_num=$(grep -c Lat selected.xyz)
        ;;
    "pynep")
        echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S") "Performing pynep sampling..." 
        (echo 2; echo 202; echo "filtered_by_box.xyz train.xyz ${pynep_sample_dist}") | gpumdkit.sh >> /dev/null
        # Check the number of structures in selected.xyz
        selected_struct_num=$(grep -c Lat selected.xyz)
        if [ $selected_struct_num -gt $max_fp_num ]; then
            (echo 2; echo 201; echo "selected.xyz uniform ${max_fp_num}") | gpumdkit.sh >> /dev/null
        fi
        mv sampled_structures.xyz selected.xyz
        selected_struct_num=$(grep -c Lat selected.xyz)
        ;;
    *)
        echo "Invalid sample_method value. Please choose 'uniform', 'random', or 'pynep'." 
        exit 1
        ;;
esac

echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S") "The number of selected structures is $selected_struct_num." 

# Continue 02.scf step
echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S") "01.select have finished. Start 02.scf:" 
mkdir ${work_dir}/02.scf
cd ${work_dir}/02.scf
mv ${work_dir}/01.select/selected.xyz .
(echo 3; echo 301; echo "${prefix_name}") | gpumdkit.sh >> /dev/null
cp ${work_dir}/{INCAR,POTCAR,KPOINTS} ./fp
submit_vasp_array scf ${selected_struct_num} ${prefix_name}
echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S") "${selected_struct_num} tasks had been submitted."

# Wait for all tasks to finish
while true; do
    finished_tasks_scf=$(grep "F=" ${prefix_name}_*/log | wc -l)
    if [ $finished_tasks_scf -eq $selected_struct_num ]; then
        break
    fi
    sleep 30
done

echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S") "All scf tasks have finished. Starting prediction ..." 
echo "---------------------------------------"
gpumdkit.sh -out2xyz .
echo "---------------------------------------"
cd NEPdataset-multiple_frames
mkdir prediction
ln -s ${work_dir}/00.modev/md/nep.txt .
ln -s ../NEP-dataset.xyz ./train.xyz
cp ${work_dir}/00.modev/md/nep.in .

if ! grep -q "prediction" nep.in; then
    cat "prediction 1" >> nep.in
fi

submit_nep_prediction
gpumdkit.sh -plt prediction save
echo $(date -d "2 second" +"%Y-%m-%d %H:%M:%S") "Prediction finished."

