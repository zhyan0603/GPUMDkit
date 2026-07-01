<div align="center">
  <h1>🔁 Active Learning Workflow</h1>
  <p style="text-align: justify;">This tutorial explains the steps and considerations involved in the <code>workflow_active_learning_dev.sh</code> script. The workflow is used to generate, select, and analyze molecular structures, followed by simulations and data collection.</p>
</div>

**Script Location:** `Scripts/workflow/`

## Overview

This script automates one iteration of NEP active learning. The typical workflow includes:

| Step | Action | Purpose |
|------|--------|---------|
| 1 | Run MD simulations | Generate candidate structures with current NEP model |
| 2 | Filter structures | Remove unphysical configurations (short distances, large boxes) |
| 3 | Sample structures | Select diverse or uncertain structures for DFT |
| 4 | Run SCF calculations | Compute DFT reference energies/forces |
| 5 | NEP prediction | Check model accuracy on new structures |

Unlike menu-driven modules, this workflow is a shell script that you customize and run directly:

```bash
bash Scripts/workflow/workflow_active_learning_dev.sh
```

The sections below explain the main blocks in the script and what each part does.

## 1. **SLURM Directives**

```
#!/bin/bash -l
#SBATCH -p intel-sc3,intel-sc3-32c
#SBATCH -q huge
#SBATCH -N 1
#SBATCH -J workflow
#SBATCH -o workflow.log
#SBATCH --ntasks-per-node=1
```

The SLURM directives are used to define how the job will be submitted to the cluster:

- `#SBATCH -p` defines the partition, in this case, it's `intel-sc3` and `intel-sc3-32c`.
- `#SBATCH -q` specifies the queue, in this case, `huge`.
- `#SBATCH -N` allocates 1 node for the job.
- `#SBATCH -J` names the job `workflow`.
- `#SBATCH -o` defines the output log file as `workflow.log`.
- `--ntasks-per-node=1` specifies that only one task should run per node.

**<u>NOTE:</u>** If your machine does not have the SLURM environment, you can also run the script directly from the command line. For example:

```
nohup bash workflow_active_learning_dev.sh &>workflow.log &
```

## 2. **Basic Setup**

```
cd $SLURM_SUBMIT_DIR 
```

Ensure the working directory is correct, and that all necessary files are present for the job.

```
source ${GPUMDkit_path}/Scripts/workflow/submit_template.sh  # Load the submit template
python_pynep=python  # Python executable for the deprecated pynep branch if needed
```

- `GPUMDkit_path` is the environment variable that stores the path for the `GPUMDkit` .
- `python_pynep` points to the Python environment used for the deprecated `pynep`-related branch if needed.

------

## 3. **Variable Definitions**

```
work_dir=${PWD}  # Set the working directory
prefix_name=LiF_iter01  # Prefix for calculations
min_dist=1.4  # Minimum atom distance
box_limit=13  # Simulation box limit
max_fp_num=50  # Maximum number of single point calculations
sample_method=pynep  # Sampling method (options: uniform, random, pynep)
pynep_sample_dist=0.01  # Sampling distance for pynep
```

You can customize:

- `prefix_name` to reflect the name of your current work.
- `min_dist`, `box_limit`, and `max_fp_num` based on your own system.
- `sample_method` (choose between `uniform`, `random`, or `pynep`).
- `pynep_sample_dist` is the sampling distance for `pynep`.

------

## 4. **Check Required Files**

```
if [ -f nep.txt ] && [ -f nep.in ] && [ -f train.xyz ] && [ "$(find . -maxdepth 1 -name 'run_*.in' | wc -l)" -ge 1 ] && [ -f INCAR ] && [ -f POTCAR ] ; then
    # Check for the required files before proceeding.
else
    echo "Please put nep.in nep.txt train.xyz run_*.in (eg. run_1.in, run_2.in, ...) INCAR POTCAR [KPOINTS] and the sample_struct.xyz in the current directory."
    exit 1
fi
```

Make sure the necessary files (`nep.txt`, `nep.in`, `train.xyz`, `run_*.in` (e.g., `run_1.in`, `run_2.in`, ...), `INCAR`, `POTCAR`, `KPOINTS`) are available in the working directory. If any of these are missing, the script will terminate. `train.xyz` and `nep.txt` are always required as prerequisites for NEP model operation, and `run_*.in` files define the simulation parameters of MD in the current iteration.

------

## 5. **File Organization**

```
mkdir 00.md common
mv ${work_dir}/{nep.txt,nep.in,*.xyz,run_*.in,INCAR,KPOINTS,POTCAR} ./common
cp ${work_dir}/common/$sample_xyz_file ${work_dir}/00.md
```

The script organizes the working files into two folders:

- `00.md`: For molecular dynamics simulation.
- `common`: For shared resources such as `nep.txt`, `run_*.in` files, and the structure files, etc.

------

## 6. **Molecular Dynamics Simulation Submission**

```
submit_gpumd_array md ${sample_struct_num}
sbatch submit.slurm
```

After preparing the input files, the script submits an array of molecular dynamics (MD) tasks using `submit_gpumd_array`.

------

## 7. **Monitoring Task Completion**

```
while true; do
    logs=$(find "${work_dir}/00.md/" -type f -name log -path "*/sample_*/log")
    finished_tasks_md=$(grep "Finished running GPUMD." $logs | wc -l)
    error_tasks_md=$(grep "Error" $logs | wc -l)

    if [ "$error_tasks_md" -ne 0 ]; then
        echo "Error: MD simulation encountered an error."
        exit 1
    fi
    if [ $finished_tasks_md -eq $sample_struct_num ]; then
        break
    fi
    sleep 30
done
```

The script continuously checks whether all MD tasks have finished by searching for a `Finished running GPUMD.` message in the logs. If an error is encountered, the script terminates.

------

## 8. **Analysis and Filtering**

```
mkdir ${work_dir}/01.select
...
```

In the `01.select` folder, all structures in the trajectory file `dump.xyz` during the MD process will be analyzed and filtered to avoid the generation of non-physical structures as much as possible. Specifically, the `get_min_dist.py`, `filter_structures_by_distance.py` and `filter_exyz_by_box.py` scripts will be employed to check whether there are situations where the distance between atoms is too close or the simulated box exceeds the limit value, and such structures will be filtered out.

------

## 9. **Sampling Methods**

```
case $sample_method in
    "uniform")
    "random")
    "pynep")
```

The script supports three sampling methods: `uniform`, `random`, and `pynep`. It selects structures according to the chosen method and checks if the number exceeds `max_fp_num`, ensuring the final structure count is within limits.

------

## 10. **SCF Calculations**

```
submit_vasp_array scf ${selected_struct_num} ${prefix_name}
```

After sampling, SCF calculations are submitted using `submit_vasp_array`. The process is similar to the MD submission, but for SCF calculations.

------

## 11. **Prediction Step**

```
submit_nep_prediction
```

The final step involves submitting the NEP prediction task to check the accuracy of the NEP model. 



---

Thank you for using `GPUMDkit`! If you have any questions or need further assistance, feel free to open an issue on our GitHub repository or contact Zihan YAN (yanzihan@westlake.edu.cn).
