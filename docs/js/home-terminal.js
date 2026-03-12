// Home Terminal Showcase Logic

// Content definitions
const mainContentText = `
          ____ ____  _   _ __  __ ____  _    _ _ 
         / ___|  _ \\| | | |  \\/  |  _ \\| | _(_) |_ 
        | |  _| |_) | | | | |\\/| | | | | |/ / | __| 
        | |_| |  __/| |_| | |  | | |_| |   <| | |_ 
         \\____|_|    \\___/|_|  |_|____/|_|\\_\\_|\\__| 
 
         GPUMDkit Version 1.5.2 (dev) (2025-03-10) 
   Core Developer: Zihan YAN (yanzihan@westlake.edu.cn) 
 
  ----------------------- GPUMD ----------------------- 
  1) Format Conversion          2) Sample Structures 
  3) Workflow                   4) Calculators 
  5) Analyzer                   6) Developing ... 
  0) Quit! 
  ------------>> 
  Input the function number: `;

const logsCommand1 = `
   Calling script by Yanzhou WANG et al. 
   Checking the convergence of OUTCARs ... 
   Total OUTCAR files: 400 
   Converged OUTCAR files: 397 
   Non-converged OUTCAR files: 3 
   Progress: [#####################] 100% (397/397) 
   Conversion complete. 
   dos2unix: converting train.xyz to Unix format... 
   Code path: Scripts/format_conversion/out2xyz.sh 
`;

const logsCommand2 = `
   Calling script by Boyi SITU. 
   Original atoms: 9, Target: 8000, Selected: 9x9x10 (7290 atoms) 
   Supercell created and saved to supercell.vasp 
   Code path: Scripts/format_conversion/replicate.py
`;

const helpContentText = ` +==================================================================================================+ 
 |                              GPUMDkit 1.5.2 (dev) (2025-03-10) Usage                             | 
 +======================================== Conversions =============================================+ 
 | -out2xyz       Convert OUTCAR to extxyz       | -pos2exyz     Convert POSCAR to extxyz           | 
 | -cif2pos       Convert cif to POSCAR          | -pos2lmp      Convert POSCAR to LAMMPS           | 
 | -cif2exyz      Convert cif to extxyz          | -lmp2exyz     Convert LAMMPS-dump to extxyz      | 
 | -addgroup      Add group label                | -addweight    Add weight to the struct in extxyz | 
 | -cp2k2xyz      Convert CP2K file to extxyz    | -traj2exyz    Convert ASE traj to extxyz         | 
 | -xdat2exyz     Convert XDATCAR to extxyz      | Developing...                                    | 
 +========================================= Analysis ===============================================+ 
 | -range         Print range of energy etc.     | -max_rmse     Get max RMSE from extxyz           | 
 | -min_dist      Get min_dist between atoms     | -min_dist_pbc Get min_dist considering PBC       | 
 | -filter_box    Filter struct by box limits    | -filter_value Filter struct by value (efs)       | 
 | -filter_dist   Filter struct by min_dist      | -analyze_comp Analyze composition of extxyz      | 
 | -pynep         Sample struct by pynep         | Developing...                                    | 
 +====================================== Misc Utilities ============================================+ 
 | -plt           Plot scripts                   | -get_frame     Extract the specified frame       | 
 | -calc          Calculators                    | -frame_range   Extract frames by fraction range  | 
 | -clean         Clear files for work_dir       | -clean_xyz     Clean extra info in XYZ file      | 
 | -time          Time consuming Analyzer        | -update        Update GPUMDkit                   | 
 +==================================================================================================+ 
 | For detailed usage and examples, use: gpumdkit.sh -<option> -h                                   | 
 +==================================================================================================+`;

const plotContentText = ` +=====================================================================================================+ 
 |                              GPUMDkit 1.5.2 (dev) (2025-03-10) Plotting Usage                       | 
 +=============================================== Plot Types ==========================================+ 
 | thermo          Plot thermo info                   | train          Plot NEP train results          | 
 | prediction      Plot NEP prediction results        | train_test     Plot NEP train and test results | 
 | msd             Plot mean square displacement      | msd_conv       Plot the convergence of MSD     | 
 | msd_all         Plot MSD of all species            | sdc            Plot self diffusion coefficient | 
 | rdf             Plot radial distribution function  | vac            Plot velocity autocorrelation   | 
 | restart         Plot parameters in nep.restart     | dimer          Plot dimer plot                 | 
 | force_errors    Plot force errors                  | des            Plot descriptors                | 
 | charge          Plot charge distribution           | lr             Plot learning rate              | 
 | doas            Plot density of atomistic states   | net_force      Plot net force distribution     | 
 | sigma           Plot Arrhenius sigma               | D              Plot Arrhenius diffusivity      | 
 | sigma_xyz       Plot directional Arrhenius sigma   | D_xyz          Plot directional Arrhenius D    | 
 | emd             Plot EMD results                   | nemd           Plot NEMD results               | 
 | hnemd           Plot HNEMD results                 | pdos           Plot VAC and PDOS               | 
 | plane-grid      Plot displacement plane grid       | parity_density Plot parity plot density        | 
 +=====================================================================================================+ 
 | For detailed usage and examples, use: gpumdkit.sh -plt <plot_type> -h                               | 
 +=====================================================================================================+ 
  See the codes in plt_scripts for more details 
  Code path: /d/Westlake/GPUMD/Gpumdkit/Scripts/plt_scripts`;

document.addEventListener('DOMContentLoaded', () => {
    // Render Main Content immediately
    const mainTermContent = document.getElementById('showcase-content-main');
    if (mainTermContent) {
        renderMain(mainTermContent);
    }

    // Render Logs Content immediately
    const logsTermContent = document.getElementById('showcase-content-logs');
    if (logsTermContent) {
        renderLogs(logsTermContent);
    }

    // Render Help Content immediately
    const helpTermContent = document.getElementById('showcase-content-help');
    if (helpTermContent) {
        renderHelp(helpTermContent);
    }

    // Render Plot Content immediately
    const plotTermContent = document.getElementById('showcase-content-plot');
    if (plotTermContent) {
        renderPlot(plotTermContent);
    }

    // Tab Switching Logic
    const tabs = document.querySelectorAll('.showcase-tab');
    tabs.forEach(tab => {
        tab.addEventListener('click', () => {
            // Remove active class from all tabs
            tabs.forEach(t => t.classList.remove('active'));
            // Add active class to clicked tab
            tab.classList.add('active');

            // Hide all content panels
            const panels = document.querySelectorAll('.showcase-panel');
            panels.forEach(p => p.classList.remove('active'));

            // Show target panel
            const targetId = tab.getAttribute('data-target');
            const targetPanel = document.getElementById(targetId);
            if (targetPanel) {
                targetPanel.classList.add('active');
            }
        });
    });
});

function renderMain(el) {
    el.innerHTML = '';
    const cmdDiv = document.createElement('div');
    cmdDiv.innerHTML = `<span style="color:#27ae60; font-weight:600;">user@gpumdkit:~$</span> <span style="font-weight:500;">gpumdkit.sh</span>`;
    el.appendChild(cmdDiv);
    
    const textNode = document.createTextNode(mainContentText);
    el.appendChild(textNode);

    // Add blinking cursor
    const cursorSpan = document.createElement('span');
    cursorSpan.className = 'cursor';
    cursorSpan.textContent = '_';
    el.appendChild(cursorSpan);
}

function renderLogs(el) {
    el.innerHTML = '';
    
    // Command 1
    const cmdDiv1 = document.createElement('div');
    cmdDiv1.innerHTML = `<span style="color:#27ae60; font-weight:600;">user@gpumdkit:~$</span> <span style="font-weight:500;">gpumdkit.sh -out2xyz .</span>`;
    el.appendChild(cmdDiv1);

    const outDiv1 = document.createElement('div');
    outDiv1.textContent = logsCommand1;
    outDiv1.style.color = 'var(--text-color)';
    outDiv1.style.marginBottom = '20px';
    el.appendChild(outDiv1);

    // Command 2
    const cmdDiv2 = document.createElement('div');
    cmdDiv2.innerHTML = `<span style="color:#27ae60; font-weight:600;">user@gpumdkit:~$</span> <span style="font-weight:500;">gpumdkit.sh -replicate unitcell.vasp supercell.vasp 8000</span>`;
    el.appendChild(cmdDiv2);

    const outDiv2 = document.createElement('div');
    outDiv2.textContent = logsCommand2;
    outDiv2.style.color = 'var(--text-color)';
    outDiv2.style.marginBottom = '20px';
    el.appendChild(outDiv2);

    // Prompt 3
    const cmdDiv3 = document.createElement('div');
    cmdDiv3.innerHTML = `<span style="color:#27ae60; font-weight:600;">user@gpumdkit:~$</span> <span class="cursor">_</span>`;
    el.appendChild(cmdDiv3);
}

function renderHelp(el) {
    el.innerHTML = '';
    const cmdDiv = document.createElement('div');
    cmdDiv.innerHTML = `<span style="color:#27ae60; font-weight:600;">user@gpumdkit:~$</span> <span style="font-weight:500;">gpumdkit.sh -h</span>`;
    el.appendChild(cmdDiv);
    
    const textNode = document.createTextNode(helpContentText);
    el.appendChild(textNode);
}

function renderPlot(el) {
    el.innerHTML = '';
    const cmdDiv = document.createElement('div');
    cmdDiv.innerHTML = `<span style="color:#27ae60; font-weight:600;">user@gpumdkit:~$</span> <span style="font-weight:500;">gpumdkit.sh -plt</span>`;
    el.appendChild(cmdDiv);
    
    const textNode = document.createTextNode(plotContentText);
    el.appendChild(textNode);
}
