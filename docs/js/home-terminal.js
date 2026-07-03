// Home Terminal Showcase Logic

// Content definitions
const mainContentText = `
          ____ ____  _   _ __  __ ____  _    _ _ 
         / ___|  _ \\| | | |  \\/  |  _ \\| | _(_) |_ 
        | |  _| |_) | | | | |\\/| | | | | |/ / | __| 
        | |_| |  __/| |_| | |  | | |_| |   <| | |_ 
         \\____|_|    \\___/|_|  |_|____/|_|\\_\\_|\\__| 

         GPUMDkit Version 1.5.6 (dev) (2026-06-17) 
   Core Developer: Zihan YAN (yanzihan@westlake.edu.cn) 

 ----------------------- GPUMD ----------------------- 
  1) Format Conversion          2) Sample Structures 
  3) Workflow                   4) Calculators 
  5) Analyzer                   6) Visualization 
  7) Utilities                  8) Help 
  0) Exit 

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

const helpContentText = ` +-------------------------------------------------------------------------------------------------------+
 |                          GPUMDkit 1.5.6 (dev) (2026-06-17) Command Help                               |
 +-------------------------------------------------------------------------------------------------------+
 |                                          MAIN FUNCTIONS                                               |
 +-------------------------------------------------------------------------------------------------------+
 | -h            Show this help table            | -plt <type>        Plot and visualization tools       |
 | -calc <type>  Calculator tools                | -time <gpumd|nep>  Time-consuming analyzer            |
 | -update       Update GPUMDkit                 | -clean             Clean extra files in current dir   |
 +-------------------------------------------------------------------------------------------------------+
 |                                         FORMAT CONVERSION                                             |
 +-------------------------------------------------------------------------------------------------------+
 | -out2xyz      OUTCAR -> extxyz (shell)        | -out2exyz          OUTCAR -> extxyz (python)          |
 | -cp2k2xyz     CP2K log -> xyz                 | -xdat2exyz         XDATCAR -> extxyz                  |
 | -cif2pos      cif -> POSCAR                   | -cif2exyz          cif -> extxyz                      |
 | -pos2exyz     POSCAR -> extxyz                | -exyz2pos          extxyz -> POSCAR                   |
 | -pos2lmp      POSCAR -> LAMMPS data           | -lmp2exyz          LAMMPS dump -> extxyz              |
 | -traj2exyz    ASE traj -> extxyz              | -replicate         Replicate structure                |
 | -addgroup     Add group labels                | -addweight         Add structure weight in extxyz     |
 | -clean_xyz    Clean extra info in extxyz      | -get_frame         Extract specific frame             |
 | -frame_range  Extract frames by range         | -dp2xyz            DeepMD npy -> extxyz               |
 +-------------------------------------------------------------------------------------------------------+
 |                                            ANALYSIS                                                   |
 +-------------------------------------------------------------------------------------------------------+
 | -range        Energy/force/virial statistics  | -analyze_comp      Analyze composition                |
 | -chem_species Analyze chemical species        | -cbc               Charge balance check               |
 | -min_dist     Min distance (no PBC)           | -min_dist_pbc      Min distance with PBC              |
 | -filter_dist  Filter by min_dist (no PBC)     | -filter_dist_pbc   Filter by min_dist (PBC)           |
 | -pda          Probability density analysis    | -filter_box        Filter by box-edge length          |
 | -pynep        Deprecated PyNEP sampling       | -nep_modifier      Modify NEP model interactively     |
 +-------------------------------------------------------------------------------------------------------+
 | Detailed usage: gpumdkit.sh -<option> -h    Plot details: gpumdkit.sh -plt <type> -h                  |
 +-------------------------------------------------------------------------------------------------------+`;

const plotContentText = ` +-----------------------------------------------------------------------------------------------+
 |                     GPUMDkit 1.5.6 (dev) (2026-06-17) PLOT & VISUALIZATION TOOLS              |
 +-----------------------------------------------------------------------------------------------+
 |  Usage: gpumdkit.sh -plt <type>                        Help: gpumdkit.sh -plt <type> -h       |
 +-----------------------------------------------------------------------------------------------+
 |                                    NEP Training & Evaluation                                  |
 +-----------------------------------------------------------------------------------------------+
 |  train          - NEP training results           prediction     - NEP prediction results      |
 |  train_test     - NEP train and test results     parity_density - Parity density plot         |
 |  train_density  - Training results density plot  restart        - Parameters in nep.restart   |
 |  charge         - Charge distribution            born_charge    - Born effective charges      |
 |  dimer          - Dimer energy/force curve       force_errors   - Force errors                |
 |  des            - Descriptors                    lr             - Learning rate for gnep      |
 +-----------------------------------------------------------------------------------------------+
 |                                     Diffusion & Transport                                     |
 +-----------------------------------------------------------------------------------------------+
 |  msd            - Mean square displacement       msd_conv       - MSD convergence             |
 |  msd_all        - MSD for all species            sdc            - Self diffusion coefficient  |
 |  msd_sdc        - MSD and SDC together           sigma          - Arrhenius ionic conductivity|
 |  D              - Arrhenius diffusivity          sigma_xyz      - Directional Arrhenius sigma |
 |  D_xyz          - Directional Arrhenius D                                                     |
 |  doas           - Density of atomistic states                                                 |
 +-----------------------------------------------------------------------------------------------+
 |                                    MD & Structural Analysis                                   |
 +-----------------------------------------------------------------------------------------------+
 |  thermo         - thermo info in thermo.out      thermo2/3      - Thermo in different styles  |
 |  rdf            - Radial distribution function   rdf_pmf        - Potential of mean force     |
 |  vac            - Velocity autocorrelation       cohesive       - Cohesive energy curve       |
 |  net_force      - Net force distribution         plane-grid     - Displacement plane grid     |
 +-----------------------------------------------------------------------------------------------+
 |                                        Heat Transport                                         |
 +-----------------------------------------------------------------------------------------------+
 |  emd            - EMD results                    nemd           - NEMD results                |
 |  hnemd          - HNEMD results                  viscosity      - Viscosity                   |
 +-----------------------------------------------------------------------------------------------+
 |                                          Phonons                                              |
 +-----------------------------------------------------------------------------------------------+
 |  pdos           - VAC and PDOS                                                                |
 +-----------------------------------------------------------------------------------------------+
   See the codes in plt_scripts for more details 
   Code path: Scripts/plt_scripts`;

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
