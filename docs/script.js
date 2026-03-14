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

const logsContentText = ` 
 $ gpumdkit.sh -out2xyz . 
   Calling script by Yanzhou WANG et al. 
   Checking the convergence of OUTCARs ... 
   Total OUTCAR files: 400 
   Converged OUTCAR files: 397 
   Non-converged OUTCAR files: 3 
   Progress: [#####################] 100% (397/397) 
   Conversion complete. 
   dos2unix: converting train.xyz to Unix format... 
   Code path: Scripts/format_conversion/out2xyz.sh 
 
 $ gpumdkit.sh -replicate unitcell.vasp supercell.vasp 8000 
   Calling script by Boyi SITU. 
   Original atoms: 9, Target: 8000, Selected: 9x9x10 (7290 atoms) 
   Supercell created and saved to supercell.vasp 
   Code path: Scripts/format_conversion/replicate.py`;

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

const analysisContentText = ` +==================================================================================================+ 
 |                              GPUMDkit Analysis Tools Details                                     | 
 +==================================================================================================+ 
 | Command: gpumdkit.sh -range                                                                      | 
 | Description: Calculates min, max, mean, and std dev of energy, force, and stress components.     | 
 |                                                                                                  | 
 | Command: gpumdkit.sh -min_dist                                                                   | 
 | Description: Computes minimum interatomic distances for all pairs in the structure.              | 
 |                                                                                                  | 
 | Command: gpumdkit.sh -max_rmse                                                                   | 
 | Description: Identifies the frame with the highest Root Mean Square Error in force prediction.   | 
 |                                                                                                  | 
 | Command: gpumdkit.sh -analyze_comp                                                               | 
 | Description: Analyzes the chemical composition of the system, useful for multi-component alloys. | 
 +==================================================================================================+`;

const miscContentText = ` +==================================================================================================+ 
 |                              GPUMDkit Miscellaneous Utilities                                    | 
 +==================================================================================================+ 
 | Command: gpumdkit.sh -clean                                                                      | 
 | Description: Removes temporary files (e.g., dump.*, log.lammps) to clean the workspace.          | 
 |                                                                                                  | 
 | Command: gpumdkit.sh -get_frame <index>                                                          | 
 | Description: Extracts a specific frame (0-based index) from a trajectory file.                   | 
 |                                                                                                  | 
 | Command: gpumdkit.sh -update                                                                     | 
 | Description: Pulls the latest changes from the remote repository to update GPUMDkit.             | 
 +==================================================================================================+`;

// DOM Elements
let terminals = [];
let contents = [];

// Get static elements
const container = document.querySelector('.container');
const scrollBtn = document.getElementById('scroll-btn');

// State for Parallax
let mouseX = 0;
let mouseY = 0;
let windowWidth = window.innerWidth;
let windowHeight = window.innerHeight;

// Navigation Function
function navigateToHome() {
    window.location.href = 'home.html';
}

// Initialization
function init() {
    // Initialize terminal array
    terminals = [
        document.getElementById('terminal-misc'),      // 0
        document.getElementById('terminal-analysis'),  // 1
        document.getElementById('terminal-plot'),      // 2
        document.getElementById('terminal-help'),      // 3
        document.getElementById('terminal-main'),      // 4
        document.getElementById('terminal-logs')       // 5
    ];

    // Initialize contents array
    contents = [
        { el: document.getElementById('content-misc'), text: miscContentText, cmd: '$ gpumdkit.sh -misc' },
        { el: document.getElementById('content-analysis'), text: analysisContentText, cmd: '$ gpumdkit.sh -ana' },
        { el: document.getElementById('content-plot'), text: plotContentText, cmd: '$ gpumdkit.sh -plt' },
        { el: document.getElementById('content-help'), text: helpContentText, cmd: '$ gpumdkit.sh -h' },
        { el: document.getElementById('content-main'), text: mainContentText, cmd: '$ gpumdkit.sh' },
        { el: document.getElementById('content-logs'), text: logsContentText, cmd: 'user@gpumdkit:~$ ' }
    ];

    // Safety check
    if (terminals.some(t => !t)) {
        console.error("Critical: Some terminal elements are missing from DOM.");
        return;
    }

    // 1. Initial State
    terminals.forEach(t => {
        t.style.opacity = '0';
        t.style.transform = 'translate(-50%, -50%) scale(0.8)';
    });

    // 2. Render Content
    contents.forEach((c, i) => {
        if (!c.el) return;
        if (i === 4 || i === 5) return; // Skip dynamic ones for now
        
        c.el.innerHTML = '';
        const cmdDiv = document.createElement('div');
        cmdDiv.innerHTML = `<span style="color:#27ae60; font-weight:600;">user@gpumdkit:~$</span> <span style="font-weight:500;">${c.cmd}</span>`;
        c.el.appendChild(cmdDiv);
        const outDiv = document.createElement('div');
        outDiv.textContent = c.text;
        c.el.appendChild(outDiv);
    });

    // 3. Start Directly (No Fan Animation)
    collapseAndShowFocus();
    
    // Show nav controls faster
    const navControls = document.querySelector('.terminal-nav-controls');
    if(navControls) {
        navControls.style.animationDelay = '0.5s';
        navControls.style.opacity = '1'; 
    }

    // 4. Setup Listeners
    setupEventListeners();
}

function setupEventListeners() {
    // Mouse Parallax (Reduced for cleaner look)
    window.addEventListener('mousemove', handleMouseMove);
    window.addEventListener('resize', () => {
        windowWidth = window.innerWidth;
        windowHeight = window.innerHeight;
    });

    // Tab Switching (Now Pill Switching)
    const pills = document.querySelectorAll('.nav-pill');
    pills.forEach(pill => {
        pill.addEventListener('click', (e) => {
            e.stopPropagation(); // Prevent conflicts
            
            // Update Active State
            pills.forEach(p => p.classList.remove('active'));
            pill.classList.add('active');

            // Switch Card
            const targetId = pill.getAttribute('data-target');
            const targetTerm = document.getElementById(targetId);
            if (targetTerm && !targetTerm.classList.contains('active-front')) {
                switchCard(targetTerm);
            }
        });
    });

    // Scroll Indicator
    if (scrollBtn) {
        scrollBtn.addEventListener('click', navigateToHome);
    }

    // Keyboard Navigation
    document.addEventListener('keydown', (e) => {
        if (e.key === 'Enter') {
            navigateToHome();
        }
    });

    // Mouse Wheel
    let wheelTimeout;
    window.addEventListener('wheel', (e) => {
        if (e.deltaY > 0) { // Scrolling down
            if (!wheelTimeout) {
                wheelTimeout = setTimeout(() => {
                    navigateToHome();
                }, 200);
            }
        }
    });
}

function switchCard(newFrontTerm) {
    // Find current front terminal
    const currentFront = document.querySelector('.terminal.active-front');
    
    // Quick swap logic without fade delays
    // 1. Move old front to stacked-right (background)
    if (currentFront && currentFront !== newFrontTerm) {
        currentFront.classList.remove('active-front');
        currentFront.classList.add('stacked-right');
        // Force style update for immediate effect if needed, but CSS transition handles move
        currentFront.style.zIndex = 40;
        currentFront.style.opacity = '1';
        currentFront.style.filter = 'brightness(0.97)';
    }

    // 2. Handle all others - keep them in back or move out of way
    terminals.forEach(t => {
        if (t === newFrontTerm || t === currentFront) return;
        
        // If it was stacked-right, move it to back or left to make room?
        // Let's just push everything else to stacked-back for cleanliness
        t.classList.remove('active-front');
        t.classList.remove('stacked-right');
        t.classList.remove('stacked-left');
        t.classList.add('stacked-back');
        
        t.style.zIndex = 10;
        t.style.opacity = '0';
    });

    // 3. Bring new terminal to front
    newFrontTerm.classList.remove('stacked-right');
    newFrontTerm.classList.remove('stacked-left');
    newFrontTerm.classList.remove('stacked-back');
    
    newFrontTerm.classList.add('active-front');
    newFrontTerm.style.zIndex = 50;
    newFrontTerm.style.opacity = '1';
    newFrontTerm.style.filter = 'brightness(1)';
    newFrontTerm.style.transform = ''; 
    
    // If it's Logs, ensure content is rendered
    if (newFrontTerm.id === 'terminal-logs' && newFrontTerm.querySelector('.terminal-content').innerHTML === '') {
        typeLogsTerminal();
    }
}

function handleMouseMove(e) {
    // Debounce or throttle check could be added here for performance
    // Simple optimization: only calculate if frame is ready
    if (!ticking) {
        window.requestAnimationFrame(() => {
            mouseX = (e.clientX - windowWidth / 2) / (windowWidth / 2); 
            mouseY = (e.clientY - windowHeight / 2) / (windowHeight / 2);
            updateParallax();
            ticking = false;
        });
        ticking = true;
    }
}

let ticking = false;

function updateParallax() {
    // Skip if on mobile/narrow screens to save battery/performance
    if (windowWidth < 768) return;

    // Apply subtle parallax only to stacked cards
    // Active Front card stays stable
    
    terminals.forEach((term) => {
        if (term.classList.contains('active-front')) {
            // Static transform for active card to avoid constant repainting
            // term.style.transform = `translate(-50%, -50%) scale(1)`; 
            // Better: Set it once in CSS class, remove inline style here
             term.style.transform = ''; 
        } else if (term.classList.contains('stacked-right')) {
            const moveX = mouseX * 5; 
            const moveY = mouseY * 5;
            term.style.transform = `translate(calc(-50% + 40px + ${moveX}px), calc(-50% + 15px + ${moveY}px)) scale(0.96)`;
        } else if (term.classList.contains('stacked-left')) {
            const moveX = mouseX * 8;
            const moveY = mouseY * 8;
            term.style.transform = `translate(calc(-50% - 40px + ${moveX}px), calc(-50% + 30px + ${moveY}px)) scale(0.92)`;
        }
    });
}

// Removed fanOutCards function

function collapseAndShowFocus() {
    // Design: 
    // Main (4) -> Center Front
    // Logs (5) -> Stacked Right (Just behind)
    // Plot (2) -> Stacked Left (Just behind)
    // Help (3) -> Hidden or deep stack
    // Analysis (1) -> Hidden or deep stack
    // Misc (0) -> Hidden or deep stack

    // Reset styles
    terminals.forEach(t => {
        t.className = 'terminal'; // Reset classes
        t.style = ''; // Reset inline styles
        // Ensure opacity is 0 initially so they don't flash, but CSS handles fade in?
        // Actually, we want them visible now.
    });

    // Main Terminal (Center)
    const mainTerm = terminals[4];
    mainTerm.classList.add('active-front');
    
    // Stacked Right Group (Logs)
    terminals[5].classList.add('stacked-right');
    
    // Stacked Left Group (Plot)
    terminals[2].classList.add('stacked-left');
    
    // Others can be hidden or just put in back (stacked-left with lower z-index)
    [terminals[0], terminals[1], terminals[3]].forEach(t => {
        t.classList.add('stacked-left');
        t.style.zIndex = 20; // Deep back
        t.style.opacity = '0'; // Hide deep ones for cleaner look
    });

    // Render Content Immediately
    typeMainTerminal();
}

function typeMainTerminal() {
    const contentEl = document.getElementById('content-main');
    if (!contentEl) return;
    
    // Reset content
    contentEl.innerHTML = '';
    contentEl.style.whiteSpace = 'pre'; // Preserve formatting
    
    // Instant render instead of typing
    const cmdDiv = document.createElement('div');
    cmdDiv.innerHTML = `<span style="color:#27ae60; font-weight:600;">user@gpumdkit:~$</span> <span style="font-weight:500;">$ gpumdkit.sh</span>`;
    contentEl.appendChild(cmdDiv);
    
    const textNode = document.createTextNode(mainContentText);
    contentEl.appendChild(textNode);
}

function typeLogsTerminal() {
    const contentEl = document.getElementById('content-logs');
    if (!contentEl) return;
    
    // Check if already typed
    if (contentEl.textContent.trim().length > 0) return;

    contentEl.innerHTML = '';
    contentEl.style.whiteSpace = 'pre';
    
    // Command 1
    const cmdDiv1 = document.createElement('div');
    cmdDiv1.innerHTML = `<span style="color:#27ae60; font-weight:600;">user@gpumdkit:~$</span> <span style="font-weight:500;">gpumdkit.sh -out2xyz .</span>`;
    contentEl.appendChild(cmdDiv1);

    const outDiv1 = document.createElement('div');
    outDiv1.textContent = logsCommand1;
    outDiv1.style.color = 'var(--text-color)';
    outDiv1.style.marginBottom = '20px';
    contentEl.appendChild(outDiv1);

    // Command 2
    const cmdDiv2 = document.createElement('div');
    cmdDiv2.innerHTML = `<span style="color:#27ae60; font-weight:600;">user@gpumdkit:~$</span> <span style="font-weight:500;">gpumdkit.sh -replicate unitcell.vasp supercell.vasp 8000</span>`;
    contentEl.appendChild(cmdDiv2);

    const outDiv2 = document.createElement('div');
    outDiv2.textContent = logsCommand2;
    outDiv2.style.color = 'var(--text-color)';
    outDiv2.style.marginBottom = '20px';
    contentEl.appendChild(outDiv2);

    // Prompt 3
    const cmdDiv3 = document.createElement('div');
    cmdDiv3.innerHTML = `<span style="color:#27ae60; font-weight:600;">user@gpumdkit:~$</span> <span class="cursor">_</span>`;
    contentEl.appendChild(cmdDiv3);
}

window.addEventListener('load', init);