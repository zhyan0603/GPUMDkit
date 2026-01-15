import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def parse_xyz_file(xyz_file):
    """Manually parse xyz file to extract force information"""
    try:
        with open(xyz_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        frames = []
        num_atoms_list = []
        i = 0
        while i < len(lines):
            # Read number of atoms
            try:
                num_atoms = int(lines[i].strip())
            except:
                i += 1
                continue
            
            # Read comment line, check for force property
            comment_line = lines[i+1].strip()
            
            # Check force property name
            force_property = None
            if "forces" in comment_line:
                force_property = "forces"
            elif "force" in comment_line:
                force_property = "force"
            
            # If no force property found, skip this frame
            if not force_property:
                i += num_atoms + 2
                continue
            
            # Extract force data
            forces = []
            for j in range(i+2, i+2+num_atoms):
                parts = lines[j].strip().split()
                # Forces are usually after position coordinates
                if len(parts) >= 6:  # element + 3 positions + 3 forces
                    try:
                        fx, fy, fz = float(parts[-3]), float(parts[-2]), float(parts[-1])
                        forces.append([fx, fy, fz])
                    except ValueError:
                        # If conversion fails, skip this atom
                        continue
            
            if forces and len(forces) == num_atoms:
                frames.append(np.array(forces))
                num_atoms_list.append(num_atoms)
            
            i += num_atoms + 2
        
        return frames, num_atoms_list
    
    except Exception as e:
        print(f" Error parsing {xyz_file}: {e}")
        return None, None

def calculate_net_forces(xyz_file):
    """Calculate net forces using manual parsing"""
    frames, num_atoms_list = parse_xyz_file(xyz_file)
    
    if frames is None or len(frames) == 0:
        return None, 0
    
    total_frames = len(frames)
    print(f" Total configurations: {total_frames}")
    
    net_forces = []
    for idx, forces in enumerate(frames):
        num_atoms = num_atoms_list[idx]
        sum_forces = np.sum(forces, axis=0)
        net_force_magnitude = np.linalg.norm(sum_forces)
        avg_net_force = net_force_magnitude / num_atoms  # Average over number of atoms
        net_forces.append(avg_net_force)
    
    return np.array(net_forces), total_frames

def plot_net_force_distribution(net_forces, output_base, save_flag):
    """
    Plot net force distribution histogram
    """
    # Simple style settings
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 10,
        'axes.labelsize': 11,
        'axes.titlesize': 11,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'lines.linewidth': 1.5,
        'axes.linewidth': 1.2,
        'figure.figsize': (4.5, 3.2)
    })
    
    fig, ax = plt.subplots()
    
    # Filter out zero values
    forces_nonzero = net_forces[net_forces > 1e-10]
    
    if len(forces_nonzero) == 0:
        print(f" Warning: All net forces are zero or near zero.")
        plt.close()
        return
    
    # Increase bin count for denser histogram
    log_bins = np.logspace(np.log10(forces_nonzero.min()), 
                           np.log10(forces_nonzero.max()), 
                           100)  # 100 bins for smoother look
    
    # Plot histogram (density normalized)
    ax.hist(forces_nonzero, bins=log_bins, 
            density=True, alpha=0.8,
            color='dodgerblue', edgecolor='black', 
            linewidth=0.5)
    
    # Set axes
    ax.set_xscale('log')
    ax.set_xlabel(r'Net Force (eV/$\mathrm{{\AA}}$)')
    ax.set_ylabel('Density')
    
    # Add grid for better readability
    ax.grid(True, which='major', linestyle='--', alpha=0.5)
    
    # Set axis limits
    ax.set_xlim(forces_nonzero.min() * 0.5, forces_nonzero.max() * 2)
    
    plt.title('Net Force Distribution')
    plt.tight_layout()
    
    if save_flag:
        # Save as both PDF and PNG
        plt.savefig(f"{output_base}.png", dpi=300)
        print(f" Plot saved as {output_base}.png")
    else:
        # Check if the current backend is non-interactive
        from matplotlib import get_backend
        if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
            print("Unable to display the plot due to the non-interactive backend.")
            print(f"The plot has been automatically saved as {output_base}.png")
            plt.savefig(f"{output_base}.png", dpi=300, bbox_inches='tight')
        else:
            plt.show()
    
    plt.close()  # Close figure to release memory

def check_convergence(xyz_file, threshold=0.001, save_flag=False):
    """
    Check convergence for a single train.xyz file
    """
    if not os.path.exists(xyz_file):
        print(f" File not found: {xyz_file}")
        return
    
    print(" "+"-"*40)
    print(f" Checking file: {xyz_file}")
    net_forces, total_frames = calculate_net_forces(xyz_file)
    
    if net_forces is not None and total_frames > 0:
        # Calculate number of unconverged configurations
        unconverged = np.sum(net_forces > threshold)
        unconverged_ratio = unconverged / total_frames * 100
        median_force = np.median(net_forces)
        
        # Plot distribution
        folder_path = os.path.dirname(xyz_file) or '.'
        output_base = os.path.join(folder_path, 'net_force_distribution')
        plot_net_force_distribution(net_forces, output_base, save_flag)
        
        # Save statistics
        stats_file = 'net_force_stats.txt'
        with open(stats_file, 'w', encoding='utf-8') as f:
            f.write(f"Net force statistics - {xyz_file}\n")
            f.write("="*50 + "\n")
            f.write(f"Total frames: {total_frames}\n")
            f.write(f"Suspicious frames (>{threshold} eV/Angstrom): {unconverged}\n")
            f.write(f"Suspicious ratio: {unconverged_ratio:.2f}%\n")
            f.write(f"Median net force: {median_force:.6f} eV/Angstrom\n")
            f.write(f"Max net force: {np.max(net_forces):.6f} eV/Angstrom\n")
            f.write(f"Min net force: {np.min(net_forces):.6f} eV/Angstrom\n")
        
        print(f" Suspicious ratio: {unconverged_ratio:.2f}%")
        print(f" Median net force: {median_force:.6f} eV/Angstrom")
        print(f" Results saved to: {stats_file}")
        print(" "+"-"*40)
    else:
        print(f" Cannot process {xyz_file} or no force data found")

# Main function
def main():
    if len(sys.argv) < 2:
        print(" Usage: python plt_net_force.py train.xyz [save]")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    save_flag = len(sys.argv) > 2 and sys.argv[2].lower() == 'save'
    
    # Check convergence with threshold 0.001 eV/Angstrom
    check_convergence(xyz_file, threshold=0.001, save_flag=save_flag)

if __name__ == "__main__":
    main()