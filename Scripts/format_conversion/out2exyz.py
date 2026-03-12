import os, sys
import numpy as np
from ase.io import read, write
from ase import Atoms, Atom
from tqdm import tqdm

def Convert_atoms(atom):
    xx,yy,zz,yz,xz,xy = -atom.calc.results['stress']*atom.get_volume() 
    atom.info['virial'] = np.array([(xx, xy, xz), (xy, yy, yz), (xz, yz, zz)])
    atom.calc.results['energy'] = atom.calc.results['free_energy']
    del atom.calc.results['stress']
    del atom.calc.results['free_energy']

def find_outcar(start_path='.'):
    result = []
    for root, dirs, files in os.walk(start_path):
        if 'OUTCAR' in files:
            result.append(os.path.join(root, 'OUTCAR'))
    return result

def is_converged(outcar_path):
    with open(outcar_path, 'r') as f:
        lines = f.readlines()
    
    nsw = None
    nelm = None
    actual_steps = 0
    has_ediff = False
    
    for line in lines:
        stripped_line = line.strip()
        
        if 'number of steps for IOM' in line:
            # Typical: "   NSW    =      0    number of steps for IOM"
            parts = line.split('=')
            if len(parts) > 1:
                nsw_str = parts[1].split()[0]
                try:
                    nsw = int(nsw_str)
                except ValueError:
                    pass
        
        if 'NELM' in line and 'of ELM steps' in line:
            # Typical: "   NELM   =     60;   NELMDL=      0     of ELM steps"
            parts = line.split('=')
            if len(parts) > 1:
                nelm_str = parts[1].split(';')[0].strip()
                try:
                    nelm = int(nelm_str)
                except ValueError:
                    pass
        
        if 'Iteration' in stripped_line:
            actual_steps += 1
        
        if 'aborting loop because EDIFF is reached' in line:
            has_ediff = True
    
    if nsw is None or nelm is None:
        return False  # Couldn't parse essential info, treat as non-converged
    
    if nsw != 0:
        return True  # Following the shell script logic: assume converged if NSW != 0
    
    # For NSW == 0
    if has_ediff:
        if actual_steps < nelm:
            return True
        else:
            return False
    else:
        return False

file_list = find_outcar(start_path=sys.argv[1])

cnum = 0     # total number of configuration
atoms_list, err_list, non_converged_list = [], [], []
for file_path in tqdm(file_list):
    if not is_converged(file_path):
        non_converged_list.append(file_path)
        continue
    
    try:
        atoms = read(file_path, format='vasp-out', index=":")
    except:
        err_list.append(file_path)
        continue
    
    for ai in range(len(atoms)):
        Convert_atoms(atoms[ai])
        atoms_list.append(atoms[ai])
    cnum += len(atoms)

write('train.xyz', atoms_list, format='extxyz')
print(' The total number of configurations is: {} \n'.format(cnum))

if err_list:
    print(" The list of failed calculation files is as follows.")
    for err_filename in err_list:
        print(err_filename)
else:
    print(" All reads are successful!")

if non_converged_list:
    print(" The list of non-converged calculation files is as follows.")
    for non_conv_filename in non_converged_list:
        print(non_conv_filename)
else:
    print(" All calculations are converged!")