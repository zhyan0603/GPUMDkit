import sys
from ase.io import read, write

input_file = sys.argv[1] if len(sys.argv) > 1 else 'train.xyz'

# Read all frames
frames = read(input_file, index=':')

# Save to POSCAR
for i, frame in enumerate(frames):
    poscar_filename = f'POSCAR_{i + 1}.vasp'
    write(poscar_filename, frame)

print(f'All frames have been converted.')