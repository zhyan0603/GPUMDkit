import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Configuration for custom file paths and legend labels (edit as needed)
FILE_CONFIG = [
    # {'path': 'des_Li_1978.npy', 'label': r'NEP$_{std}$'},
    # {'path': 'des_Li_802.npy', 'label': r'NEP$_{802}$'},
    # {'path': 'des_Li_25.npy', 'label': r'NEP$_{25}$'},
]

def load_descriptors(file_paths):
    """Load .npy files and return descriptors and labels."""
    descriptors = []
    labels = []
    for path in file_paths:
        if not os.path.isfile(path):
            print(f"Error: File '{path}' not found.")
            sys.exit(1)
        print(f"Loading '{path}'...")
        data = np.load(path)
        descriptors.append(data)
        # Use custom label if provided, else use file name without .npy
        label = next((item['label'] for item in FILE_CONFIG if item['path'] == path), os.path.splitext(os.path.basename(path))[0])
        labels.append(label)
    return descriptors, labels

def reduce_dimensionality(descriptors, method='pca'):
    """Apply PCA or UMAP to reduce descriptors to 2D."""
    # Concatenate all descriptors
    combined = np.concatenate(descriptors, axis=0)
    
    print(f"Applying {method.upper()} reduction...")
    if method.lower() == 'pca':
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        # Standardize data for PCA
        scaler = StandardScaler()
        combined_scaled = scaler.fit_transform(combined)
        reducer = PCA(n_components=2)
        reduced = reducer.fit_transform(combined_scaled)
    elif method.lower() == 'umap':
        import umap
        reducer = umap.UMAP(n_components=2, random_state=42)
        reduced = reducer.fit_transform(combined)
    else:
        print("Error: Method must be 'pca' or 'umap'.")
        sys.exit(1)
    
    # Split reduced data back into original groups
    reduced_splits = []
    start = 0
    for data in descriptors:
        end = start + data.shape[0]
        reduced_splits.append(reduced[start:end])
        start = end
    return reduced_splits

def plot_reduced_data(reduced_splits, labels, method):
    """Create and save a 2D scatter plot of reduced data."""
    plt.figure(figsize=(6, 4.5), dpi=100)
    #colors = plt.cm.tab10(np.linspace(0, 1, len(reduced_splits)))
    colors = ['#9BBBE1', '#EAB883', '#A9CA70', '#DD7C4F', '#F09BA0', '#B58C9A'] 
    
    for reduced, label, color in zip(reduced_splits, labels, colors):
        plt.scatter(reduced[:, 0], reduced[:, 1], label=label, color=color, alpha=0.6, s=15)
    
    #plt.title(f'{method.upper()} Analysis')
    plt.xlabel(f'{method.upper()} Component 1')
    plt.ylabel(f'{method.upper()} Component 2')
    plt.legend()
    plt.tight_layout()
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'descriptors.png'.")
        plt.savefig('descriptors.png', dpi=300)
    else:
        plt.show()

def main():
    # Check command-line arguments
    if len(sys.argv) < 3:
        print("Usage: python scripts.py <method> <file1.npy> <file2.npy> ...")
        print("Method: 'pca' or 'umap'")
        sys.exit(1)
    
    method = sys.argv[1]
    file_paths = sys.argv[2:]
    
    # Use FILE_CONFIG if provided, else use command-line files
    if FILE_CONFIG:
        file_paths = [item['path'] for item in FILE_CONFIG]
    
    # Load descriptors
    descriptors, labels = load_descriptors(file_paths)
    
    # Reduce dimensionality
    reduced_splits = reduce_dimensionality(descriptors, method)
    
    # Plot and save
    plot_reduced_data(reduced_splits, labels, method)

if __name__ == "__main__":
    main()