import sys
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('loss.out')

epochs = np.arange(1, len(data) + 1) 
learning_rates = data[:, -2]  

print(f"Ploting the learning rate during GNEP training process...")

plt.figure(figsize=(5, 3.5))
plt.plot(epochs, learning_rates, c='C0', label='Learning Rate')
plt.xlabel('Epoch')
plt.ylabel('Learning Rate')
plt.legend()
plt.tight_layout()

if len(sys.argv) > 1 and sys.argv[1] == 'save':
    plt.savefig('learning_rate.png', dpi=300)
else:
    # Handle saving or displaying the plot
    from matplotlib import get_backend
    if get_backend().lower() in ['agg', 'cairo', 'pdf', 'ps', 'svg']:
        print("Unable to display the plot due to the non-interactive backend.")
        print("The plot has been automatically saved as 'learning_rate.png'.")
        plt.savefig('learning_rate.png', dpi=300)
    else:
        plt.show()