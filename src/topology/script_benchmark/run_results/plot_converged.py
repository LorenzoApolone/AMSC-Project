import matplotlib.pyplot as plt
import os
from collections import defaultdict

files = ["strong_dim32_v1.txt", "strong_dim64_v1.txt", "strong_dim128_v1.txt"]

# data[topology][dim] = list of converged values
data = defaultdict(lambda: defaultdict(list))

for file_path in files:
    if not os.path.exists(file_path):
        print(f"Attenzione: {file_path} non trovato.")
        continue
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('RESULT'):
                parts = line.strip().split(',')
                if parts[1] == 'not_converged':
                    continue
                
                topology = parts[1]
                dim = int(parts[2])
                converged = int(parts[9])
                
                data[topology][dim].append(converged)

# Compute means
mean_data = defaultdict(lambda: {'dim': [], 'converged': []})
all_dims = set()

for top, dim_dict in data.items():
    for d, vals in sorted(dim_dict.items()):
        all_dims.add(d)
        avg_conv = sum(vals) / len(vals)
        mean_data[top]['dim'].append(d)
        mean_data[top]['converged'].append(avg_conv)

# Plotting
plt.style.use('default')

style_map = {
    'classic': {'color': '#1f77b4', 'marker': 'o', 'label': 'Classic'},
    'random': {'color': '#ff7f0e', 'marker': 's', 'label': 'Random'},
    'scale_free': {'color': '#2ca02c', 'marker': '^', 'label': 'Scale Free'},
    'small_world': {'color': '#d62728', 'marker': 'D', 'label': 'Small World'}
}

sorted_dims = sorted(list(all_dims))

plt.figure(figsize=(10, 6))

# Define an order so the legend matches the image
plot_order = ['classic', 'random', 'scale_free', 'small_world']

for top in plot_order:
    if top in mean_data:
        vals = mean_data[top]
        st = style_map[top]
        plt.plot(vals['dim'], vals['converged'], marker=st['marker'], markersize=6, linewidth=1.5, label=st['label'], color=st['color'])

plt.title('Converged Functions vs Problem Dimension (Strong Scaling)')
plt.xlabel('Problem Dimension')
plt.ylabel('Average Number of Converged Functions (out of 45)')
plt.xticks(sorted_dims)
plt.grid(True)

# Legend matching the style exactly (no title, standard frame)
plt.legend()
plt.tight_layout()
plt.savefig('converged_vs_dim.png', dpi=100)
plt.close()

print("Plot generato con successo: converged_vs_dim.png")
