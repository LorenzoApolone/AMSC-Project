import matplotlib.pyplot as plt
import os
import re
from collections import defaultdict

files = ["strong_dim32_v1.txt", "strong_dim64_v1.txt", "strong_dim128_v1.txt"]

for file_path in files:
    if not os.path.exists(file_path):
        print(f"Attenzione: {file_path} non trovato.")
        continue
    
    data = defaultdict(lambda: defaultdict(list))
    current_cores = 0
    n_points = None
    max_iter = None
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('=> ESECUZIONE COMANDO:'):
                if 'topology_serial' in line:
                    current_cores = 0
                else:
                    match = re.search(r'-np\s+(\d+)', line)
                    if match:
                        current_cores = int(match.group(1))
            elif line.startswith('RESULT'):
                parts = line.split(',')
                if parts[1] == 'not_converged':
                    continue
                
                # Assicuriamoci di estrarre numero di particelle e max iterazioni dalla prima riga utile
                if n_points is None:
                    n_points = parts[3]
                    max_iter = parts[4]
                
                topology = parts[1]
                if topology == 'classic':
                    continue
                
                if current_cores > 0:
                    total_time = float(parts[7])
                    overhead_time = float(parts[8])
                    if total_time > 0:
                        overhead_percent = (overhead_time / total_time) * 100.0
                    else:
                        overhead_percent = 0.0
                    data[topology][current_cores].append(overhead_percent)
    
    mean_data = defaultdict(lambda: {'cores': [], 'overhead': []})
    all_cores = set()
    
    for top, core_dict in data.items():
        for cores, vals in sorted(core_dict.items()):
            all_cores.add(cores)
            avg_ovh = sum(vals) / len(vals)
            mean_data[top]['cores'].append(cores)
            mean_data[top]['overhead'].append(avg_ovh)
            
    if not all_cores:
        print(f"Nessun dato valido in {file_path}")
        continue
        
    sorted_cores = sorted(list(all_cores))
    
    plt.style.use('default')
    
    style_map = {
        'random': {'color': '#ff7f0e', 'marker': 's', 'label': 'Random'},
        'scale_free': {'color': '#2ca02c', 'marker': '^', 'label': 'Scale Free'},
        'small_world': {'color': '#d62728', 'marker': 'D', 'label': 'Small World'}
    }
    
    plt.figure(figsize=(10, 6))
    
    plot_order = ['random', 'scale_free', 'small_world']
    
    for top in plot_order:
        if top in mean_data:
            vals = mean_data[top]
            st = style_map[top]
            plt.plot(vals['cores'], vals['overhead'], marker=st['marker'], markersize=6, linewidth=1.5, label=st['label'], color=st['color'])
    
    dim_match = re.search(r'dim(\d+)', file_path)
    dim_str = dim_match.group(1) if dim_match else "Unknown"
    
    plt.title(f'Communication Overhead Percentage vs Number of Cores\n(Strong Scaling, Dim={dim_str}, Part={n_points}, MaxIter={max_iter})')
    plt.xlabel('Number of Cores')
    plt.ylabel('Average Overhead Percentage (%)')
    plt.xticks(sorted_cores)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    
    out_name = f'overhead_dim{dim_str}.png'
    plt.savefig(out_name, dpi=100)
    plt.close()
    
    print(f"Generato plot: {out_name}")
