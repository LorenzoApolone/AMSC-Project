import matplotlib.pyplot as plt
import os
import re
from collections import defaultdict

# Definizione dei test e dei file associati
tests = {
    "test1": {
        "file": "weak_dim64_20260403_120757.txt",
        "title_suffix": "Test 1 (Dim 64, Part 32*NP)"
    },
    "test2": {
        "file": "weak_dim_x_np.txt",
        "title_suffix": "Test 2 (Dim 8*NP, Part 224)"
    },
    "test3": {
        "file": "weak_local_constant_20260406_012538.txt",
        "title_suffix": "Test 3 (Dim 4*NP, Part 32*NP)"
    }
}

style_map = {
    'classic': {'color': '#1f77b4', 'marker': 'o', 'label': 'Classic'},
    'random': {'color': '#ff7f0e', 'marker': 's', 'label': 'Random'},
    'scale_free': {'color': '#2ca02c', 'marker': '^', 'label': 'Scale Free'},
    'small_world': {'color': '#d62728', 'marker': 'D', 'label': 'Small World'}
}
plot_order = ['classic', 'random', 'scale_free', 'small_world']

plt.style.use('default')

for test_name, test_info in tests.items():
    file_path = test_info["file"]
    if not os.path.exists(file_path):
        print(f"File non trovato: {file_path}")
        continue
    
    data = defaultdict(lambda: defaultdict(lambda: {'total_time': [], 'converged': [], 'overhead_pct': [], 'comp_time': []}))
    current_cores = -1
    
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
                topology = parts[1]
                # Prendiamo solo i run paralleli (cores > 0) per i grafici
                if current_cores > 0:
                    t_time = float(parts[7])
                    o_time = float(parts[8])
                    conv = int(parts[9])
                    
                    ovh_pct = (o_time / t_time * 100.0) if t_time > 0 else 0.0
                    comp_time = t_time - o_time
                    
                    data[topology][current_cores]['total_time'].append(t_time)
                    data[topology][current_cores]['converged'].append(conv)
                    data[topology][current_cores]['overhead_pct'].append(ovh_pct)
                    data[topology][current_cores]['comp_time'].append(comp_time)
                    
    # Calcolo delle medie
    mean_data = defaultdict(lambda: {'cores': [], 'time': [], 'converged': [], 'overhead_pct': [], 'comp_time': []})
    all_cores = set()
    
    for top, core_dict in data.items():
        for cores, metrics in sorted(core_dict.items()):
            all_cores.add(cores)
            mean_time = sum(metrics['total_time']) / len(metrics['total_time'])
            mean_conv = sum(metrics['converged']) / len(metrics['converged'])
            mean_ovh = sum(metrics['overhead_pct']) / len(metrics['overhead_pct'])
            mean_comp = sum(metrics['comp_time']) / len(metrics['comp_time'])
            
            mean_data[top]['cores'].append(cores)
            mean_data[top]['time'].append(mean_time)
            mean_data[top]['converged'].append(mean_conv)
            mean_data[top]['overhead_pct'].append(mean_ovh)
            mean_data[top]['comp_time'].append(mean_comp)
            
    if not all_cores:
        continue
        
    sorted_cores = sorted(list(all_cores))
    
    # =========================================================================
    # 1. Plot Execution Time
    # =========================================================================
    plt.figure(figsize=(10, 6))
    for top in plot_order:
        if top in mean_data:
            vals = mean_data[top]
            st = style_map[top]
            plt.plot(vals['cores'], vals['time'], marker=st['marker'], markersize=6, linewidth=1.5, label=st['label'], color=st['color'])
    
    plt.title(f'Execution Time vs Number of Cores\nWeak Scaling: {test_info["title_suffix"]}')
    plt.xlabel('Number of Cores')
    plt.ylabel('Execution Time (s)')
    plt.xticks(sorted_cores)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{test_name}_time.png', dpi=100)
    plt.close()
    
    # =========================================================================
    # 2. Plot Weak Scaling Efficiency
    # =========================================================================
    plt.figure(figsize=(10, 6))
    for top in plot_order:
        if top in mean_data:
            vals = mean_data[top]
            st = style_map[top]
            
            t1 = None
            for c, t in zip(vals['cores'], vals['time']):
                if c == 1:
                    t1 = t
                    break
            
            if t1 is not None:
                eff = [(t1 / t) * 100 for t in vals['time']]
                plt.plot(vals['cores'], eff, marker=st['marker'], markersize=6, linewidth=1.5, label=st['label'], color=st['color'])
                
    plt.axhline(y=100, color='black', linestyle='--', linewidth=1.5, label='Ideal (100%)')
    
    plt.title(f'Weak Scaling Efficiency\n{test_info["title_suffix"]}')
    plt.xlabel('Number of Cores')
    plt.ylabel('Efficiency (%)')
    plt.xticks(sorted_cores)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{test_name}_eff.png', dpi=100)
    plt.close()
    
    # =========================================================================
    # 3. Plot Convergence (Solo Test 2 e Test 3)
    # =========================================================================
    if test_name in ["test2", "test3"]:
        plt.figure(figsize=(10, 6))
        for top in plot_order:
            if top in mean_data:
                vals = mean_data[top]
                st = style_map[top]
                plt.plot(vals['cores'], vals['converged'], marker=st['marker'], markersize=6, linewidth=1.5, label=st['label'], color=st['color'])
        
        plt.title(f'Convergence vs Number of Cores (Increasing Dim)\nWeak Scaling: {test_info["title_suffix"]}')
        plt.xlabel('Number of Cores (Proxy for Problem Dimension)')
        plt.ylabel('Average Converged Functions (out of 45)')
        plt.xticks(sorted_cores)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{test_name}_conv.png', dpi=100)
        plt.close()

    # =========================================================================
    # 4. Plot Overhead Percentage
    # =========================================================================
    plt.figure(figsize=(10, 6))
    for top in plot_order:
        if top in mean_data:
            vals = mean_data[top]
            st = style_map[top]
            plt.plot(vals['cores'], vals['overhead_pct'], marker=st['marker'], markersize=6, linewidth=1.5, label=st['label'], color=st['color'])
    
    plt.title(f'Communication Overhead Percentage vs Number of Cores\nWeak Scaling: {test_info["title_suffix"]}')
    plt.xlabel('Number of Cores')
    plt.ylabel('Overhead Percentage (%)')
    plt.xticks(sorted_cores)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{test_name}_ovh.png', dpi=100)
    plt.close()

    # =========================================================================
    # 5. Plot Computation Time (Total - Overhead)
    # =========================================================================
    plt.figure(figsize=(10, 6))
    for top in plot_order:
        if top in mean_data:
            vals = mean_data[top]
            st = style_map[top]
            plt.plot(vals['cores'], vals['comp_time'], marker=st['marker'], markersize=6, linewidth=1.5, label=st['label'], color=st['color'])
    
    plt.title(f'Computation Time (Without Overhead) vs Number of Cores\nWeak Scaling: {test_info["title_suffix"]}')
    plt.xlabel('Number of Cores')
    plt.ylabel('Computation Time (s)')
    plt.xticks(sorted_cores)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{test_name}_comp.png', dpi=100)
    plt.close()

print("Generazione dei grafici di Weak Scaling completata con successo, inclusi overhead e comp_time!")
