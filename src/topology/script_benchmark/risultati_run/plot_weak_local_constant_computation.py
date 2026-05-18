import sys
import os
import re
import matplotlib.pyplot as plt
from collections import defaultdict

def parse_file(filepath):
    data = defaultdict(lambda: defaultdict(list))
    current_np = None
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('=> ESECUZIONE COMANDO:'):
                if './topology_serial' in line:
                    current_np = 0
                else:
                    m = re.search(r'-np (\d+)', line)
                    if m:
                        current_np = int(m.group(1))
            elif line.startswith('RESULT,'):
                parts = line.split(',')
                if parts[1] == 'not_converged' or len(parts) < 11:
                    continue
                
                topology = parts[1]
                try:
                    # Index -4 is total time, index -3 is communication overhead
                    total_time = float(parts[-4])
                    comm_time = float(parts[-3])
                    comp_time = total_time - comm_time
                    
                    if current_np is not None and current_np > 0: # Only consider parallel runs (cores >= 1)
                        data[current_np][topology].append(comp_time)
                except ValueError:
                    pass
    return data

def main():
    filepath = 'weak_local_constant_20260406_012538.txt'
    if not os.path.exists(filepath):
        print(f"Errore: {filepath} non trovato.")
        sys.exit(1)
        
    data = parse_file(filepath)
    
    topologies = set()
    for np_val in data:
        for top in data[np_val]:
            topologies.add(top)
            
    nps = sorted(list(data.keys()))
    
    plt.figure(figsize=(10, 6))
    markers = ['o', 's', '^', 'D', 'v', 'p']
    
    for i, top in enumerate(sorted(topologies)):
        avg_comp_times = []
        plot_nps = []
        for np_val in nps:
            if top in data[np_val]:
                times = data[np_val][top]
                avg_comp_times.append(sum(times) / len(times))
                plot_nps.append(np_val)
        
        plt.plot(plot_nps, avg_comp_times, marker=markers[i % len(markers)], label=top.replace('_', ' ').title())
        
    plt.xlabel('Cores')
    plt.ylabel('Computation Time (Total - Overhead) (s)')
    plt.title('Computation Time vs Cores (Weak Scaling Local Constant)')
    plt.grid(True)
    plt.legend()
    
    output_png = 'weak_local_constant_computation_plot.png'
    plt.savefig(output_png, dpi=300)
    print(f"Plot generated and saved to: {output_png}")

if __name__ == '__main__':
    main()
