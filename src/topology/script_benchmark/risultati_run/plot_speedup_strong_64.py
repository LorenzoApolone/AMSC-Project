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
                if parts[1] == 'not_converged' or len(parts) < 8:
                    continue
                
                topology = parts[1]
                try:
                    time = float(parts[7])
                    if current_np is not None:
                        data[current_np][topology].append(time)
                except ValueError:
                    pass
    return data

def main():
    filepath = 'strong_dim64_20260403_120756.txt'
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
        avg_times = []
        plot_nps = []
        for np_val in nps:
            if top in data[np_val]:
                times = data[np_val][top]
                avg_times.append(sum(times) / len(times))
                plot_nps.append(np_val)
                
        if len(plot_nps) == 0:
            continue
            
        if 0 in plot_nps:
            base_idx = plot_nps.index(0)
        elif 1 in plot_nps:
            base_idx = plot_nps.index(1)
        else:
            continue
            
        base_time = avg_times[base_idx]
        speedups = [base_time / t for t in avg_times]
        
        plt.plot(plot_nps, speedups, marker=markers[i % len(markers)], label=top.replace('_', ' ').title())
        
    plt.xlabel('Numero di Core (0 = Seriale)')
    plt.ylabel('Speedup (T_base / T_np)')
    plt.title('Speedup Benchmark Strong Scaling (Dim=64, Particelle=256)')
    plt.grid(True)
    plt.legend()
    
    output_png = 'strong_dim64_speedup.png'
    plt.savefig(output_png, dpi=300)
    print(f"Plot generato e salvato in: {output_png}")

if __name__ == '__main__':
    main()
