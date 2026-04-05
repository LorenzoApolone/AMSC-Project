import sys
import os
import re
import matplotlib.pyplot as plt
from collections import defaultdict

filepath = '/home/lorenzo/AMSC/AMSC-Project/src/topology/script_benchmark/risultati_run/weak_dim64_20260403_120757.txt'
output_png = '/home/lorenzo/AMSC/AMSC-Project/src/topology/script_benchmark/risultati_run/plot_weak_dim64.png'

data = defaultdict(lambda: defaultdict(list))
current_np = None

if not os.path.exists(filepath):
    print('File not found')
    sys.exit(1)

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
    plt.plot(plot_nps, avg_times, marker=markers[i % len(markers)], label=top.replace('_', ' ').title())

plt.xlabel('Numero di Core (0 = Seriale)')
plt.ylabel('Tempo Medio (s)')
plt.title('Benchmark Weak Scaling (Dim=64, Part/Proc=32)')
plt.grid(True)
plt.legend()
plt.savefig(output_png, dpi=300)
print('PLOT_COMPLETATO: ' + output_png)