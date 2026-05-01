import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import io

csv_data = """tool,db_load_sec,classify_sec,wall_time_sec,total_reads,throughput,peak_rss_kb
TDKC,7.48,194.68,298.31,494039772,2537713,2937828
TDKC-A,21.73,196.50,202.85,494039772,2514186,5195484
TDKC-D,311.17,281.40,592.57,494039772,1755649,40076752
TDKC-AD,324.11,296.20,620.31,494039772,1667936,42313316
Full-Taxon,392,914.533,1307,494039772,540209,98721404
"""

df = pd.read_csv(io.StringIO(csv_data))

df['Speed (Mreads/min)'] = (df['throughput'] * 60) / 1e6

df['Memory usage (GB)'] = df['peak_rss_kb'] / (1024**2)

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'axes.linewidth': 1.0,
    'figure.dpi': 400 
})

fig, ax = plt.subplots(figsize=(6.5, 4))

bar_width = 0.35
x = np.arange(len(df['tool']))

color_speed = '#E81F28'
color_mem = '#2A2676'

bar1 = ax.bar(x - bar_width/2, df['Speed (Mreads/min)'], bar_width, 
              label='Processing speed', color=color_speed, zorder=3)
bar2 = ax.bar(x + bar_width/2, df['Memory usage (GB)'], bar_width, 
              label='Memory usage', color=color_mem, zorder=3)

ax.set_ylabel('Speed (Mreads/min) and\nmemory usage (GB)', linespacing=1.2)
ax.set_xticks(x)
ax.set_xticklabels(df['tool'])

ax.set_ylim(0, max(df['Memory usage (GB)'].max(), df['Speed (Mreads/min)'].max()) * 1.15)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.grid(axis='y', linestyle='-', linewidth=0.5, color='#CCCCCC', zorder=0)

ax.legend(
    loc='upper right',
    bbox_to_anchor=(1, 1.1),
    frameon=True,
    edgecolor='black',
    fancybox=False,
    framealpha=1
)
plt.tight_layout()
plt.savefig('performance_benchmark.png', format='png', bbox_inches='tight')
plt.show()