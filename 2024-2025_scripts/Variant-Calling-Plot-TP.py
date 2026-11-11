# This code is the percentage of TRUE POSITIVE values in the variant calling --> Variants from the Reference Variant Calling that are also found in the MAF TCGA file / Total Number of Variants in the Reference Variant Calling File 

import matplotlib.pyplot as plt
import numpy as np 

donors = ['TCGA-FD-A3N6', 'TCGA-4Z-AA86','TCGA-BT-A20V', 'TCGA-FD-A3B4']
# Values are pangenome then linear proportion 
TCGA_FD_A3N6_values = [30, 31.25]
TCGA_4Z_AA86_values = [9.95, 42.57]
TCGA_BT_A20V_values = [37.65,32.44]
TCGA_FD_A3B4_values = [23.38, 20.34]

data_values = [TCGA_FD_A3N6_values, TCGA_4Z_AA86_values, TCGA_BT_A20V_values, TCGA_FD_A3B4_values]

# Organizes x-axis positions for donors
x = np.arange(len(donors)) * 2

bar_width = 0.5

fig, ax = plt.subplots(figsize = (10,10))

for i, values in enumerate(data_values):
    ax.bar(x[i] - bar_width/2, values[0], bar_width, label=f'{donors[i]} - Pangenome Percentage - {values[0]}%', color = 'lightskyblue')
    ax.bar(x[i] + bar_width/2, values[1], bar_width, label=f'{donors[i]} - Linear Percentage - {values[1]}%', color = 'royalblue')

ax.set_xlabel('Donor ID')
ax.set_ylabel('Percentage of Variants Found in Reference and MAF File')
ax.set_title('Percentage of Variants Found in Each Reference and MAF File')
ax.set_xticks(x)
ax.set_xticklabels(donors)
ax.legend(loc='upper left', bbox_to_anchor=(0,1), fontsize =7)

plt.savefig('TP_Percentage_of_Variants_Found_in_Each_Reference_and_MAF_File.pdf')

plt.show