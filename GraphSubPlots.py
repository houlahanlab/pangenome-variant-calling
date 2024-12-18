import matplotlib.pyplot as plt
import re
import numpy as np  

# This function will read the txt stat files and extract the headings sent for each specific file (the total number of sequences, the properly paired, and total mapped, and mapping quality values)
def extract_data_from_files(filenames, line_headings):
    extracted_data = []
    for filename in filenames: 
        data = {}
        with open(filename, 'r') as file:
            for line in file: 
                normalized_line = line.strip().lower()
                for heading in line_headings:
                    if heading.lower() in normalized_line:
                            number = re.findall(r'\d+', line)
                            if number: 
                                data[heading] = int(number[0])
        extracted_data.append(data)
    return extracted_data

# This function will add bar value numbers to each bar and centers these values 
def add_value_labels(ax, bars, bar_values, font_size):
     for i, bar in enumerate(bars):
          yval = bar_values[i]
          ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(), f"{yval:.2f}", ha='center', va='bottom' if yval >=0 else 'top', fontsize = font_size)

# These are lists with the files of each ancestry and donor (the order is blood than bladder stats file for each donor)
linear_reference_aligned_files = ['AA_Patients/TCGA-FD-A5BU/2024_10_31_TCGA-FD-A5BU-10A-01D-A26K-08_linear_aligned_stats.txt','AA_Patients/TCGA-FD-A5BU/2024_10_31_TCGA-FD-A5BU-01A-31D-A26M-08_linear_aligned_stats.txt','AA_Patients/TCGA-FD-A3N6/GRCh38_linear/2024_10_25_TCGA-FD-A3N6-10A-01D-A21A-08_linear_aligned_stats.txt','AA_Patients/TCGA-FD-A3N6/GRCh38_linear/2024_10_31_TCGA-FD-A3N6-01A-11D-A21A-08_linear_aligned_stats.txt','AA_Patients/TCGA-FD-A3B7/2024_10_31_TCGA-FD-A3B7-10A-01D-A20D-08_linear_aligned_stats.txt','AA_Patients/TCGA-FD-A3B7/2024_10_31_TCGA-FD-A3B7-01A-31D-A20D-08_linear_aligned_stats.txt','AA_Patients/TCGA-FD-A43P/2024_10_31_TCGA-FD-A43P-10A-01D-A23U-08_linear_aligned_stats.txt','AA_Patients/TCGA-FD-A43P/2024_10_31_TCGA-FD-A43P-01A-31D-A23U-08_linear_aligned_stats.txt','AA_Patients/TCGA-FD-A3N5/2024_10_31_TCGA-FD-A3N5-10A-01D-A21A-08_linear_aligned_stats.txt', 'AA_Patients/TCGA-FD-A3N5/2024_10_31_TCGA-FD-A3N5-01A-11D-A21A-08_linear_aligned_stats.txt','AA_Patients/TCGA-FD-A3B5/2024_10_31_TCGA-FD-A3B5-10A-01D-A20D-08_linear_aligned_stats.txt','AA_Patients/TCGA-FD-A3B5/2024_10_31_TCGA-FD-A3B5-01A-11D-A20D-08_linear_aligned_stats.txt','AA_Patients/TCGA-FD-A3B4/2024_10_22_TCGA-FD-A3B4-10A-01D-A202-08_linear_aligned_stats.txt', 'AA_Patients/TCGA-FD-A3B4/2024_10_22_TCGA-FD-A3B4-01A-12D-A202-08_linear_aligned_stats.txt','AA_Patients/TCGA-BT-A20V/2024_10_08_TCGA-BT-A20V-11A-11D-A14W-08_linear_aligned_stats.txt', 'AA_Patients/TCGA-BT-A20V/2024_10_08_TCGA-BT-A20V-01A-11D-A14W-08_linear_aligned_stats.txt','AA_Patients/TCGA-4Z-AA7R/2024_10_11_TCGA-4Z-AA7R-10A-01D-A394-08_linear_aligned_stats.txt','AA_Patients/TCGA-4Z-AA7R/2024_10_11_TCGA-4Z-AA7R-01A-11D-A391-08_linear_aligned_stats.txt','AA_Patients/TCGA-4Z-AA86/2024_10_11_TCGA-4Z-AA86-10A-01D-A394-08_linear_aligned_stats.txt','AA_Patients/TCGA-4Z-AA86/2024_10_11_TCGA-4Z-AA86-01A-11D-A391-08_linear_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7Y/2024_10_30_TCGA-4Z-AA7Y-10A-01D-A394-08_linear_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7Y/2024_10_30_TCGA-4Z-AA7Y-01A-11D-A391-08_linear_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7Q/2024_10_30_TCGA-4Z-AA7Q-10A-01D-A394-08_linear_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7Q/2024_10_30_TCGA-4Z-AA7Q-01A-11D-A391-08_linear_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7M/2024_10_30_TCGA-4Z-AA7M-10A-01D-A394-08_linear_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7M/2024_10_30_TCGA-4Z-AA7M-01A-11D-A391-08_linear_aligned_stats.txt','EA_Patients/TCGA-2F-A9KT/2024_10_28_TCGA-2F-A9KT-10A-01D-A38J-08_linear_aligned_stats.txt','EA_Patients/TCGA-2F-A9KT/2024_10_30_TCGA-2F-A9KT-01A-11D-A38G-08_linear_aligned_stats.txt','EA_Patients/TCGA-2F-A9KW/2024_10_30_TCGA-2F-A9KW-10A-01D-A38J-08_linear_aligned_stats.txt','EA_Patients/TCGA-2F-A9KW/2024_10_30_TCGA-2F-A9KW-01A-11D-A38G-08_linear_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7N/2024_10_30_TCGA-4Z-AA7N-10A-01D-A394-08_linear_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7N/2024_10_30_TCGA-4Z-AA7N-01A-11D-A391-08_linear_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7W/2024_10_30_TCGA-4Z-AA7W-10A-01D-A394-08_linear_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7W/2024_10_30_TCGA-4Z-AA7W-01A-11D-A391-08_linear_aligned_stats.txt','EA_Patients/TCGA-2F-A9KO/2024_10_18_TCGA-2F-A9KO-11A-12D-A38J-08_linear_aligned_stats.txt','EA_Patients/TCGA-2F-A9KO/2024_10_18_TCGA-2F-A9KO-01A-11D-A38G-08_linear_aligned_stats.txt','EA_Patients/TCGA-2F-A9KP/2024_10_18_TCGA-2F-A9KP-10A-01D-A38J-08_linear_aligned_stats.txt','EA_Patients/TCGA-2F-A9KP/2024_10_21_TCGA-2F-A9KP-01A-11D-A38G-08_linear_aligned_stats.txt','EA_Patients/TCGA-2F-A9KQ/2024_10_21_TCGA-2F-A9KQ-11A-11D-A38J-08_aligned_linear_stats.txt','EA_Patients/TCGA-2F-A9KQ/2024_10_22_TCGA-2F-A9KQ-01A-11D-A38G-08_linear_aligned_stats.txt', 'EAA_Patients/TCGA-CF-A1HR/2024_10_22_TCGA-CF-A1HR-10A-01D-A13W-08_linear_aligned_stats.txt', 'EAA_Patients/TCGA-CF-A1HR/2024_10_22_TCGA-CF-A1HR-01A-11D-A13W-08_linear_aligned_stats.txt', 'EAA_Patients/TCGA-CF-A1HS/2024_10_22_TCGA-CF-A1HS-10A-01D-A13W-08_linear_aligned_stats.txt', 'EAA_Patients/TCGA-CF-A1HS/2024_10_22_TCGA-CF-A1HS-01A-11D-A13W-08_linear_aligned_stats.txt', 'EAA_Patients/TCGA-CF-A27C/2024_10_22_TCGA-CF-A27C-10A-01D-A16O-08_linear_aligned_stats.txt', 'EAA_Patients/TCGA-CF-A27C/2024_10_22_TCGA-CF-A27C-01A-11D-A16O-08_linear_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MG/2024_10_30_TCGA-CF-A3MG-10A-01D-A20D-08_linear_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MG/2024_10_30_TCGA-CF-A3MG-01A-11D-A20D-08_linear_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MI/2024_10_25_TCGA-CF-A3MI-10A-01D-A20D-08_linear_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MI/2024_10_29_TCGA-CF-A3MI-01A-11D-A20D-08_linear_aligned_stats.txt','EAA_Patients/TCGA-CF-A47T/2024_10_31_TCGA-CF-A47T-10A-01D-A23U-08_linear_aligned_stats.txt','EAA_Patients/TCGA-CF-A47T/2024_10_31_TCGA-CF-A47T-01A-11D-A23U-08_linear_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MF/2024_10_30_TCGA-CF-A3MF-10A-01D-A21A-08_linear_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MF/2024_10_30_TCGA-CF-A3MF-01A-12D-A21A-08_linear_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MH/2024_10_31_TCGA-CF-A3MH-10A-01D-A20D-08_linear_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MH/2024_10_31_TCGA-CF-A3MH-01A-11D-A20D-08_linear_aligned_stats.txt','EAA_Patients/TCGA-CF-A47S/2024_10_25_TCGA-CF-A47S-10A-01D-A23U-08_linear_aligned_stats.txt','EAA_Patients/TCGA-CF-A47S/2024_10_30_TCGA-CF-A47S-01A-11D-A23U-08_linear_aligned_stats.txt','EAA_Patients/TCGA-CF-A47V/2024_10_31_TCGA-CF-A47V-10A-01D-A23U-08_linear_aligned_stats.txt','EAA_Patients/TCGA-CF-A47V/2024_10_31_TCGA-CF-A47V-01A-11D-A23U-08_linear_aligned_stats.txt']
pangenome_aligned_files = ['AA_Patients/TCGA-FD-A5BU/2024_10_31_TCGA-FD-A5BU-10A-01D-A26K-08_CHM13_aligned_stats.txt','AA_Patients/TCGA-FD-A5BU/2024_10_31_TCGA-FD-A5BU-01A-31D-A26M-08_CHM13_aligned_stats.txt','AA_Patients/TCGA-FD-A3N6/T2T-CHM13_pangenome/2024_11_13_TCGA-FD-A3N6-10A-01D-A21A-08_fragment_length_150_CHM13v1.1_aligned_stats.txt','AA_Patients/TCGA-FD-A3N6/T2T-CHM13_pangenome/2024_11_13_TCGA-FD-A3N6-01A-11D-A21A-08_fragment_length_150_CHM13v1.1_aligned_stats.txt','AA_Patients/TCGA-FD-A3B7/2024_10_31_TCGA-FD-A3B7-10A-01D-A20D-08_CHM13_aligned_stats.txt','AA_Patients/TCGA-FD-A3B7/2024_10_31_TCGA-FD-A3B7-01A-31D-A20D-08_CHM13_aligned_stats.txt','AA_Patients/TCGA-FD-A43P/2024_10_31_TCGA-FD-A43P-10A-01D-A23U-08_CHM13_aligned_stats.txt','AA_Patients/TCGA-FD-A43P/2024_10_31_TCGA-FD-A43P-01A-31D-A23U-08_CHM13_aligned_stats.txt','AA_Patients/TCGA-FD-A3N5/2024_11_13_TCGA-FD-A3N5-10A-01D-A21A-08_fragment_length_150_CHM13v1.1_aligned_stats.txt','AA_Patients/TCGA-FD-A3N5/2024_10_31_TCGA-FD-A3N5-01A-11D-A21A-08_CHM13_aligned_stats.txt','AA_Patients/TCGA-FD-A3B5/2024_10_31_TCGA-FD-A3B5-10A-01D-A20D-08_CHM13_aligned_stats.txt','AA_Patients/TCGA-FD-A3B5/2024_10_31_TCGA-FD-A3B5-01A-11D-A20D-08_CHM13_aligned_stats.txt','AA_Patients/TCGA-FD-A3B4/2024_10_22_TCGA-FD-A3B4-10A-01D-A202-08_CHM13_aligned_stats.txt', 'AA_Patients/TCGA-FD-A3B4/2024_10_22_TCGA-FD-A3B4-01A-12D-A202-08_CHM13_aligned_stats.txt','AA_Patients/TCGA-BT-A20V/2024_10_07_TCGA-BT-A20V-11A-11D-A14W-08_pangenome_aligned_stats.txt', 'AA_Patients/TCGA-BT-A20V/2024_10_08_TCGA-BT-A20V-01A-11D-A14W-08_pangenome_aligned_stats.txt','AA_Patients/TCGA-4Z-AA7R/2024_10_10_TCGA-4Z-AA7R-10A-01D-A394-08_pangenome_aligned_stats.txt','AA_Patients/TCGA-4Z-AA7R/2024_10_11_TCGA-4Z-AA7R-01A-11D-A391-08_pangenome_aligned_stats.txt','AA_Patients/TCGA-4Z-AA86/2024_10_15_TCGA-4Z-AA86-10A-01D-A394-08_pangenome_aligned_stats.txt','AA_Patients/TCGA-4Z-AA86/2024_10_10_TCGA-4Z-AA86-01A-11D-A391-08_pangenome_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7Y/2024_10_30_TCGA-4Z-AA7Y-10A-01D-A394-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7Y/2024_10_30_TCGA-4Z-AA7Y-01A-11D-A391-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7Q/2024_10_30_TCGA-4Z-AA7Q-10A-01D-A394-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7Q/2024_10_30_TCGA-4Z-AA7Q-01A-11D-A391-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7M/2024_10_30_TCGA-4Z-AA7M-10A-01D-A394-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7M/2024_10_30_TCGA-4Z-AA7M-01A-11D-A391-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-2F-A9KT/2024_10_28_TCGA-2F-A9KT-10A-01D-A38J-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-2F-A9KT/2024_10_30_TCGA-2F-A9KT-01A-11D-A38G-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-2F-A9KW/2024_10_30_TCGA-2F-A9KW-10A-01D-A38J-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-2F-A9KW/2024_10_30_TCGA-2F-A9KW-01A-11D-A38G-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7N/2024_10_30_TCGA-4Z-AA7N-10A-01D-A394-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7N/2024_10_30_TCGA-4Z-AA7N-01A-11D-A391-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7W/2024_10_30_TCGA-4Z-AA7W-10A-01D-A394-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-4Z-AA7W/2024_10_30_TCGA-4Z-AA7W-01A-11D-A391-08_CHM13_aligned_stats.txt', 'EA_Patients/TCGA-2F-A9KO/2024_10_18_TCGA-2F-A9KO-11A-12D-A38J-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-2F-A9KO/2024_10_18_TCGA-2F-A9KO-01A-11D-A38G-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-2F-A9KP/2024_10_18_TCGA-2F-A9KP-10A-01D-A38J-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-2F-A9KP/2024_10_21_TCGA-2F-A9KP-01A-11D-A38G-08_CHM13_aligned_stats.txt','EA_Patients/TCGA-2F-A9KQ/2024_10_21_TCGA-2F-A9KQ-11A-11D-A38J-08_CHM13_aligned_stats.txt', 'EA_Patients/TCGA-2F-A9KQ/2024_10_22_TCGA-2F-A9KQ-01A-11D-A38G-08_CHM13_aligned_stats.txt', 'EAA_Patients/TCGA-CF-A1HR/2024_10_22_TCGA-CF-A1HR-10A-01D-A13W-08_CHM13_aligned_stats.txt', 'EAA_Patients/TCGA-CF-A1HR/2024_10_22_TCGA-CF-A1HR-01A-11D-A13W-08_CHM13_aligned_stats.txt', 'EAA_Patients/TCGA-CF-A1HS/2024_10_22_TCGA-CF-A1HS-10A-01D-A13W-08_CHM13_aligned_stats.txt', 'EAA_Patients/TCGA-CF-A1HS/2024_10_22_TCGA-CF-A1HS-01A-11D-A13W-08_CHM13_aligned_stats.txt', 'EAA_Patients/TCGA-CF-A27C/2024_10_22_TCGA-CF-A27C-10A-01D-A16O-08_CHM13_aligned_stats.txt', 'EAA_Patients/TCGA-CF-A27C/2024_10_22_TCGA-CF-A27C-01A-11D-A16O-08_CHM13_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MG/2024_10_30_TCGA-CF-A3MG-10A-01D-A20D-08_CHM13_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MG/2024_10_30_TCGA-CF-A3MG-01A-11D-A20D-08_CHM13_aligned_stats.txt', 'EAA_Patients/TCGA-CF-A3MI/2024_10_25_TCGA-CF-A3MI-10A-01D-A20D-08_CHM13-T2T_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MI/2024_10_29_TCGA-CF-A3MI-01A-11D-A20D-08_CHM13-T2T_aligned_stats.txt','EAA_Patients/TCGA-CF-A47T/2024_10_31_TCGA-CF-A47T-10A-01D-A23U-08_CHM13_aligned_stats.txt','EAA_Patients/TCGA-CF-A47T/2024_10_31_TCGA-CF-A47T-01A-11D-A23U-08_CHM13_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MF/2024_10_30_TCGA-CF-A3MF-10A-01D-A21A-08_CHM13_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MF/2024_10_30_TCGA-CF-A3MF-01A-12D-A21A-08_CHM13_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MH/2024_10_31_TCGA-CF-A3MH-10A-01D-A20D-08_CHM13_aligned_stats.txt','EAA_Patients/TCGA-CF-A3MH/2024_10_31_TCGA-CF-A3MH-01A-11D-A20D-08_CHM13_aligned_stats.txt','EAA_Patients/TCGA-CF-A47S/2024_10_25_TCGA-CF-A47S-10A-01D-A23U-08_CHM13_aligned_stats.txt','EAA_Patients/TCGA-CF-A47S/2024_10_30_TCGA-CF-A47S-01A-11D-A23U-08_CHM13_aligned_stats.txt', 'EAA_Patients/TCGA-CF-A47V/2024_10_31_TCGA-CF-A47V-10A-01D-A23U-08_CHM13_aligned_stats.txt','EAA_Patients/TCGA-CF-A47V/2024_10_31_TCGA-CF-A47V-01A-11D-A23U-08_CHM13_aligned_stats.txt']

# Lists with the specific values being looked for 
linear_reference_aligned_headings = ['raw total sequence', 'reads mapped:', 'reads properly paired']
pangenome_aligned_headings = ['total alignments', 'total aligned', 'total properly paired', 'mapping quality: mean']

linear_reference_aligned_data = extract_data_from_files(linear_reference_aligned_files, linear_reference_aligned_headings)

pangenome_aligned_data = extract_data_from_files(pangenome_aligned_files, pangenome_aligned_headings)

# Calculates the reads aligned for the pangenome and linear alignment for each donor and sample 
reads_aligned_linear = [
    ((data.get('raw total sequence', 0) - data.get('reads mapped:', 0)) /data.get('raw total sequence', 1) if data.get('raw total sequence', 1) > 0 else 0) * 100
    for data in linear_reference_aligned_data
]
print(reads_aligned_linear)

pangenome_aligned = [
    ((data.get('total alignments', 0) - data.get('total aligned', 0)) / data.get('total alignments', 1) if data.get('total alignments', 1) > 0 else 0) * 100
    for data in pangenome_aligned_data
]
print(pangenome_aligned)

bar_values_aligned_linear = [
    data.get('raw total sequence', 0) - data.get('reads mapped:', 0)
    for data in linear_reference_aligned_data
]

bar_values_aligned_pangenome = [
    data.get('total alignments', 0) - data.get('total aligned', 0)
    for data in pangenome_aligned_data
]

bar_width = 3.5
bar_spacing = 5

AA_patient_ID = ['TCGA-FD-A5BU','TCGA-FD-A3N6','TCGA-FD-A3B7','TCGA-FD-A43P','TCGA-FD-A3N5','TCGA-FD-A3B5','TCGA-FD-A3B4','TCGA-BT-A20V', 'TCGA-4Z-AA7R', 'TCGA-4Z-AA86']
EA_patient_ID = ['TCGA-4Z-AA7Y','TCGA-4Z-AA7Q','TCGA-4Z-AA7M','TCGA-2F-A9KT','TCGA-2F-A9KW','TCGA-4Z-AA7N','TCGA-4Z-AA7W','TCGA-2F-A9KO', 'TCGA-2F-A9KP', 'TCGA-2F-A9KQ']
EAA_patient_ID = ['TCGA-CF-A1HR', 'TCGA-CF-A1HS', 'TCGA-CF-A27C','TCGA-CF-A3MG','TCGA-CF-A3MI','TCGA-CF-A47T','TCGA-CF-A3MF','TCGA-CF-A3MH','TCGA-CF-A47S','TCGA-CF-A47V']
total_patient_ID = [AA_patient_ID + EA_patient_ID + EAA_patient_ID]
X_AA = np.arange(len(AA_patient_ID)) * bar_spacing
X_EA = np.arange(len(EA_patient_ID)) * bar_spacing + X_AA[-1] + bar_spacing
X_EAA = np.arange(len(EAA_patient_ID)) * bar_spacing + X_EA[-1] + bar_spacing
total_patient_ID = AA_patient_ID + EA_patient_ID + EAA_patient_ID 

x_total = np.concatenate([X_AA, X_EA, X_EAA])

# These lists carry the information of each donor and find the delta of the unaligned proportion reads by subtracting the unaligned value in the linear reference from the unaligned value in the pangenome reference 
# The tumour data is the odd numbers since they are the second txt files in the linear_reference_aligned_files and pangenome_aligned_files lists 
AA_normal_aligned_data = [pangenome_aligned[0] - reads_aligned_linear[0], pangenome_aligned[2] - reads_aligned_linear[2], pangenome_aligned[4] -reads_aligned_linear[4], pangenome_aligned[6] - reads_aligned_linear[6], pangenome_aligned[8] - reads_aligned_linear[8],pangenome_aligned[10] - reads_aligned_linear[10],pangenome_aligned[12] - reads_aligned_linear[12], pangenome_aligned[14] - reads_aligned_linear[14],pangenome_aligned[16] - reads_aligned_linear[16], pangenome_aligned[18] - reads_aligned_linear[18]]
AA_tumour_aligned_data = [pangenome_aligned[1] - reads_aligned_linear[1], pangenome_aligned[3] - reads_aligned_linear[3], pangenome_aligned[5] - reads_aligned_linear[5], pangenome_aligned[7] - reads_aligned_linear[7], pangenome_aligned[9] - reads_aligned_linear[9],pangenome_aligned[11] - reads_aligned_linear[11],pangenome_aligned[13] - reads_aligned_linear[13],pangenome_aligned[15] - reads_aligned_linear[15], pangenome_aligned[17] - reads_aligned_linear[17],pangenome_aligned[19] - reads_aligned_linear[19]]
EA_normal_aligned_data =[pangenome_aligned[20] - reads_aligned_linear[20], pangenome_aligned[22] - reads_aligned_linear[22], pangenome_aligned[24] - reads_aligned_linear[24], pangenome_aligned[26] - reads_aligned_linear[26], pangenome_aligned[28] - reads_aligned_linear[28],pangenome_aligned[30] - reads_aligned_linear[30],pangenome_aligned[32] - reads_aligned_linear[32],pangenome_aligned[34] - reads_aligned_linear[34],pangenome_aligned[36] - reads_aligned_linear[36],pangenome_aligned[38] - reads_aligned_linear[38]]
EA_tumour_aligned_data =[pangenome_aligned[21] - reads_aligned_linear[21], pangenome_aligned[23] - reads_aligned_linear[23], pangenome_aligned[25] - reads_aligned_linear[25],pangenome_aligned[27] - reads_aligned_linear[27],pangenome_aligned[29] - reads_aligned_linear[29],pangenome_aligned[31] - reads_aligned_linear[31],pangenome_aligned[33] - reads_aligned_linear[33],pangenome_aligned[35] - reads_aligned_linear[35],pangenome_aligned[37] - reads_aligned_linear[37],pangenome_aligned[39] - reads_aligned_linear[39]]
EAA_normal_aligned_data= [pangenome_aligned[40] - reads_aligned_linear[40], pangenome_aligned[42] - reads_aligned_linear[42],pangenome_aligned[44] - reads_aligned_linear[44],pangenome_aligned[46] - reads_aligned_linear[46],pangenome_aligned[48] - reads_aligned_linear[48],pangenome_aligned[50] - reads_aligned_linear[50],pangenome_aligned[52] - reads_aligned_linear[52],pangenome_aligned[54] - reads_aligned_linear[54],pangenome_aligned[56] - reads_aligned_linear[56],pangenome_aligned[58] - reads_aligned_linear[58]]
EAA_tumour_aligned_data = [pangenome_aligned[41] - reads_aligned_linear[41], pangenome_aligned[43] - reads_aligned_linear[43],pangenome_aligned[45] - reads_aligned_linear[45],pangenome_aligned[47] - reads_aligned_linear[47],pangenome_aligned[49] - reads_aligned_linear[49],pangenome_aligned[51] - reads_aligned_linear[51],pangenome_aligned[53] - reads_aligned_linear[53],pangenome_aligned[55] - reads_aligned_linear[55],pangenome_aligned[57] - reads_aligned_linear[57],pangenome_aligned[59] - reads_aligned_linear[59]]   

AA_bar_aligned_values = [bar_values_aligned_pangenome[1] - bar_values_aligned_linear[1], bar_values_aligned_pangenome[3] - bar_values_aligned_linear[3], bar_values_aligned_pangenome[5] - bar_values_aligned_linear[5], bar_values_aligned_pangenome[7] - bar_values_aligned_linear[7], bar_values_aligned_pangenome[9]- bar_values_aligned_linear[9], bar_values_aligned_pangenome[11] - bar_values_aligned_linear[11], bar_values_aligned_pangenome[13]- bar_values_aligned_linear[13], bar_values_aligned_pangenome[15] - bar_values_aligned_linear[15], bar_values_aligned_pangenome[17] - bar_values_aligned_linear[17], bar_values_aligned_pangenome[19] - bar_values_aligned_linear[19], bar_values_aligned_pangenome[0] - bar_values_aligned_linear[0], bar_values_aligned_pangenome[2] - bar_values_aligned_linear[2], bar_values_aligned_pangenome[4] - bar_values_aligned_linear[4], bar_values_aligned_pangenome[6] - bar_values_aligned_linear[6], bar_values_aligned_pangenome[8] - bar_values_aligned_linear[8], bar_values_aligned_pangenome[10] - bar_values_aligned_linear[10], bar_values_aligned_pangenome[12] - bar_values_aligned_linear[12], bar_values_aligned_pangenome[14] - bar_values_aligned_linear[14], bar_values_aligned_pangenome[16] - bar_values_aligned_linear[16], bar_values_aligned_pangenome[18] - bar_values_aligned_linear[18]]
EA_bar_aligned_values = [bar_values_aligned_pangenome[21] - bar_values_aligned_linear[21], bar_values_aligned_pangenome[23] - bar_values_aligned_linear[23], bar_values_aligned_pangenome[25] - bar_values_aligned_linear[25], bar_values_aligned_pangenome[27] - bar_values_aligned_linear[27], bar_values_aligned_pangenome[29] - bar_values_aligned_linear[29], bar_values_aligned_pangenome[31] - bar_values_aligned_linear[31], bar_values_aligned_pangenome[33] - bar_values_aligned_linear[33], bar_values_aligned_pangenome[35] - bar_values_aligned_linear[35], bar_values_aligned_pangenome[37] - bar_values_aligned_linear[37], bar_values_aligned_pangenome[39] - bar_values_aligned_linear[39], bar_values_aligned_pangenome[20] - bar_values_aligned_linear[20], bar_values_aligned_pangenome[22] - bar_values_aligned_linear[22], bar_values_aligned_pangenome[24] - bar_values_aligned_linear[24], bar_values_aligned_pangenome[26] - bar_values_aligned_linear[26], bar_values_aligned_pangenome[28] - bar_values_aligned_linear[28], bar_values_aligned_pangenome[30] - bar_values_aligned_linear[30], bar_values_aligned_pangenome[32] - bar_values_aligned_linear[32], bar_values_aligned_pangenome[34] - bar_values_aligned_linear[34], bar_values_aligned_pangenome[36] - bar_values_aligned_linear[36],bar_values_aligned_pangenome[38] - bar_values_aligned_linear[38]]
EAA_bar_aligned_values = [bar_values_aligned_pangenome[41] - bar_values_aligned_linear[41],bar_values_aligned_pangenome[43] - bar_values_aligned_linear[43], bar_values_aligned_pangenome[45] - bar_values_aligned_linear[45], bar_values_aligned_pangenome[47] - bar_values_aligned_linear[47], bar_values_aligned_pangenome[49] - bar_values_aligned_linear[49], bar_values_aligned_pangenome[51] - bar_values_aligned_linear[51], bar_values_aligned_pangenome[53] - bar_values_aligned_linear[53], bar_values_aligned_pangenome[55] - bar_values_aligned_linear[55], bar_values_aligned_pangenome[57] - bar_values_aligned_linear[57], bar_values_aligned_pangenome[59] - bar_values_aligned_linear[59], bar_values_aligned_pangenome[40] - bar_values_aligned_linear[40], bar_values_aligned_pangenome[42] - bar_values_aligned_linear[42],bar_values_aligned_pangenome[44] - bar_values_aligned_linear[44], bar_values_aligned_pangenome[46] - bar_values_aligned_linear[46], bar_values_aligned_pangenome[48] - bar_values_aligned_linear[48], bar_values_aligned_pangenome[50] - bar_values_aligned_linear[50], bar_values_aligned_pangenome[52] - bar_values_aligned_linear[52], bar_values_aligned_pangenome[54] - bar_values_aligned_linear[54], bar_values_aligned_pangenome[56] - bar_values_aligned_linear[56], bar_values_aligned_pangenome[58] - bar_values_aligned_linear[58]]

fig, (tumour_aligned_graph, normal_aligned_graph) = plt.subplots(2, 1, sharex = True, figsize = (70,50)) 
colors = ['royalblue','violet', 'yellowgreen'] 

# Adds the tumour bar graphs 
AA_tumour_aligned_bars = tumour_aligned_graph.bar(X_AA, AA_tumour_aligned_data, width = bar_width, color= colors[0], bottom = 0)
EA_tumour_aligned_bars = tumour_aligned_graph.bar(X_EA, EA_tumour_aligned_data, width = bar_width, color=colors[1], bottom = 0)
EAA_tumour_aligned_bars = tumour_aligned_graph.bar(X_EAA, EAA_tumour_aligned_data, width = bar_width, color=colors[2], bottom = 0)

tumour_aligned_graph.set_title('Bladder/Tumour', fontsize =40)
tumour_aligned_graph.set_ylabel('Percentage of Unaligned Reads', fontsize = 40)
tumour_aligned_graph.tick_params(axis = 'y', labelsize=30)
# Organizes the bars so they go to the edge rather than being centered on the plot (makes them a bit more spaced out)
tumour_aligned_graph.set_xlim(x_total[0] - 3, x_total[-1] + 5)

# Adds the bar values on the tumour bar graphs 
tumour_aligned_graph.axhline(0, color = 'black', linewidth=1.0, linestyle='--')
add_value_labels(tumour_aligned_graph, AA_tumour_aligned_bars, AA_bar_aligned_values[:10], font_size = 20)
add_value_labels(tumour_aligned_graph, EA_tumour_aligned_bars, EA_bar_aligned_values [:10], font_size = 20)
add_value_labels(tumour_aligned_graph, EAA_tumour_aligned_bars, EAA_bar_aligned_values [:10], font_size = 20)

AA_normal_aligned_bars = normal_aligned_graph.bar(X_AA, AA_normal_aligned_data, width = bar_width, color= colors[0])
EA_normal_aligned_bars = normal_aligned_graph.bar(X_EA, EA_normal_aligned_data, width = bar_width, color= colors[1])
EAA_normal_aligned_bars = normal_aligned_graph.bar(X_EAA, EAA_normal_aligned_data, width = bar_width, color= colors[2])

normal_aligned_graph.axhline(0, color = 'black', linewidth=1.0, linestyle='--')
add_value_labels(normal_aligned_graph, AA_normal_aligned_bars, AA_bar_aligned_values[10:], font_size = 20)
add_value_labels(normal_aligned_graph, EA_normal_aligned_bars, EA_bar_aligned_values [10:], font_size = 20)
add_value_labels(normal_aligned_graph, EAA_normal_aligned_bars, EAA_bar_aligned_values [10:], font_size = 20)

normal_aligned_graph.set_title('Blood', fontsize=40)
normal_aligned_graph.set_ylabel('Percentage of Unaligned Reads', fontsize=40)
normal_aligned_graph.set_xlabel("Patient", fontsize=40)

# normal_aligned_graph.tick_params(axis ='x', labelsize = 20)
normal_aligned_graph.tick_params(axis = 'y', labelsize=30)
# normal_aligned_graph.set_xticks(x_total)
# normal_aligned_graph.set_xticklabels(total_patient_ID)

patient_legend = [
    plt.Rectangle((0,0), 1, 3, color = 'royalblue', label = 'AA Patients'), 
    plt.Rectangle((0,0), 1, 3, color = 'violet', label='EA Patients'),
    plt.Rectangle((0,0), 1, 3, color = 'yellowgreen', label='EAA Patients')
]
tumour_aligned_graph.legend(handles=patient_legend, loc='lower right', title='Patient Ancestry', title_fontsize = 35,bbox_to_anchor=(1,0), fontsize=35)


plt.tight_layout()

# Saves output to pdf, if png needed can just change ending to .png
plt.savefig('Percentage_of_Unaligned_Reads_in_Pangenome_vs_Linear_Reference_Genome.pdf')

print("The plot has been successfully saved as Percentage_of_Unaligned_Reads_in_Pangenome_vs_Linear_Reference_Genome.pdf")

reads_properly_paired = [
    data.get('reads properly paired', 0) / data.get('raw total sequence', 1) if data.get('raw total sequence', 1) > 0 else 0 
    for data in linear_reference_aligned_data
]
print(reads_properly_paired)

pangenome_properly_paired = [
    data.get('total properly paired', 0) / data.get('total alignments', 1) if data.get('total alignments', 1) > 0 else 0 
    for data in pangenome_aligned_data
]
print (pangenome_properly_paired)

bar_values_paired_linear = [
    data.get('reads properly paired', 0) / data.get('raw total sequence', 1) if data.get('raw total sequence', 1) > 0 else 0 
    for data in linear_reference_aligned_data
]

bar_values_paired_pangenome = [
    data.get('total properly paired', 0) / data.get('total alignments', 1) if data.get('total alignments', 1) > 0 else 0 
    for data in pangenome_aligned_data
]

# These lists carry the information of each donor and find the delta of the percentage of properly paired values by subtracting the properly paired value in the linear reference from the properly paired value in the pangenome reference 
# The tumour data is the odd numbers since they are the second txt files in the linear_reference_aligned_files and pangenome_aligned_files lists 
AA_normal_paired_data = [pangenome_properly_paired[0] - reads_properly_paired[0], pangenome_properly_paired[2] - reads_properly_paired[2], pangenome_properly_paired[4] -reads_properly_paired[4], pangenome_properly_paired[6] - reads_properly_paired[6], pangenome_properly_paired[8] - reads_properly_paired[8],pangenome_properly_paired[10] - reads_properly_paired[10],pangenome_properly_paired[12] - reads_properly_paired[12], pangenome_properly_paired[14] - reads_properly_paired[14],pangenome_properly_paired[16] - reads_properly_paired[16], pangenome_properly_paired[18] - reads_properly_paired[18]]
AA_tumour_paired_data = [pangenome_properly_paired[1] - reads_properly_paired[1], pangenome_properly_paired[3] - reads_properly_paired[3], pangenome_properly_paired[5] - reads_properly_paired[5], pangenome_properly_paired[7] - reads_properly_paired[7], pangenome_properly_paired[9] - reads_properly_paired[9],pangenome_properly_paired[11] - reads_properly_paired[11],pangenome_properly_paired[13] - reads_properly_paired[13],pangenome_properly_paired[15] - reads_properly_paired[15], pangenome_properly_paired[17] - reads_properly_paired[17],pangenome_properly_paired[19] - reads_properly_paired[19]]
EA_normal_paired_data =[pangenome_properly_paired[20] - reads_properly_paired[20], pangenome_properly_paired[22] - reads_properly_paired[22], pangenome_properly_paired[24] - reads_properly_paired[24], pangenome_properly_paired[26] - reads_properly_paired[26], pangenome_properly_paired[28] - reads_properly_paired[28],pangenome_properly_paired[30] - reads_properly_paired[30],pangenome_properly_paired[32] - reads_properly_paired[32],pangenome_properly_paired[34] - reads_properly_paired[34],pangenome_properly_paired[36] - reads_properly_paired[36],pangenome_properly_paired[38] - reads_properly_paired[38]]
EA_tumour_paired_data =[pangenome_properly_paired[21] - reads_properly_paired[21], pangenome_properly_paired[23] - reads_properly_paired[23], pangenome_properly_paired[25] - reads_properly_paired[25],pangenome_properly_paired[27] - reads_properly_paired[27],pangenome_properly_paired[29] - reads_properly_paired[29],pangenome_properly_paired[31] - reads_properly_paired[31],pangenome_properly_paired[33] - reads_properly_paired[33],pangenome_properly_paired[35] - reads_properly_paired[35],pangenome_properly_paired[37] - reads_properly_paired[37],pangenome_properly_paired[39] - reads_properly_paired[39]]
EAA_normal_paired_data= [pangenome_properly_paired[40] - reads_properly_paired[40], pangenome_properly_paired[42] - reads_properly_paired[42],pangenome_properly_paired[44] - reads_properly_paired[44],pangenome_properly_paired[46] - reads_properly_paired[46],pangenome_properly_paired[48] - reads_properly_paired[48],pangenome_properly_paired[50] - reads_properly_paired[50],pangenome_properly_paired[52] - reads_properly_paired[52],pangenome_properly_paired[54] - reads_properly_paired[54],pangenome_properly_paired[56] - reads_properly_paired[56],pangenome_properly_paired[58] - reads_properly_paired[58]]
EAA_tumour_paired_data = [pangenome_properly_paired[41] - reads_properly_paired[41], pangenome_properly_paired[43] - reads_properly_paired[43],pangenome_properly_paired[45] - reads_properly_paired[45],pangenome_properly_paired[47] - reads_properly_paired[47],pangenome_properly_paired[49] - reads_properly_paired[49],pangenome_properly_paired[51] - reads_properly_paired[51],pangenome_properly_paired[53] - reads_properly_paired[53],pangenome_properly_paired[55] - reads_properly_paired[55],pangenome_properly_paired[57] - reads_properly_paired[57],pangenome_properly_paired[59] - reads_properly_paired[59]]   

AA_bar_paired_values = [bar_values_paired_pangenome[1] - bar_values_paired_linear[1], bar_values_paired_pangenome[3] - bar_values_paired_linear[3], bar_values_paired_pangenome[5] - bar_values_paired_linear[5], bar_values_paired_pangenome[7] - bar_values_paired_linear[7], bar_values_paired_pangenome[9]- bar_values_paired_linear[9], bar_values_paired_pangenome[11] - bar_values_paired_linear[11], bar_values_paired_pangenome[13]- bar_values_paired_linear[13], bar_values_paired_pangenome[15] - bar_values_paired_linear[15], bar_values_paired_pangenome[17] - bar_values_paired_linear[17], bar_values_paired_pangenome[19] - bar_values_paired_linear[19], bar_values_paired_pangenome[0] - bar_values_paired_linear[0], bar_values_paired_pangenome[2] - bar_values_paired_linear[2], bar_values_paired_pangenome[4] - bar_values_paired_linear[4], bar_values_paired_pangenome[6] - bar_values_paired_linear[6], bar_values_paired_pangenome[8] - bar_values_paired_linear[8], bar_values_paired_pangenome[10] - bar_values_paired_linear[10], bar_values_paired_pangenome[12] - bar_values_paired_linear[12], bar_values_paired_pangenome[14] - bar_values_paired_linear[14], bar_values_paired_pangenome[16] - bar_values_paired_linear[16], bar_values_paired_pangenome[18] - bar_values_paired_linear[18]]
EA_bar_paired_values = [bar_values_paired_pangenome[21] - bar_values_paired_linear[21], bar_values_paired_pangenome[23] - bar_values_paired_linear[23], bar_values_paired_pangenome[25] - bar_values_paired_linear[25], bar_values_paired_pangenome[27] - bar_values_paired_linear[27], bar_values_paired_pangenome[29] - bar_values_paired_linear[29], bar_values_paired_pangenome[31] - bar_values_paired_linear[31], bar_values_paired_pangenome[33] - bar_values_paired_linear[33], bar_values_paired_pangenome[35] - bar_values_paired_linear[35], bar_values_paired_pangenome[37] - bar_values_paired_linear[37], bar_values_paired_pangenome[39] - bar_values_paired_linear[39], bar_values_paired_pangenome[20] - bar_values_paired_linear[20], bar_values_paired_pangenome[22] - bar_values_paired_linear[22], bar_values_paired_pangenome[24] - bar_values_paired_linear[24], bar_values_paired_pangenome[26] - bar_values_paired_linear[26], bar_values_paired_pangenome[28] - bar_values_paired_linear[28], bar_values_paired_pangenome[30] - bar_values_paired_linear[30], bar_values_paired_pangenome[32] - bar_values_paired_linear[32], bar_values_paired_pangenome[34] - bar_values_paired_linear[34], bar_values_paired_pangenome[36] - bar_values_paired_linear[36],bar_values_paired_pangenome[38] - bar_values_paired_linear[38]]
EAA_bar_paired_values = [bar_values_paired_pangenome[41] - bar_values_paired_linear[41],bar_values_paired_pangenome[43] - bar_values_paired_linear[43], bar_values_paired_pangenome[45] - bar_values_paired_linear[45], bar_values_paired_pangenome[47] - bar_values_paired_linear[47], bar_values_paired_pangenome[49] - bar_values_paired_linear[49], bar_values_paired_pangenome[51] - bar_values_paired_linear[51], bar_values_paired_pangenome[53] - bar_values_paired_linear[53], bar_values_paired_pangenome[55] - bar_values_paired_linear[55], bar_values_paired_pangenome[57] - bar_values_paired_linear[57], bar_values_paired_pangenome[59] - bar_values_paired_linear[59], bar_values_paired_pangenome[40] - bar_values_paired_linear[40], bar_values_paired_pangenome[42] - bar_values_paired_linear[42],bar_values_paired_pangenome[44] - bar_values_paired_linear[44], bar_values_paired_pangenome[46] - bar_values_paired_linear[46], bar_values_paired_pangenome[48] - bar_values_paired_linear[48], bar_values_paired_pangenome[50] - bar_values_paired_linear[50], bar_values_paired_pangenome[52] - bar_values_paired_linear[52], bar_values_paired_pangenome[54] - bar_values_paired_linear[54], bar_values_paired_pangenome[56] - bar_values_paired_linear[56], bar_values_paired_pangenome[58] - bar_values_paired_linear[58]]

fig, (tumour_paired_graph, normal_paired_graph) = plt.subplots(2, 1, sharex = True, figsize = (70,50)) 

AA_tumour_paired_bars = tumour_paired_graph.bar(X_AA, AA_tumour_paired_data, width = bar_width, color= colors[0], bottom = 0)
EA_tumour_paired_bars = tumour_paired_graph.bar(X_EA, EA_tumour_paired_data, width = bar_width, color=colors[1], bottom = 0)
EAA_tumour_paired_bars = tumour_paired_graph.bar(X_EAA, EAA_tumour_paired_data, width = bar_width, color=colors[2], bottom = 0)

tumour_paired_graph.set_title('Bladder/Tumour', fontsize =40)
tumour_paired_graph.set_ylabel('Proportion of Properly Paired Reads', fontsize = 40)
tumour_paired_graph.tick_params(axis = 'y', labelsize=30)
tumour_paired_graph.set_xlim(x_total[0] - 3, x_total[-1] + 5)

tumour_paired_graph.axhline(0, color = 'black', linewidth=1.0, linestyle='--')
add_value_labels(tumour_paired_graph, AA_tumour_paired_bars, AA_bar_paired_values[:10], font_size = 20)
add_value_labels(tumour_paired_graph, EA_tumour_paired_bars, EA_bar_paired_values [:10], font_size = 20)
add_value_labels(tumour_paired_graph, EAA_tumour_paired_bars, EAA_bar_paired_values [:10], font_size = 20)

AA_normal_paired_bars = normal_paired_graph.bar(X_AA, AA_normal_paired_data, width = bar_width, color= colors[0], bottom = 0)
EA_normal_paired_bars = normal_paired_graph.bar(X_EA, EA_normal_paired_data, width = bar_width, color= colors[1], bottom = 0)
EAA_normal_paired_bars = normal_paired_graph.bar(X_EAA, EAA_normal_paired_data, width = bar_width, color= colors[2], bottom = 0)

normal_paired_graph.axhline(0, color = 'black', linewidth=1.0, linestyle='--')
add_value_labels(normal_paired_graph, AA_normal_paired_bars, AA_bar_paired_values[10:], font_size = 20)
add_value_labels(normal_paired_graph, EA_normal_paired_bars, EA_bar_paired_values [10:], font_size = 20)
add_value_labels(normal_paired_graph, EAA_normal_paired_bars, EAA_bar_paired_values [10:], font_size = 20)

normal_paired_graph.set_title('Blood', fontsize=40)
normal_paired_graph.set_ylabel('Proportion of Properly Paired Reads', fontsize=40)
normal_paired_graph.set_xlabel("Patient", fontsize=40)

# normal_paired_graph.tick_params(axis ='x', labelsize = 20)
normal_paired_graph.tick_params(axis = 'y', labelsize=30)
# normal_paired_graph.set_xticks(x_total)
# normal_paired_graph.set_xticklabels(total_patient_ID)

tumour_paired_graph.legend(handles=patient_legend, loc='lower right', title='Patient Ancestry', title_fontsize = 35,bbox_to_anchor=(1,0), fontsize=35)

plt.tight_layout()
# Saves output to pdf, if png needed can just change ending to .png
plt.savefig('Proportion_of_Properly_Paired_Reads_in_Pangenome_vs_Linear_Reference_Genome.pdf')


print("The plot has been successfully saved as Proportion_of_Properly_Paired_Reads_in_Pangenome_vs_Linear_Reference_Genome.pdf")

# For the linear reference, samtools does not provide a mapping quality value in the stats txt output, these are the outputs that correspond to each DONOR:SAMPLE (for example DONOR1:BLOOD, DONOR1:BLADDER, etc) mapping quality found from running the mapping quality command 
AA_mapq = [53.735, 53.9429, 49.7885, 53.4471, 53.2826, 53.5952, 53.8863, 54.025, 50.0626, 52.7562, 53.1873, 53.7629, 53.5683, 53.578, 54.1198, 54.1529, 52.6802, 53.0379, 52.7718, 52.8663]
EA_mapq = [52.8471, 53.058, 52.7379, 53.1318, 52.5516, 52.6087, 53.1323, 53.1617, 52.9405, 53.0686, 52.673, 53.9567, 52.6269, 52.5575, 52.992, 52.8285, 53.1357, 53.1127, 52.6633, 53.0461]
EAA_mapq = [53.8758, 54.0445, 53.606, 53.1472, 52.8453, 53.1594, 53.8524, 53.8524, 49.6688, 52.4882, 54.3576, 54.2719, 52.7017, 53.5434, 53.3416, 53.6627, 54.082, 53.6495, 54.1814, 53.2378]
mapping_quality = [
    data.get('mapping quality: mean', 0)
    for data in pangenome_aligned_data
]
print(mapping_quality)

# These lists carry the information of each donor and finds the delta of the mapping quality by subtracting the mapping quality in the linear reference from the mapping quality in the pangenome reference 
# The tumour data is the odd numbers since they are the second txt files in the linear_reference_aligned_files and pangenome_aligned_files lists 
AA_mapq_normal_data = [mapping_quality[0] - AA_mapq[0], mapping_quality[2] - AA_mapq[2], mapping_quality[4] -AA_mapq[4], mapping_quality[6] - AA_mapq[6], mapping_quality[8] - AA_mapq[8],mapping_quality[10] - AA_mapq[10],mapping_quality[12] - AA_mapq[12], mapping_quality[14] - AA_mapq[14],mapping_quality[16] - AA_mapq[16], mapping_quality[18] - AA_mapq[18]]
AA_mapq_tumour_data = [mapping_quality[1] - AA_mapq[1], mapping_quality[3] - AA_mapq[3], mapping_quality[5] - AA_mapq[5], mapping_quality[7] - AA_mapq[7], mapping_quality[9] - AA_mapq[9],mapping_quality[11] - AA_mapq[11],mapping_quality[13] - AA_mapq[13],mapping_quality[15] - AA_mapq[15], mapping_quality[17] - AA_mapq[17],mapping_quality[19] - AA_mapq[19]]
EA_mapq_normal_data =[mapping_quality[20] - EA_mapq[0], mapping_quality[22] - EA_mapq[2], mapping_quality[24] - EA_mapq[4], mapping_quality[26] - EA_mapq[6], mapping_quality[28] - EA_mapq[8],mapping_quality[30] - EA_mapq[10],mapping_quality[32] - EA_mapq[12],mapping_quality[34] - EA_mapq[14],mapping_quality[36] - EA_mapq[16],mapping_quality[38] - EA_mapq[18]]
EA_mapq_tumour_data =[mapping_quality[21] - EA_mapq[1], mapping_quality[23] - EA_mapq[3], mapping_quality[25] - EA_mapq[5],mapping_quality[27] - EA_mapq[7],mapping_quality[29] - EA_mapq[9],mapping_quality[31] - EA_mapq[11],mapping_quality[33] - EA_mapq[13],mapping_quality[35] - EA_mapq[15],mapping_quality[37] - EA_mapq[17],mapping_quality[39] - EA_mapq[19]]
EAA_mapq_normal_data= [mapping_quality[40] - EAA_mapq[0], mapping_quality[42] - EAA_mapq[2],mapping_quality[44] - EAA_mapq[4],mapping_quality[46] - EAA_mapq[6],mapping_quality[48] - EAA_mapq[8],mapping_quality[50] - EAA_mapq[10],mapping_quality[52] - EAA_mapq[12],mapping_quality[54] - EAA_mapq[14],mapping_quality[56] - EAA_mapq[16],mapping_quality[58] - EAA_mapq[18]]
EAA_mapq_tumour_data = [mapping_quality[41] - EAA_mapq[1], mapping_quality[43] - EAA_mapq[3],mapping_quality[45] - EAA_mapq[5],mapping_quality[47] - EAA_mapq[7],mapping_quality[49] - EAA_mapq[9],mapping_quality[51] - EAA_mapq[11],mapping_quality[53] - EAA_mapq[13],mapping_quality[55] - EAA_mapq[15],mapping_quality[57] - EAA_mapq[17],mapping_quality[59] - EAA_mapq[19]]   

AA_bar_mapq_values = [mapping_quality[1] - AA_mapq[1], mapping_quality[3] - AA_mapq[3], mapping_quality[5] - AA_mapq[5], mapping_quality[7] - AA_mapq[7], mapping_quality[9]- AA_mapq[9], mapping_quality[11] - AA_mapq[11], mapping_quality[13]- AA_mapq[13], mapping_quality[15] - AA_mapq[15], mapping_quality[17] - AA_mapq[17], mapping_quality[19] - AA_mapq[19], mapping_quality[0] - AA_mapq[0], mapping_quality[2] - AA_mapq[2], mapping_quality[4] - AA_mapq[4], mapping_quality[6] - AA_mapq[6], mapping_quality[8] - AA_mapq[8], mapping_quality[10] - AA_mapq[10], mapping_quality[12] - AA_mapq[12], mapping_quality[14] - AA_mapq[14], mapping_quality[16] - AA_mapq[16], mapping_quality[18] - AA_mapq[18]]
EA_bar_mapq_values = [mapping_quality[21] - EA_mapq[1], mapping_quality[23] - EA_mapq[3], mapping_quality[25] - EA_mapq[5], mapping_quality[27] - EA_mapq[7], mapping_quality[29] - EA_mapq[9], mapping_quality[31] - EA_mapq[11], mapping_quality[33] - EA_mapq[13], mapping_quality[35] - EA_mapq[15], mapping_quality[37] - EA_mapq[17], mapping_quality[39] - EA_mapq[19], mapping_quality[20] - EA_mapq[0], mapping_quality[22] - EA_mapq[2], mapping_quality[24] - EA_mapq[4], mapping_quality[26] - EA_mapq[6], mapping_quality[28] - EA_mapq[8], mapping_quality[10] - EA_mapq[10], mapping_quality[32] - EA_mapq[12], mapping_quality[34] - EA_mapq[14], mapping_quality[36] - EA_mapq[16],mapping_quality[38] - EA_mapq[18]]
EAA_bar_mapq_values = [mapping_quality[41] - EAA_mapq[1],mapping_quality[43] - EAA_mapq[3], mapping_quality[45] - EAA_mapq[5], mapping_quality[47] - EAA_mapq[7], mapping_quality[49] - EAA_mapq[9], mapping_quality[51] - EAA_mapq[11], mapping_quality[53] - EAA_mapq[13], mapping_quality[55] - EAA_mapq[15], mapping_quality[57] - EAA_mapq[17], mapping_quality[59] - EAA_mapq[19], mapping_quality[40] - EAA_mapq[0], mapping_quality[42] - EAA_mapq[2],mapping_quality[44] - EAA_mapq[4], mapping_quality[46] - EAA_mapq[6], mapping_quality[48] - EAA_mapq[8], mapping_quality[50] - EAA_mapq[10], mapping_quality[52] - EAA_mapq[12], mapping_quality[54] - EAA_mapq[14], mapping_quality[56] - EAA_mapq[16], mapping_quality[58] - EAA_mapq[18]]

fig, (tumour_mapq_graph, normal_mapq_graph) = plt.subplots(2, 1, sharex = True, figsize = (70,50)) 

AA_mapq_tumour_bars = tumour_mapq_graph.bar(X_AA, AA_mapq_tumour_data, width = bar_width, color= colors[0], bottom = 0)
EA_mapq_tumour_bars = tumour_mapq_graph.bar(X_EA, EA_mapq_tumour_data, width = bar_width, color=colors[1], bottom = 0)
EAA_mapq_tumour_bars = tumour_mapq_graph.bar(X_EAA, EAA_mapq_tumour_data, width = bar_width, color=colors[2], bottom = 0)

tumour_mapq_graph.set_title('Bladder/Tumour', fontsize =40)
tumour_mapq_graph.set_ylabel('Mapping Quality', fontsize = 40)
tumour_mapq_graph.tick_params(axis = 'y', labelsize=30)
tumour_mapq_graph.set_xlim(x_total[0] - 3, x_total[-1] + 5)

tumour_mapq_graph.axhline(0, color = 'black', linewidth=1.0, linestyle='--')
add_value_labels(tumour_mapq_graph, AA_mapq_tumour_bars, AA_bar_mapq_values[:10], font_size = 20)
add_value_labels(tumour_mapq_graph, EA_mapq_tumour_bars, EA_bar_mapq_values [:10], font_size = 20)
add_value_labels(tumour_mapq_graph, EAA_mapq_tumour_bars, EAA_bar_mapq_values [:10], font_size = 20)

AA_mapq_normal_bars = normal_mapq_graph.bar(X_AA, AA_mapq_normal_data, width = bar_width, color= colors[0], bottom = 0)
EA_mapq_normal_bars = normal_mapq_graph.bar(X_EA, EA_mapq_normal_data, width = bar_width, color= colors[1], bottom = 0)
EAA_mapq_normal_bars = normal_mapq_graph.bar(X_EAA, EAA_mapq_normal_data, width = bar_width, color= colors[2], bottom = 0)

normal_mapq_graph.axhline(0, color = 'black', linewidth=1.0, linestyle='--')
add_value_labels(normal_mapq_graph, AA_mapq_normal_bars, AA_bar_mapq_values[10:], font_size = 20)
add_value_labels(normal_mapq_graph, EA_mapq_normal_bars, EA_bar_mapq_values [10:], font_size = 20)
add_value_labels(normal_mapq_graph, EAA_mapq_normal_bars, EAA_bar_mapq_values [10:], font_size =20)

normal_mapq_graph.set_title('Blood', fontsize=40)
normal_mapq_graph.set_ylabel('Mapping Quality', fontsize=40)
normal_mapq_graph.set_xlabel("Patient", fontsize=40)

# normal_mapq_graph.tick_params(axis ='x', labelsize = 20)
normal_mapq_graph.tick_params(axis = 'y', labelsize=30)
# normal_mapq_graph.set_xticks(x_total)
# normal_mapq_graph.set_xticklabels(total_patient_ID)

tumour_mapq_graph.legend(handles=patient_legend, loc='lower right', title='Patient Ancestry', title_fontsize = 35,bbox_to_anchor=(1,0), fontsize=35)

plt.tight_layout()
# Saves output to pdf, if png needed can just change ending to .png
plt.savefig('Mapping_Quality_in_Pangenome_vs_Linear_Reference_Genome.pdf')


print("The plot has been successfully saved as Mapping_Quality_in_Pangenome_vs_Linear_Reference_Genome.pdf")