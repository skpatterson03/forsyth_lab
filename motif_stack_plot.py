import pandas as pd
import logomaker
import matplotlib.pyplot as plt

# need to check that this does the same thing as the R file 

# Filepath to the motif data
filepath = '../data/csv_files/RCGDAD_allmotifs.csv'

# Read in the list of RCGDAD motifs
motifs = pd.read_csv(filepath)

# Convert motifs into a frequency matrix (assuming sequences are in a single column)
def create_freq_matrix(motif_sequences):
    counts_df = logomaker.alignment_to_matrix(motif_sequences, to_type='probability')
    return counts_df

# Generate frequency matrix
freq_matrix = create_freq_matrix(motifs.iloc[:, 0])  # Assuming sequences are in the first column

# Create the sequence logo
fig, ax = plt.subplots(figsize=(10, 4))
logomaker.Logo(freq_matrix, ax=ax)
ax.set_ylabel("Probability")
ax.set_xlabel("Position")
plt.show()
