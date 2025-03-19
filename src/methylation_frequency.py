import pandas as pd
import matplotlib.pyplot as plt
import re

# Function to extract "IPDRatio" values from attributes column
def extract_ipd_ratio(attributes):
    match = re.search(r'IPDRatio=(\d+\.\d+)', attributes)
    return float(match.group(1)) if match else None

# Read GFF file
file_path = "../data/gff_files/m6A_calls.gff"
df = pd.read_csv(file_path, sep='\t', comment='#', header=None)
df.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]

# Extract IPDRatio values
df["IPDRatio"] = df["attributes"].apply(extract_ipd_ratio)

# Filter out missing values
df = df.dropna(subset=["IPDRatio"])

# Define threshold and filter data
threshold = 3
df_high_ipd = df[df["IPDRatio"] > threshold]

# Plot histogram
plt.hist(df_high_ipd["IPDRatio"], bins=30, color='blue', edgecolor='black')
plt.xlabel("IPDRatio")
plt.ylabel("Frequency")
plt.title("Histogram of IPDRatio Values Greater than 3")
plt.show()
