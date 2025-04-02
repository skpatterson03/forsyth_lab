import gffpandas.gffpandas as gffpd
import pandas as pd
import numpy as np
import argparse

def load_data():
    """Load all required data files"""
    data = {}
    
    # Load GFF files
    data['ctrl_ph7'] = gffpd.read_gff3('../data/sample1_modification.gff')  # Wild Type pH 7
    data['arss_ph7'] = gffpd.read_gff3('../data/sample2_modification.gff')  # ArsS deletion pH 7, cannot detect pH change
    data['ctrl_ph5'] = gffpd.read_gff3('../data/sample5_modification.gff')  # wild type pH 5
    data['arss_ph5'] = gffpd.read_gff3('../data/sample6_modification.gff')  # ArsS deletion pH 5
    
    # Load annotation and TSS files
    data['annot'] = pd.read_csv("../data/hpy_annot.tsv", sep='\t')  # hpy annotation file
    data['annot'] = data['annot'].drop(columns=["seqname", "source", "feature", "score", "frame", "Prot_Size(est)", "attribute"])  # clean
    
    data['tss'] = pd.read_csv("../data/TSS.csv")
    
    # Preprocessing GFF files
    for key in ['ctrl_ph7', 'ctrl_ph5', 'arss_ph7', 'arss_ph5']:
        data[key] = data[key].attributes_to_columns()
        data[key] = data[key][data[key]['IPDRatio'] >= '2']
    
    # Create merged dataframe
    frames = [data['ctrl_ph7'], data['ctrl_ph5'], data['arss_ph7'], data['arss_ph5']]
    keys = ['ctrl pH7', 'ctrl pH5', 'arsS pH7', 'arsS pH5']
    data['merged'] = pd.concat(frames, keys=keys)
    data['merged'] = data['merged'].reset_index(level=0).rename(columns={'level_0': 'sample'})
    
    return data

def create_ctrl5_csv(data):
    """Create ctrl5.csv file"""
    print("Creating ctrl5.csv...")
    motif_counts = data['ctrl_ph5']['motif'].value_counts()
    data['ctrl_ph5'].to_csv('../data/ctrl5.csv')
    print("ctrl5.csv created successfully")

def create_gene_counts_csv(data):
    """Create gene_counts.csv file"""
    print("Creating gene_counts.csv...")
    counts = pd.DataFrame(columns=['gene', '1\' pH7', '1\' pH5', '16\' pH7', '16\' pH5'])
    counts['gene'] = data['annot']['HP_number']
    
    for ind in counts.index:
        counts.loc[ind, '1\' pH7'] = ((data['ctrl_ph7']['start'].between(data['annot'].loc[ind, 'start'], data['annot'].loc[ind, 'end'] + 1))).sum()
        counts.loc[ind, '1\' pH5'] = ((data['ctrl_ph5']['start'].between(data['annot'].loc[ind, 'start'], data['annot'].loc[ind, 'end'] + 1))).sum()
        counts.loc[ind, '16\' pH7'] = ((data['arss_ph7']['start'].between(data['annot'].loc[ind, 'start'], data['annot'].loc[ind, 'end'] + 1))).sum()
        counts.loc[ind, '16\' pH5'] = ((data['arss_ph5']['start'].between(data['annot'].loc[ind, 'start'], data['annot'].loc[ind, 'end'] + 1))).sum()
    
    for ind in counts.index:
        gene_length = data['annot'].loc[ind, 'end'] - data['annot'].loc[ind, 'start']
        counts.loc[ind, '1\' pH7 density'] = counts.loc[ind, '1\' pH7'] / gene_length
        counts.loc[ind, '16\' pH7 density'] = counts.loc[ind, '16\' pH7'] / gene_length
        counts.loc[ind, '1\' pH5 density'] = counts.loc[ind, '1\' pH5'] / gene_length
        counts.loc[ind, '16\' pH5 density'] = counts.loc[ind, '16\' pH5'] / gene_length
    
    counts.to_csv('../data/gene_counts.csv')
    print("gene_counts.csv created successfully")

def create_promoter_methyl_counts_csv(data):
    """Create promoter_methyl_counts.csv file"""
    print("Creating promoter_methyl_counts.csv...")
    primary = data['tss'][(data['tss']['Primary'] == 1)]
    
    def calculate_adjusted_tss(row):
        if row['Strand'] == '+':
            return row['TSS'] - 50  # Subtract 50 bp for the + strand
        elif row['Strand'] == '-':
            return row['TSS'] + 50  # Add 50 bp for the - strand
        else:
            return row['TSS']
    
    primary = primary.copy()  # Create a copy to avoid the SettingWithCopyWarning
    primary['+/- 50 TSS'] = primary.apply(calculate_adjusted_tss, axis=1)
    new_df = primary[['TSS', '+/- 50 TSS', 'Strand', 'Locus_tag']]
    
    p_counts = pd.DataFrame(columns=['gene', '-Val', 'strand', '1\' pH7',
                                     '1\' pH5', '16\' pH5', '16\' pH7'])
    p_counts['gene'] = new_df['Locus_tag']
    p_counts['strand'] = new_df['Strand']
    p_counts['TSS'] = new_df['TSS']
    p_counts['-Val'] = new_df['+/- 50 TSS']
    
    for ind in p_counts.index:
        tss_value = int(p_counts.loc[ind, 'TSS'])
        val_value = int(p_counts.loc[ind, '-Val'])
        strand = p_counts.loc[ind, 'strand']
        # Determine the range based on strand
        if strand == '-':
            range_start = tss_value
            range_end = val_value
        elif strand == '+':
            range_start = val_value
            range_end = tss_value
        # Update p_counts with the counts from the ranges
        p_counts.loc[ind, "1' pH7"] = (data['ctrl_ph7']['start'].between(range_start, range_end)).sum()
        p_counts.loc[ind, "1' pH5"] = (data['ctrl_ph5']['start'].between(range_start, range_end)).sum()
        p_counts.loc[ind, "16' pH7"] = (data['arss_ph7']['start'].between(range_start, range_end)).sum()
        p_counts.loc[ind, "16' pH5"] = (data['arss_ph5']['start'].between(range_start, range_end)).sum()
    
    p_counts.to_csv('../data/promoter_methyl_counts.csv', index=False)
    print("promoter_methyl_counts.csv created successfully")
    
    return p_counts

def find_promoter_id(start_val, p_counts):
    """Find the promoter ID for a given genomic position"""
    match = p_counts[((p_counts['-Val'] <= start_val) & (p_counts['TSS'] >= start_val) & (p_counts['strand'] == '+')) |
                     ((p_counts['TSS'] <= start_val) & (p_counts['-Val'] >= start_val) & (p_counts['strand'] == '-'))]
    # Return the first matching gene, or None if no match
    if not match.empty:
        return match.iloc[0]['gene']
    return None


def find_coding_reg_id(start_val, annot):
    """Find the coding region ID for a given genomic position"""
    match = annot[((annot['start'] <= start_val) & (annot['end'] >= start_val))]
    # Return the first matching gene, or None if no match
    if not match.empty:
        return match.iloc[0]['HP_number']
    return None


def create_methylation_df_csv(data, p_counts):
    """Create methylation_df.csv file with genome classification"""
    print("Creating methylation_df.csv...")
    
    # Combine promoter DataFrames, + and -
    prom_minus = p_counts[p_counts['strand'] == '-']
    prom_plus = p_counts[p_counts['strand'] == '+']
    
    prom_combined = pd.concat([prom_plus, prom_minus], ignore_index=True)
    
    # Extract the necessary arrays for + and - strand promoters
    val_arr_plus = prom_plus['-Val'].to_numpy() if not prom_plus.empty else np.array([0])
    tss_arr_plus = prom_plus['TSS'].to_numpy() if not prom_plus.empty else np.array([0])
    
    val_arr_minus = prom_minus['-Val'].to_numpy() if not prom_minus.empty else np.array([0])
    tss_arr_minus = prom_minus['TSS'].to_numpy() if not prom_minus.empty else np.array([0])
    
    # Coding region arrays
    start_arr = data['annot']['start'].to_numpy()
    end_arr = data['annot']['end'].to_numpy()
    
    # Determine the final genome length based on both promoter and coding region data
    end_promoters = max(
        max(val_arr_plus) if len(val_arr_plus) > 0 else 0, 
        max(val_arr_minus) if len(val_arr_minus) > 0 else 0, 
        max(tss_arr_plus) if len(tss_arr_plus) > 0 else 0, 
        max(tss_arr_minus) if len(tss_arr_minus) > 0 else 0
    )
    end_coding = max(end_arr)
    final_genome_length = max(end_promoters, end_coding) + 100  # Adding padding
    
    # Create an array initialized with zeros (intergenic regions)
    genome_array = np.zeros(final_genome_length + 1, dtype=int)
    
    # First pass: Mark promoter regions as 1
    for index, row in prom_combined.iterrows():
        tss_value = int(row['TSS'])
        val_value = int(row['-Val'])
        strand = row['strand']  # Get strand information
        
        if strand == '+':
            # Mark positions from -Val to TSS for + strand as promoter
            genome_array[val_value:tss_value+1] = 1  # Set promoter region to 1
        elif strand == '-':
            # Mark positions from TSS to -Val for - strand as promoter
            genome_array[tss_value:val_value+1] = 1  # Set promoter region to 1
    
    # Second pass: Mark coding regions as 2, or as 3 if they overlap with promoters
    for index, row in data['annot'].iterrows():
        start = row['start']
        end = row['end']
        
        # If the region overlaps with a promoter (1), mark it as both (3)
        # Otherwise, mark it as a coding region (2)
        for pos in range(start, end+1):
            if pos < len(genome_array):
                if genome_array[pos] == 1:
                    genome_array[pos] = 3  # Both promoter and coding region
                elif genome_array[pos] == 0:
                    genome_array[pos] = 2  # Mark as coding region
    
    def classify_coding_region(row):
        start_location = row['start']
        if start_location < len(genome_array):
            if genome_array[start_location] == 1:
                return 'Promoter Region'
            elif genome_array[start_location] == 2:
                return 'Coding Region'
            elif genome_array[start_location] == 3:
                return 'Promoter and Coding Region'
            else:
                return 'Intergenic Region'
        else:
            return np.nan
    
    # Create a copy of merged data to avoid modifying the original
    methylation_df = data['merged'].copy()
    
    # Apply classification
    methylation_df['Classification'] = methylation_df.apply(classify_coding_region, axis=1)
    
    # Add promoter IDs
    methylation_df['promoter_ID'] = methylation_df['start'].apply(
        find_promoter_id, p_counts=p_counts
    )
    
    # Add coding region IDs
    methylation_df['CodingRegion_ID'] = methylation_df['start'].apply(
        find_coding_reg_id, annot=data['annot']
    )
    
    # Try to add gene categories if the file exists
    try:
        gene_categories = pd.read_csv('../data/hpy_annot.csv')
        
        # Merge coding region gene types
        methylation_df = pd.merge(
            methylation_df,
            gene_categories[['HP_number', 'gene_type']],
            left_on='CodingRegion_ID',
            right_on='HP_number',
            how='left'
        )
        methylation_df.rename(columns={'gene_type': 'CodingRegion_gene_type'}, inplace=True)
        
        # Merge promoter gene types
        methylation_df = pd.merge(
            methylation_df,
            gene_categories[['HP_number', 'gene_type']],
            left_on='promoter_ID',
            right_on='HP_number',
            how='left'
        )
        methylation_df.rename(columns={'gene_type': 'Promoter_gene_type'}, inplace=True)
        
        # Clean up column duplicates and format gene types
        methylation_df = methylation_df.loc[:, ~methylation_df.columns.duplicated()]
        
        # Format gene types
        for col in ['CodingRegion_gene_type', 'Promoter_gene_type']:
            if col in methylation_df.columns:
                methylation_df[col] = methylation_df[col].astype(str)
                methylation_df[col] = methylation_df[col].str.replace(r'\r\n', ' ', regex=True)
        
        # Handle any specific type replacements
        if 'CodingRegion_gene_type' in methylation_df.columns:
            methylation_df['CodingRegion_gene_type'] = methylation_df['CodingRegion_gene_type'].replace({
                'CELLULAR PROCESSE': 'CELLULAR PROCESSES',
                'HYPOTHETICAL PROTEIEN': 'HYPOTHETICAL PROTEIN'
            })
        
        # Drop any remnant HP_number columns
        methylation_df.drop(columns=['HP_number_x', 'HP_number_y'], 
                          inplace=True, errors='ignore')
            
    except Exception as e:
        print(f"Warning: Could not add gene categories: {e}")
    
    # Save the file
    methylation_df.to_csv('../data/methylation_df.csv')
    print("methylation_df.csv created successfully")
    
    return methylation_df


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Create analysis files.')
    parser.add_argument('--files', nargs='+', choices=['ctrl5', 'gene_counts', 'promoter_methyl_counts', 'methylation_df', 'all'],
                        default=['all'], help='Specify which files to create')
    
    args = parser.parse_args()
    files_to_create = args.files
    
    # If 'all' is in the list, create all files
    if 'all' in files_to_create:
        files_to_create = ['ctrl5', 'gene_counts', 'promoter_methyl_counts', 'methylation_df']
    
    # Load data
    print("Loading data...")
    data = load_data()
    print("Data loaded successfully")
    
    # Create files based on arguments
    p_counts = None
    
    if 'ctrl5' in files_to_create:
        create_ctrl5_csv(data)
    
    if 'gene_counts' in files_to_create:
        create_gene_counts_csv(data)
    
    if 'promoter_methyl_counts' in files_to_create or 'methylation_df' in files_to_create:
        p_counts = create_promoter_methyl_counts_csv(data)
    
    if 'methylation_df' in files_to_create:
        create_methylation_df_csv(data, p_counts)
    
    print("All requested files have been created.")

if __name__ == "__main__":
    main()