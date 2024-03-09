import csv

import numpy as np
from Bio.Align import substitution_matrices
from Bio import pairwise2, SeqIO
from tqdm import tqdm
import concurrent.futures


def load_sequences_from_csv(file_path):
    """Load protein sequences from a CSV file."""
    df = pd.read_csv(file_path)
    return df
def filter_critical_sequences(dataset):
    filtered_dataset = dataset[(dataset['Extracellular'] == 1) |
                               (dataset['Cell membrane'] == 1) |
                               (dataset['Mitochondrion'] == 1) |
                               (dataset['Plastid'] == 1) |
                               (dataset['Endoplasmic reticulum'] == 1) |
                               (dataset['Lysosome/Vacuole'] == 1) |
                               (dataset['Golgi apparatus'] == 1) |
                               (dataset['Peroxisome'] == 1)]
    return filtered_dataset

def remove_duplicates(dataset1, dataset2):
    # Step 1: Merge the DataFrames based on the key column
    merged_df = pd.concat([dataset1, dataset2],axis=0)

    # Step 2: Identify Duplicates
    duplicates = merged_df[merged_df.duplicated(subset=['Sequence'], keep=False)]
    # Step 3: Filter Duplicate Entries
    if not duplicates.empty:
        # Remove duplicates from dataset1
        dataset1_cleaned = dataset1[~dataset1['Sequence'].isin(duplicates['Sequence'])]
        return dataset1_cleaned
    else:
        return dataset1  # No duplicates found, return dataset1 as is


substitution_matrix = substitution_matrices.load('BLOSUM62')  # use matrix names to load
import pandas as pd
from multiprocessing import cpu_count
# Load protein sequences from CSV files
dataset1_sequences = filter_critical_sequences(load_sequences_from_csv("Dataset/deeploc_data.csv"))
reference_dataset = load_sequences_from_csv("Dataset/Swissprot_Train_Validation_dataset.csv")
filtered_dataset1 = remove_duplicates(dataset1_sequences, reference_dataset)
# Divide the dataset into 5 equal parts
dataset_parts = np.array_split(filtered_dataset1, 5)
dataset_parts[0].to_csv("Dataset/deeploc_data_sub_part_1.csv", index=False)
dataset_parts[1].to_csv("Dataset/deeploc_data_sub_part_2.csv", index=False)
dataset_parts[2].to_csv("Dataset/deeploc_data_sub_part_3.csv", index=False)
dataset_parts[3].to_csv("Dataset/deeploc_data_sub_part_4.csv", index=False)
dataset_parts[4].to_csv("Dataset/deeploc_data_sub_part_5.csv", index=False)