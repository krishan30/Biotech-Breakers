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

substitution_matrix = substitution_matrices.load('BLOSUM62')  # use matrix names to load
import pandas as pd
from multiprocessing import cpu_count
# Load protein sequences from CSV files
dataset1_sequences = filter_critical_sequences(load_sequences_from_csv("Dataset/deeploc_data.csv"))
# Skip the first 200 entries
dataset1_sequences = dataset1_sequences.iloc[200:]
# Divide the dataset into 5 equal parts
dataset_parts = np.array_split(dataset1_sequences, 5)
dataset_parts[0].to_csv("Dataset/deeploc_data_sub_part_1.csv", index=False)
dataset_parts[1].to_csv("Dataset/deeploc_data_sub_part_2.csv", index=False)
dataset_parts[2].to_csv("Dataset/deeploc_data_sub_part_3.csv", index=False)
dataset_parts[3].to_csv("Dataset/deeploc_data_sub_part_4.csv", index=False)
dataset_parts[4].to_csv("Dataset/deeploc_data_sub_part_5.csv", index=False)