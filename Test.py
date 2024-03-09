import csv
from Bio.Align import substitution_matrices
from Bio import pairwise2, SeqIO
from tqdm import tqdm
import concurrent.futures

substitution_matrix = substitution_matrices.load('BLOSUM62')  # use matrix names to load
import pandas as pd
from multiprocessing import cpu_count


def calculate_sequence_identity(seq1, seq2):
    alignments = pairwise2.align.globalds(seq1, seq2, substitution_matrix, -0.5, -0.3)
    best_alignment = alignments[0]  # Assuming you want the best alignment
    aligned_seq1, aligned_seq2, score, start, end = best_alignment

    # Calculate sequence identity
    seq_len = len(aligned_seq1)
    identical_positions = sum(1 for i in range(seq_len) if aligned_seq1[i] == aligned_seq2[i])
    identity_percent = (100 * identical_positions / ((len(aligned_seq1) + len(aligned_seq2)) / 2))

    return identity_percent


def load_sequences_from_csv(file_path):
    """Load protein sequences from a CSV file."""
    df = pd.read_csv(file_path)
    return df


def check_duplicates(dataset1, dataset2):
    # Step 1: Merge the DataFrames based on the key column
    merged_df = pd.merge(dataset1, dataset2, on='Sequence', how='inner')

    # Step 2: Identify Duplicates
    duplicates = merged_df[merged_df.duplicated(subset=['Sequence'], keep=False)]

    # Step 3: Filter Duplicate Entries
    if not duplicates.empty:
        return True
    else:
        return False


def remove_duplicates(dataset1, dataset2):
    # Step 1: Merge the DataFrames based on the key column
    merged_df = pd.merge(dataset1, dataset2, on='Sequence', how='inner')

    # Step 2: Identify Duplicates
    duplicates = merged_df[merged_df.duplicated(subset=['Sequence'], keep=False)]
    # Step 3: Filter Duplicate Entries
    if not duplicates.empty:
        # Remove duplicates from dataset1
        dataset1_cleaned = dataset1[~dataset1['Sequence'].isin(duplicates['Sequence'])]
        return dataset1_cleaned
    else:
        return dataset1  # No duplicates found, return dataset1 as is


def check_sequence_validity(seq1, dataset2_sequences, threshold):
    for seq2 in dataset2_sequences:
        score = calculate_sequence_identity(seq1, seq2)
        if score > threshold:
            return False
    return True


def check_sequence_validity_parallel(seq1, dataset2_sequences, threshold):
    num_chunks = cpu_count() - 1 # Number of chunks equal to number of CPU cores
    chunk_size = len(dataset2_sequences) // num_chunks
    dataset2_chunks = [dataset2_sequences[i:i + chunk_size] for i in range(0, len(dataset2_sequences), chunk_size)]
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_chunks) as executor:
        # Submit tasks to the executor
        futures = [executor.submit(check_sequence_validity, seq1, chunk, threshold) for chunk in dataset2_chunks]

        result = True
        for future in concurrent.futures.as_completed(futures):
            result = result and (future.result())  # Store the result in the results dictionary
            if not result:
                return False
    return result


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


def main():
    # Load protein sequences from CSV files
    dataset1_sequences = filter_critical_sequences(load_sequences_from_csv("Dataset/deeploc_data_sub_part_1.csv")) # change the file name
    reference_dataset = load_sequences_from_csv("Dataset/Swissprot_Train_Validation_dataset.csv")
    filtered_dataset1 = remove_duplicates(dataset1_sequences, reference_dataset)
    dataset1_seq = filtered_dataset1['Sequence'].tolist()
    dataset2_sequences = load_sequences_from_csv("Dataset/hpa_testset.csv")
    dataset2_seq = dataset2_sequences['fasta'].tolist()
    new_sequences = pd.DataFrame(columns=filtered_dataset1.columns)
    # Calculate overall similarity between datasets
    for index, row in enumerate(tqdm(dataset1_seq)):
        seq1_index = index
        seq1 = row
        if check_sequence_validity_parallel(seq1, dataset2_seq, 30):
            print("Entry", seq1_index, "is unique")
            new_sequences = pd.concat([new_sequences, pd.DataFrame([filtered_dataset1.iloc[seq1_index]])],
                                      ignore_index=True)
            new_sequences.to_csv("Dataset/new_entries_subpart_1.csv", index=False)
    new_sequences.to_csv("Dataset/new_entries_subpart_1.csv", index=False)


if __name__ == "__main__":
    main()
