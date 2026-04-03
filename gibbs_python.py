import sys
import copy

def calculate_q_matrix(site_sequences, pseudocount, symbol_alphabet):
    number_of_sequences = len(site_sequences)
    number_of_symbols = len(symbol_alphabet)
    site_width = len(site_sequences[0])

    count_matrix = []
    count_row = [0] * site_width

    for i in range(number_of_symbols):
        count_matrix.append(count_row.copy())

    q_matrix = copy.deepcopy(count_matrix)

    # Calculate count matrix
    for position in range(site_width):
        for symbol_index in range(number_of_symbols):
            count = 0
            for index in range(number_of_sequences):
                if site_sequences[index][position] == symbol_alphabet[symbol_index]:
                    count += 1
            count_matrix[symbol_index][position] = count

    # Calculate q matrix with pseudocount
    for position in range(site_width):
        for symbol_index in range(number_of_symbols):
            q_matrix[symbol_index][position] = (
                count_matrix[symbol_index][position] + pseudocount
            ) / (number_of_sequences + number_of_symbols * pseudocount)

    return count_matrix, q_matrix

# ------------------------
if __name__ == "__main__":
    # 1. Read sequences from file
    filename = sys.argv[1]

    with open(filename) as f:
        site_sequences = [line.strip() for line in f if line.strip()]

    # 2. DEBUG: show how many sequences and their length
    print("Number of sequences:", len(site_sequences))
    print("Original sequence length:", len(site_sequences[0]))

    # 3. Extract motif window (6 bp, e.g. positions 5-10)
    motif_width = 6
    site_sequences = [seq[5:5+motif_width] for seq in site_sequences]

    print("Motif window sequences:", site_sequences)

    pseudocount = 1
    symbol_alphabet = ['A', 'C', 'G', 'T']

    # 4. Calculate matrices
    count_matrix, q_matrix = calculate_q_matrix(
        site_sequences,
        pseudocount,
        symbol_alphabet
    )

    # 5. Print results nicely
    print("\nCount Matrix (A,C,G,T rows):")
    for row in count_matrix:
        print(row)

    print("\nQ Matrix (probabilities, A,C,G,T rows):")
    for row in q_matrix:
        print([round(x, 3) for x in row])
