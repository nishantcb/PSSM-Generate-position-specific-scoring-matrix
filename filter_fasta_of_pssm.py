import os
import sys

def usage():
    print("Usage:")
    print("python filter_fasta_by_pssm.py -i <input_fasta> -d <pssm_directory> -o <output_fasta>")
    print("-i: Path to the input FASTA file containing all sequences.")
    print("-d: Path to the directory containing PSSM files.")
    print("-o: Path to the output FASTA file containing sequences with generated PSSM metrics.")

def read_fasta(fasta_path):
    """Reads a FASTA file and returns a dictionary with sequence IDs as keys and sequences as values."""
    fasta_dict = {}
    with open(fasta_path, "r") as f:
        seq_id = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:  # Save the previous sequence
                    fasta_dict[seq_id] = "".join(seq)
                seq_id = line[1:]  # Remove ">" from the header
                seq = []  # Start a new sequence
            else:
                seq.append(line)
        if seq_id:  # Save the last sequence
            fasta_dict[seq_id] = "".join(seq)
    return fasta_dict

def filter_sequences(fasta_dict, pssm_dir):
    """Filters the sequences based on the presence of PSSM files."""
    filtered_dict = {}
    # Get the base filenames of all PSSM files (without extensions)
    pssm_files = {os.path.splitext(f)[0] for f in os.listdir(pssm_dir) if f.endswith(".pssm")}
    for seq_id, sequence in fasta_dict.items():
        if seq_id in pssm_files:
            filtered_dict[seq_id] = sequence
    return filtered_dict

def write_fasta(fasta_dict, output_path):
    """Writes a dictionary of sequences to a FASTA file."""
    with open(output_path, "w") as f:
        for seq_id, sequence in fasta_dict.items():
            f.write(f">{seq_id}\n")
            f.write(f"{sequence}\n")

if __name__ == "__main__":
    import getopt

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:d:o:h", ["input=", "directory=", "output=", "help"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)

    input_fasta = ""
    pssm_directory = ""
    output_fasta = ""

    for opt, arg in opts:
        if opt in ("-i", "--input"):
            input_fasta = arg
        elif opt in ("-d", "--directory"):
            pssm_directory = arg
        elif opt in ("-o", "--output"):
            output_fasta = arg
        elif opt in ("-h", "--help"):
            usage()
            sys.exit()

    if not input_fasta or not pssm_directory or not output_fasta:
        print("Missing required arguments!")
        usage()
        sys.exit(2)

    # Read input FASTA file
    try:
        fasta_data = read_fasta(input_fasta)
    except Exception as e:
        print(f"Error reading input FASTA file: {e}")
        sys.exit(1)

    # Filter sequences based on PSSM files
    try:
        filtered_sequences = filter_sequences(fasta_data, pssm_directory)
    except Exception as e:
        print(f"Error filtering sequences: {e}")
        sys.exit(1)

    # Write filtered sequences to output FASTA file
    try:
        write_fasta(filtered_sequences, output_fasta)
        print(f"Filtered FASTA file created: {output_fasta}")
    except Exception as e:
        print(f"Error writing output FASTA file: {e}")
        sys.exit(1)

