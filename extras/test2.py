import os
import sys
import codecs

def remove_bom(input_file, output_file):
    with codecs.open(input_file, 'r', encoding='utf-8-sig') as f_in:
        fasta_data = f_in.read()

    with open(output_file, 'w', encoding='utf-8') as f_out:
        f_out.write(fasta_data)

    print("BOM removed and saved to", output_file)

def process_folder(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for filename in os.listdir(input_folder):
        if filename.endswith(".txt"):
            input_file = os.path.join(input_folder, filename)
            output_file = os.path.join(output_folder, filename)
            remove_bom(input_file, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python remove_bom_from_fasta_folder.py <input_folder_with_bom> <output_folder_without_bom>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    process_folder(input_folder, output_folder)
