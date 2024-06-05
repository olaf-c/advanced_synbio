# Challenge One
## Advanced Synthetic and Systems Biology
By: Olaf
Github: https://github.com/olaf-c/advanced_synbio-challenge_one 

### Code Overview
##### Input
*  `load_in_files()` 
Files are read in via the `load_in_files()` function. It saves the fasta file as a dictionary and the csv as a dataframe. The function does not currently take custom file names and does not contain error handling. It depends on os to locate the input_files directory, SeqIO to read the fasta file, and pandas to read the csv. 

##### Processing
* `process_files`
  Wrapper function for the subsequent processing functions. 
* `split_luby`
Uses regular expressions to extract patterns from the pandas dataframe. Lambda functions are used to apply this pattern across the dataframe and split the dataframe from one column to 4.
Columns are named: 'OriginalValue', 'DropletNumber', and 'BlockNumbers'. BlockNumbers are isolated in two steps to create a list within the dataframe column that does not have extra '[' or ']' This method is dependent on the pandas and regex modules.
* `convert_seq_to_binary`
Takes the immutable droplet sequence dict (generated from the fasta file) and creates a new dictionary containing the droplet_id as the key and dictionary value that itself contains the original sequence under the key 'Seq' and the calls upon `convert_str_to_binary` to create the equivalent binary string. 
* `convert_str_to_binary`
Takes a string and replaces A, G, T, and C's with binary strings. It relies on a for loop and the following encoding scheme `{"A": '00', "G": '01', "C": '10', "T": '11' }`
* `rename_keys_and_split_fasta`
Takes the binary dictionary and renames the keys to just the number to standardize between the fasta files and the block files. `split_binary_droplet_string` is called to to identify the 'LubyIndex', 'ErrorCorrection', and 'DropletMessage' message segments which are then saved in the dictionary. Renaming the keys depends on regex.
* `split_binary_droplet_string`
Takes a string and slices it at the following locations: 0:16, 272:end, and 16:272 which are expected to be the luby seed, error correction, and message respectively.

##### Reverse Luby
* `reverse_luby`
Identifies droplets containing only one block and saves their respective droplet binary sequences as known. Then loops through the remainder of the droplets solving droplets with only one unknown block by calls to the `bitwise_xor` function to calculate the bitwise XOR value. The droplets are sorted by number of block components to reduce looping time. The resulting strings are concatenated and returned. Depends on pandas. 
* `bitwise_xor`
Takes two strings and iterates through each character. It performs a XOR operation on the respective characters. It obeys the following truth table:

 Input1 | Input2 | Output 
---|---|---
0| 0 | 0  
0| 1   | 1 
1| 0 | 1 
1 | 1   | 0 

Since XOR operations are reversible this should result in the final unknown string.

##### Convert to ASCII
* `convert_to_ascii`
Converts the string into an int and returns the ASCII encoded by the binary without dependencies.

##### Test Functions
Rudimentary unit tests can be found under the tests folder. 

### The String Flag
> In the year 3074, humanity has transcended its organic roots, embracing the boundless possibilities of synthetic biology. This new era, dubbed the "Bioforge Epoch," is characterized by the fusion of technology and biology, where genetic design and artificial life forms are the backbone of civilization. The economic landscape is dominated by Bio-Conglomerates, colossal corporations that wield the power to shape planets and engineer life itself.
Our saga begins on the fringe world of Terranova, a planet teeming with bio-engineered marvels and the frontier of synthetic biological exploration. Here, amidst the towering neon jungles and sprawling megacities, a young bio-engineer named Kaelin discovers an ancient genetic sequence hidden within the DNA of a seemingly mundane plant. This discovery could revolutionize the field of synthetic biology, offering powers beyond imagination, or it could unleash a catastrophe of cosmic proportions.
As Kaelin delves deeper into the mystery, she finds herself entangled in a web of intrigue and danger. The Bio-Conglomerates, fearing the potential disruption to their dominion, deploy their most sophisticated bio-constructs to seize the genetic code. Meanwhile, a faction of radical bio-hackers views Kaelin's discovery as the key to liberating humanity from corporate tyranny, willing to go to any lengths to obtain it.
With the fate of Terranova and the entire Bioforge Epoch hanging in the balance, Kaelin must navigate a world where alliances are as mutable as genetic code, and where her greatest tool is her intellect. Alongside a diverse crew of outcasts and renegades, she embarks on an epic journey across star systems, through the depths of synthetic worlds, and into the heart of what it means to be human in an age where biology is the ultimate technology. ```

### Code:
```python
from Bio import SeqIO
import pandas as pd
import regex as re
import os

def load_in_files(droplet_sequence_fasta = "droplet_sequences.fasta", luby_block_csv = "luby_blocks.csv"):
    # Doing it with a list. Wouldn't work for larger files. Use SeqUI.index or BioSQL instead. E.G.:
    current_directory = os.path.dirname(__file__)
    
    fasta_path = os.path.join(current_directory, '..', 'input_files', droplet_sequence_fasta)
    csv_path = os.path.join(current_directory, '..', 'input_files', luby_block_csv)
    with open(fasta_path) as handle:
        droplet_sequence_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    
    luby_blocks = pd.read_csv(csv_path, header=None)
    return droplet_sequence_dict, luby_blocks

def split_luby(luby_blocks):
    luby_blocks.columns = ["OriginalValue"]
    # Convert string blocks to list of integers
    luby_blocks['DropletNumber'] = luby_blocks['OriginalValue'].str.extract(r'n(\d+)').apply(lambda x: list(map(int, x)))
    luby_blocks['BlockNumbers'] = luby_blocks['OriginalValue'].apply(lambda x: x[x.find('[')+1:x.find(']')] if '[' in x else None)
    luby_blocks['BlockNumbers'] = luby_blocks['BlockNumbers'].str.split(',',expand=False).apply(lambda x: list(map(int, x)))
    return luby_blocks

def split_binary_droplet_string(droplet_string):
    luby = droplet_string[:16]
    error = droplet_string[272:]
    message = droplet_string[16:272]
    return luby, message, error

def split_fasta(binary_sequence_dict):
    pattern = re.compile(r'droplet_n(\d+)_.*')
    binary_sequence_dict = {pattern.sub(r'\1', key): value for key, value in binary_sequence_dict.items()}
    for droplet in binary_sequence_dict:
        luby_index_string, droplet_message_string, error_correction_string = split_binary_droplet_string(binary_sequence_dict[droplet]['Binary'])
        binary_sequence_dict[droplet]['LubyIndex'] = luby_index_string
        binary_sequence_dict[droplet]['ErrorCorrection'] = error_correction_string
        binary_sequence_dict[droplet]['DropletMessage'] = droplet_message_string
    binary_sequence_dict = {int(key): value for key, value in binary_sequence_dict.items()}
    return binary_sequence_dict

def convert_str_to_binary(seq_string):
    encoding_scheme = {"A": '00', "G": '01', "C": '10', "T": '11' }
    binary = ''
    for nucleotide in seq_string:
        binary += encoding_scheme.get(nucleotide, nucleotide)
    return binary

def convert_seq_to_binary(droplet_sequence_dict):
    binary_droplet_dict = {} #original is immutable
    for droplet_id in droplet_sequence_dict:
        binary = convert_str_to_binary(droplet_sequence_dict[droplet_id].seq)  
        binary_droplet_dict[droplet_id] = {"Seq": droplet_sequence_dict[droplet_id].seq, "Binary": binary}
    return binary_droplet_dict

def bitwise_xor(binary_message_1, binary_message_2):
    result_binary = ''
    for x in range(len(binary_message_1)):
        if (binary_message_1[x] == '1' or binary_message_2[x] == '1') and not ((binary_message_1[x] == '1' and binary_message_2[x] == '1')):
            result_binary = result_binary + '1'
        else:
            result_binary = result_binary + '0'
    return result_binary

def reverse_luby(message_droplet_dict, block_to_droplet_df):
    #XOR operation
    binary_string = ''
    solved_block_dict = {}
    #sort the luby_block pandas file by droplets containing the fewest number of block to the most
    block_to_droplet_df = block_to_droplet_df.iloc[block_to_droplet_df['BlockNumbers'].apply(len).argsort()] #works
    max_block_per_drop = [max(x) for x in block_to_droplet_df['BlockNumbers']] #calculate the number of unique blocks
    number_of_blocks = max(max_block_per_drop)+1 #blocks are zero indexed you numpty
    #droplets that don't need XOR saved to one_to_one // solved_block_dict
    one_to_one_droplets = block_to_droplet_df[(block_to_droplet_df['BlockNumbers'].apply(len) == 1)] #works
    for index, droplet in one_to_one_droplets.iterrows(): 
        solved_block_dict[droplet['BlockNumbers'][0]] = message_droplet_dict[droplet['DropletNumber']]['DropletMessage'] #Works
    while len(solved_block_dict)< number_of_blocks:
        for index, droplet in block_to_droplet_df.iterrows():
            #find droplets with all but one block known
            known_list = []
            unknown_list = []
            for block in droplet['BlockNumbers']:
                if block in solved_block_dict:
                    known_list.append(block)
                else:
                    unknown_list.append(block)
            if len(unknown_list) == 1:
                calculated_code = "0"*256
                for block in known_list:
                    calculated_code = bitwise_xor(calculated_code, solved_block_dict[block])
                calculated_code = bitwise_xor(calculated_code, message_droplet_dict[droplet['DropletNumber']]['DropletMessage'] )
        
                solved_block_dict[unknown_list[0]] = calculated_code
    for i in range(number_of_blocks):
        binary_string = binary_string + solved_block_dict[i]
    return binary_string

def convert_to_ascii(binary_string):
    binary_message = int(binary_string, 2)
    ascii_message = binary_message.to_bytes((binary_message.bit_length() + 7) // 8, 'big').decode()
    return ascii_message

def process_files(droplet_sequence_dict, luby_blocks):
    luby_blocks = split_luby(luby_blocks)
    binary_droplet_dict = convert_seq_to_binary(droplet_sequence_dict)
    binary_droplet_dict = split_fasta(binary_droplet_dict)
    return binary_droplet_dict, luby_blocks

def main():
    droplet_sequence_dict, luby_blocks = load_in_files()
    binary_droplet_dict, luby_blocks = process_files(droplet_sequence_dict, luby_blocks)
    binary_string = reverse_luby(binary_droplet_dict, luby_blocks)
    print(convert_to_ascii(binary_string))
if __name__ == "__main__":
    main()

```