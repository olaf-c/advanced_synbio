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

def split_binary_droplet_string(droplet_string):
    luby = droplet_string[:16]
    error = droplet_string[272:]
    message = droplet_string[16:272]
    return luby, message, error

def rename_keys_and_split_fasta(binary_sequence_dict):
    pattern = re.compile(r'droplet_n(\d+)_.*')
    binary_sequence_dict = {pattern.sub(r'\1', key): value for key, value in binary_sequence_dict.items()}
    for droplet in binary_sequence_dict:
        luby_index_string, droplet_message_string, error_correction_string = split_binary_droplet_string(binary_sequence_dict[droplet]['Binary'])
        binary_sequence_dict[droplet]['LubyIndex'] = luby_index_string
        binary_sequence_dict[droplet]['ErrorCorrection'] = error_correction_string
        binary_sequence_dict[droplet]['DropletMessage'] = droplet_message_string
    binary_sequence_dict = {int(key): value for key, value in binary_sequence_dict.items()}
    return binary_sequence_dict
    
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
    binary_droplet_dict = rename_keys_and_split_fasta(binary_droplet_dict)
    return binary_droplet_dict, luby_blocks

def main():
    droplet_sequence_dict, luby_blocks = load_in_files()
    binary_droplet_dict, luby_blocks = process_files(droplet_sequence_dict, luby_blocks)
    binary_string = reverse_luby(binary_droplet_dict, luby_blocks)
    print(convert_to_ascii(binary_string))
if __name__ == "__main__":
    main()
