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

 
# '0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000'
# '0110011001110101011000110110001000110000011100110111010101101101011110100111101101110100001100000111110101111001011100010111011001111101011110110111000101100100011111010111101001111011001110000011000001111010011000110011000001111101011001000011000001110010'
# '0110011001110101011000110110001000110000011100110111010101101101011110100111101101110100001100000111110101111001011100010111011001111101011110110111000101100100011111010111101001111011001110000011000001111010011000110011000001111101011001000011000001110010'
#
# '0110011001110101011000110110001000110000011100110111010101101101011110100111101101110100001100000111110101111001011100010111011001111101011110110111000101100100011111010111101001111011001110000011000001111010011000110011000001111101011001000011000001110010'
# '0110110100110010000011100111010100111101001001110111101001101000011000100111001000100001010100000110010100110000011110100110001101110001010010100100110001101000001101000100110101111111011000000010101000110000011001110111110101100101011011110110111101100111'
# '0000101101000111011011010001011100001101010101000000111100000101000110000000100101010101011000000001100001001001000010110001010100001100001100010011110100001100010010010011011100000100010110000001101001001010000001000100110100011000000010110101111100010101'
# '0000101101000111011011010001011100001101010101000000111100000101000110000000100101010101011000000001100001001001000010110001010100001100001100010011110100001100010010010011011100000100010110000001101001001010000001000100110100011000000010110101111100010101'                   
 

def bad_luby(binary_droplet_dict, luby_blocks):
    binary_string = []
    for droplet in luby_blocks:
        binary_message_list = []
        xor_product=[]
        for block in droplet:
            binary_message_list.append(binary_droplet_dict[block])
        for digit in len(binary_message_list[0]):
            number_true = 0
            for block in binary_message_list:
                if block[digit] == 1:
                    number_true = number_true +1
            xor_product.append(number_true%2)
        binary_string.append(xor_product)  

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
