import unittest
from package.challenge_one_dna_fountain import load_in_files, convert_seq_to_binary, convert_str_to_binary, split_luby, split_fasta, split_binary_droplet_string, reverse_luby, convert_to_ascii, bitwise_xor

def contains_no_integers(string):
    for char in string:
        if char.isdigit():
            return False
    return True

def contains_no_letters(string):
    for char in string:
        if char.isalpha():
            return False
    return True

class TestDNAFountain(unittest.TestCase):
    def test_load_in_files(self):
        # Read sequences from temporary FASTA file
        sequences, blocks = load_in_files()

        # Check if file isn't empty
        self.assertTrue(sequences)
        self.assertFalse(blocks.empty)
    def test_convert_seq_to_binary(self):
        sequences, _ = load_in_files()
        binary_dict = convert_seq_to_binary(sequences)
        self.assertTrue(binary_dict)
        for entry in binary_dict:
            self.assertTrue(contains_no_integers(binary_dict[entry]['Seq']))
            self.assertTrue(contains_no_letters(binary_dict[entry]['Binary']))
    def test_split_fasta(self):
        sequences, _ = load_in_files()
        binary_dict = convert_seq_to_binary(sequences)
        split_dict = split_fasta(binary_dict)
        self.assertTrue(split_dict)
        for entry in split_dict:
            self.assertEqual(len(split_dict[entry]['LubyIndex']), 16)
            self.assertEqual(len(split_dict[entry]['DropletMessage']), 256)
            self.assertEqual(len(split_dict[entry]['ErrorCorrection']), 16)
    def test_reverse_luby(self):
        sequences, blocks = load_in_files()
        blocks = split_luby(blocks)
        binary_dict = convert_seq_to_binary(sequences)
        split_dict = split_fasta(binary_dict)
        luby_reversed_message = reverse_luby(split_dict, blocks)
        self.assertIsInstance(luby_reversed_message, str)
    def test_convert_to_ascii(self):
        binary_string = "01000001011100110010000001001001001000000111011101100001011011000110101100100000011101000110100001110010011011110111010101100111011010000010000001110100011010000110010100100000011101100110000101101100011011000110010101111001001000000110111101100110001000000111010001101000011001010010000001110011011010000110000101100100011011110111011100100000011011110110011000100000011001000110010101100001011101000110100000100000010010010010000001110100011000010110101101100101001000000110000100100000011011000110111101101111011010110010000001100001011101000010000001101101011110010010000001101100011010010110011001100101001000000110000101101110011001000010000001110010011001010110000101101100011010010111101001100101001000000111010001101000011001010111001001100101001001110111001100100000011011100110111101110100011010000110100101101110001001110010000001101100011001010110011001110100"
        string = convert_to_ascii(binary_string)
        self.assertTrue(contains_no_integers(string))
        self.assertEqual(string, "As I walk through the valley of the shadow of death I take a look at my life and realize there's nothin' left")
    def test_bitwise_xor(self):
        binary_1 = "101"
        binary_2 = "110"
        binary_3 = "001"
        output1 = bitwise_xor(binary_1, binary_2) #011
        output2 = bitwise_xor(binary_2, binary_3) #111
        output3 = bitwise_xor(binary_3, binary_1) #100
        output4 = bitwise_xor(binary_2, binary_1) #011
        output5 = bitwise_xor(output1, binary_3)  #010
        self.assertEqual(output1, output4)
        self.assertTrue(output1 == '011')
        self.assertTrue(output2 == '111')
        self.assertTrue(output3 == '100')
        self.assertTrue(output5 == '010')

    def test_convert_str_to_binary(self):
        nucleotide_string = 'AAAAAAAAGACGGCTCACAAGTGAGCCAGCGGACAAGTCGGCGGGCAGGTACACAAATATATAAATGTATGAACTAACAAGCCAGTGGGCTGGCAGGCTCGCCGGTGAGTCGACAAGCCAGCAGGTATACAAGTGAGGCTATAT'
        binary = convert_str_to_binary(nucleotide_string)
        self.assertTrue(contains_no_letters(binary))
        nucleotide_string_two = 'AAGATCATTGGC'
        expected_binary_two = '000001001011001010010111'
        output_binary_two = convert_str_to_binary(nucleotide_string_two)
        self.assertEqual(output_binary_two, expected_binary_two)

    def test_split_binary_droplet_string(self):
        binary_string = '000000000000000001001101011110110011000001100100011111000111010100110000011011010111010101110001011000110011000000100010001000000010011000100100001110000011000001111100011001010111100101110001011110110111110101100100011011010011000001111100011100010110001000110000011001000101111000100010'
        binary = convert_str_to_binary(binary_string)
        _, _, message = split_binary_droplet_string(binary)
        self.assertEqual(len(message), 256)
        self.assertEqual(message, '0100110101111011001100000110010001111100011101010011000001101101011101010111000101100011001100000010001000100000001001100010010000111000001100000111110001100101011110010111000101111011011111010110010001101101001100000111110001110001011000100011000001100100')

    def test_block_0(self):
        string = 'AAAAAAAAGACGGCTCACAAGTGAGCCAGCGGACAAGTCGGCGGGCAGGTACACAAATATATAAATGTATGAACTAACAAGCCAGTGGGCTGGCAGGCTCGCCGGTGAGTCGACAAGCCAGCAGGTATACAAGTGAGGCTATAT'
        binary = convert_str_to_binary(string)
        _, _, message = split_binary_droplet_string(binary)
        ascii_string = convert_to_ascii(message)
        print(ascii_string)

if __name__ == '__main__':
    unittest.main()

