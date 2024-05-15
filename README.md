# Challenge One
## Advanced Synthetic and Systems Biology
By: Olaf

### Code Structure
##### Input
*  `load_in_files()` 
Files are read in via the `load_in_files()` function. It saves the fasta file as a dictionary and the csv as a dataframe. The function does not current take custom file names and does not contain error handling. It depends on os to locate the input_files director, SeqIO to read the fasta file, and pandas to read the csv. 
##### Processing
* `process_files`
  Wrapper function for the subsequent processing functions. 
* `split_luby`
Uses regular expressions to extract patterns from the pandas dataframe. Lambda functions are used to apply this pattern across the dataframe and split the dataframe from one column to 4.
Columns are named: 'OriginalValue', 'DropletNumber', and 'BlockNumbers'. BlockNumbers is isolated in two steps to create a list within the dataframe column that does not have extra '[' or ']'
* `convert to binary`
Takes the droplet sequence dict (generated from the fasta file) and replaces nucleotide sequence with binary strings. It relies on a forloop and the following encoding scheme `{"A": '00', "G": '01', "T": '10', "C": '11', }`
* `split_fasta`
Takes the binary dictionary and renames the keys to just the number. Each key's value is replaced with a list containing 'LubyIndex', 'ErrorCorrection', and 'DropletMessage.' This split is created by predetermined slices 0:16, 272:end, and 16:272 respectively.
##### Reverse Luby
* `reverse_luby`
Identifies droplets containing only one block. Then loops through the remainder of the droplets solving droplets with only one unknown block by calls the `bitwise_xor` function to calculate the bitwise XOR value. The droplets are sorted by number of block components to reduce looping time. The resulting strings are concatenated and returned. 
* `bitwise_xor`
Takes two strings and iterates through each character. It performs a XOR operation on the respective digits. It obeys the following truth table:
  | Input1         | Input2         | Output                      |
| :----------- | :--------------: | -------------------------: |
| 0| 0 | 0  |
| 0    | 1   | 1 |
| 1| 0 | 1  |
| 1    | 1   | 0 |
##### Convert to ASCII





