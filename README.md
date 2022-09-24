# SCaFoSpy
..Description coming

## Preparing files
Input directory should contain any alignment files. Usually .ali or .fasta.
Sequence names of one sample in alignment files should either be named exactly the same 
or the deviating part should be marked by '@': 'sample@specification' (e.g. Lepidosiren_paradoxa@Lpar85533.E.lc and Lepidosiren_paradoxa@Lpar103035.E.lc)

## Running script

```python3 SCaFoS.py --dir <Path to input directory> --out <Path to output directory> --type <type of sequences> --mode <mode of alignment comparison>```

`--type` Defines type of sequences. Select between 'protein' and 'nucleotide'

`--mode` Defines type of alignment comparison. Select between 'ambiguity', 'longest_sequence' and 'lowest_diversion'

## Modes
The different modes define how the script evaluates positions in an alignment where different protein or nuclotide characters appear in the respective sequences. Thus the consensus sequence will look different at this position for each mode.

`ambiguity` Here in the consensus sequence '?' is put in place for ambiguous positions.

`longest_sequence` Here in the consensus sequence the char of the sequence with the most information (number of chars that are not missing data) is put in place for ambiguous positions.

`lowest_divergence` Here in the consensus sequence the most common char from ALL sequences of the alignment is put in place for ambiguous positions.

## Output
All files are put with their original name in a new folder now containing a consensus sequence for each sample.
