# SCaFoSpy
..Description coming

## Preparing files
Input directory should contain any alignment files. Usually .ali or .fasta.
Sequence names of one sample in alignment files should either be named exactly the same 
or the deviating part should be marked by '@': 'sample@specification' (e.g. Lepidosiren_paradoxa@Lpar85533.E.lc and Lepidosiren_paradoxa@Lpar103035.E.lc)

## Running script

```python3 SCaFoS.py --dir <Path to input directory> --out <Path to output directory> --type <type of sequences> ```

`--type` Defines type of sequences. Select between 'protein' and 'nucleotide'

## Output
All files are put with their original name in a new folder now containing a consensus sequence for each sample.
