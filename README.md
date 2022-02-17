# dna-optimization
## Overview
This small command line interface tool allows you to run the following DNA optimizations on a given DNA sequence:  
| Type  | Description | Development Stage |
| ------------- | ------------- | ------------- |
| Codon optimization  | TBC | WIP |
| Splicing removal  | Identify splicing sites using [this tool](https://www.fruitfly.org/cgi-bin/seq_tools/splice.pl) and remove them  | Supported |
| GGGG and CCCC removal | Identify GGGG and CCCC codons and replace them | Supported |
| Increase WCRH/RCY | TBC | WIP |

The optimizations preserve the amino acids representation of the sequence and only optimize the codon selection.

## Next Phases
- [ ] Support codon optimization and increase WCRH/RCY rules.
- [ ] Support optimization of a sub-sequence of the input (based on indices).
- [ ] Support multiple, sequential, optimizations. For example, run GGGG removal, then CCCC removal, etc.
- [ ] Support composition of rules with sophisticated optimization algorithm (simplex?).
- [ ] Order logs and add colors.
- [ ] Add validations of input.  
- [ ] Support optimization of multiple DNA sequences at once.
- [ ] Support multiple expression host.

## Integration with Replit

### Run 
Press on "Run" button.  

### Get Updates
In the left-hand bar choose "Version control" and press "Pull".

## Local Development Environment

### Python Virtual Environment
```shell
# Create venv (only first time)
python3 -m venv venv
# Activate venv
source venv/bin/activate
# Install Requirements (only first time)
python3 -m pip install -r requirements.txt
# Deactivate venv
deactivate
```

## References

### Biology
* Expression host: [link](https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=10090)
* Frequency Table: [link](https://www.genscript.com/tools/codon-frequency-table)

### Mathematics
* Simplex: [link](https://github.com/mmolteratx/Simplex/blob/master/LinearModel_2-PhaseSimplex.ipynb)

