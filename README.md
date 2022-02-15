# dna-optimization

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

Reference
Expression host: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=10090
https://www.genscript.com/tools/codon-frequency-table

simplex https://github.com/mmolteratx/Simplex/blob/master/LinearModel_2-PhaseSimplex.ipynb

https://en.wikipedia.org/wiki/List_of_genetic_codes

1. Codon optimization
2. Splicing removal
3. GGGG and CCCC removal
4. More WCRH
5. More RCY
