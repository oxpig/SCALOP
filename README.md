# SCALOP

## Description
 
SCALOP (Python 3)
Sequence-based antibody Canonical LOoP structure annotation

## Installation

The easiest way to install SCALOP dependencies is using `conda`:

```python
git clone https://github.com/oxpig/SCALOP.git
conda create -n scalop-env python=3.8 -y
conda activate scalop-env
conda install -c bioconda numpy pandas hmmer biopython -y
pip install SCALOP/
```

**Input format**
* Single complete antibody sequence 
* Complete heavy and light chain sequences separated by '/'
* Fasta file 

**Example** Command line tool
```bash
SCALOP -i VKLLEQSGAEVKKPGASVKVSCKASGYSFTSYGLHWVRQAPGQRLEWMGWISAGTGNTKYSQKFRGRVTFTRDTSATTAYMGLSSLRPEDTAVYYCARDPYGGGKSEFDYWGQGTLVTVSS/ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIKRTV --scheme imgt --definition north
```

**Example** Python module
```python
from scalop.predict import assign
input='VKLLEQSGAEVKKPGASVKVSCKASGYSFTSYGLHWVRQAPGQRLEWMGWISAGTGNTKYSQKFRGRVTFTRDTSATTAYMGLSSLRPEDTAVYYCARDPYGGGKSEFDYWGQGTLVTVSS/ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIKRTV'
assign(input)
```

**Dependencies**        
* HMMER 3.1b2 (February 2015); http://hmmer.org/
* numpy v. 1.10+
* pandas
		
## Authors
 
2018
Wing Ki Wong
Dr Jinwoo Leem
Prof Charlotte M. Deane - Oxford Protein Informatics group.

**Contact** opig@stats.ox.ac.uk