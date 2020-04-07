# sequenceBinder
Tool for binding counterions onto charged protein residues and discussed in detail within the research paper below:  
[Model for Counterion Binding and Charge Reversal on Protein Surfaces](https://doi-org.manchester.idm.oclc.org/10.1021/acs.molpharmaceut.9b01047).

<p align="center">
<img src="https://github.com/jkalayan/sequenceBinder/blob/master/abstract_image5.png" width="400">
</p>

Example: 
```
python sequenceBinder.py --fileName file.pdb --pH 7 --conc_neg 0.1 --neg_charge 3 --pos_charge 1
``` 
**Note:** Only set the concentration of the positive *or* negative counterion. The concentration of the other ion is calculated within the tool so that concentrations are charge balanced.

```
optional arguments:

  -h, --help            show this help message and exit
  
  -f FILENAME, --fileName FILENAME
                        input pdb list file name (default: pdbs.list)
                        
  -pH PH, --pH PH       input pH (default: 7)
  
  -cn CONC_NEG, --conc_neg CONC_NEG
                        input concentration (Molar) of negative ion (default:
                        nul)
                        
  -cp CONC_POS, --conc_pos CONC_POS
                        input concentration (Molar) of positive ion (default:
                        nul)
                        
  -n NEG_CHARGE, --neg_charge NEG_CHARGE
                        charge of negative ion (default: 1)
                        
  -p POS_CHARGE, --pos_charge POS_CHARGE
                        charge of positive ion (default: 1)
                        
  -r, --res_counts      include flag if residue counts are required (default:
                        False)
```

Output files:
1. .csv file containing binding information such as protein start charge, end charge, etc.
2. If the -r flag is used, then a .csv containing number of each residue type in protein is also outputted.
