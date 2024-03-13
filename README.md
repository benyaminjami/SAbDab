**PDB data**

Please download all the structure data of antibodies from the [download page of SAbDab](http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/search/?all=true). Please enter the *Downloads* tab on the left of the web page and download the tsv summary file and archived zip file for the structures, then decompress it:

```bash
wget https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/summary/all/ -O summaries/sabdab_summary.tsv
wget https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/archive/all/ -O all_structures.zip
unzip all_structures.zip all_structures/imgt/* -d .
```

You should get a folder named *all_structures* with the following hierarchy:

```
├── all_structures
│   ├── imgt
```

The imgt subfolder contains the pdb files renumbered with the imgt scheme.

**Data**

To preprocess the raw data, we need to first generate summaries for each benchmark in json format, then split the datasets into train/validation/test sets, and finally transform the pdb data to python objects. We have provided the script for all these procedures in `scripts/data_preprocess.sh`. Suppose the IMGT-renumbered pdb data are located at `all_structures/imgt/`, and that you want to store the processed data (~5G) at `all_data`, you can simply run:

```bash
bash data_preprocess.sh all_structures/imgt all_data
```
which takes about 1 hour to process SAbDab, RAbD, Igfold test set, and SKEMPI V2.0. It is normal to see reported errors in this process because some antibody structures are wrongly annotated or have wrong format, which will be dropped out in the data cleaning phase.
