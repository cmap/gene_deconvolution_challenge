lastemperor
==============================

Project with the 8 place model in CMAP competition

Project Organization
------------

    ├── README.md          <- The top-level README for developers using this project.
    ├── barcode_to_gene_map.txt <- Barcode for gens
    └── main.py            <- Script for data processing and predicting


--------
## Prerequisites:
Linux, Python 3.7, pandas

Manual setup:
```
pip install pandas
```

## Running
0. Put original .csv files
1. Run python main.py with arguments: --dspath PATH/TO/INPUT/DATASET/FOLDER --out /PATH/TO/OUTPUT/FOLDER  --create_subdir 0  --plate TEST_SET_NAME

## Results
TEST_SET_NAME.csv file is in /PATH/TO/OUTPUT/FOLDER folder