This preprocessing pipeline requires a directory that includes sorted raw gene counts table in .csv format.

Once you run the './main.py', a Window GUI prompt asks for the directory to be selected.

If you have different set of DLP barcodes, you need to edit the 'process_data()' in './utils.py', especially 'lookuptable', 'drug_name', and '__get_lookup_table__()'. 
