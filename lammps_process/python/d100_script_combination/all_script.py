#!/usr/bin/env python3
import d02_intermediate.out_pickle as out_pickle
import d02_intermediate.out_npy as out_npy
import d02_intermediate.rename_variable as rename_variable
import d03_processing.out_new_variable as out_new_variable
import d03_processing.out_std as out_std
import d07_visualization.output_figures as output_figures
import d08_latex.output_latex as output_latex
import d09_sync.sync as sync

def all_script_without_sync():
    out_pickle.main()
    out_npy.main()
    rename_variable.main()
    out_new_variable.main()
    out_std.main()
    output_figures.main()
    output_latex.main()
    
def all_script():
    all_script_without_sync()
    sync.main()

def main():
    all_script()

if __name__ == "__main__":
    main()