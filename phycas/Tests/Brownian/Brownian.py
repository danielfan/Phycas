from phycas import *

if __name__ == '__main__':
    phycas = Phycas()
    phycas.log_file_name = 'logfile.txt'
    phycas.brownian_input_tree_file = 'poecilia2.tre'
    phycas.brownian()
    