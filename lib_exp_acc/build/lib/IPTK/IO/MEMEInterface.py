#!/usr/bin/env Python 
"""The module contains functions to to call meme software via a system call.   
"""
# load the module 
from Bio import motifs 
import subprocess as sp 
import os
from IPTK.Utils.UtilityFunction import append_to_calling_string
## define the interface functions 
def is_meme_callable()->bool:
    """
    :return: True if meme is callable, False otherwise. 
    :rtype: bool
    """
    try: 
        res=sp.run(['meme'],stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        return True
    except FileNotFoundError: 
        return False

def get_meme_help()->None:
    """Print the command line help interface for the meme tool

    :raises FileNotFoundError: if meme is not callable 
    """
    if is_meme_callable():
        sp.run(['meme','-h'])
    else: 
        raise FileNotFoundError(f"meme is either not installed or not part of the PATH ==> {os.environ['PATH']}")

def call_meme(input_fasta_file:str, output_dir:str, verbose: bool=True, objfunc: str='classic', test: str = 'mhg',
              use_llr: bool = False, shuf: int = 2, hsfrac: float = 0.5, cefrac: float = 0.25, 
              searchsize: int= -1, maxsize: int = -1 , norand: bool = False, csites: int = -1, seed: int = -1,
              mod: str='oops', nmotifs: int= -1, evt: float = -1.0, time: int = -1, nsite: int = -1, 
              minsites: int = -1, maxsite: int = -1, nsites: int = -1, w: int = -1, minw: int = -1, 
              maxw: int = -1, nomatrim: bool = False, wg: int = -1, ws: int = -1, noendgaps: bool = False,
              maxiter: int = -1, prior: str = 'dirichlet', b: int = -1, p: int = -1 )->None:
    """warper for making a system call to meme software for sequence motif finding \
        for the reset of the function parameters use the function \
            **get_meme_help** defined in the module IO, submodule MEMEInterface. \ 
            
    :param input_fasta_file: The path to input FASTA files. 
    :type input_fasta_file: str
    :param output_dir: the output dir to write the results, **IT WILL OVERWRITE EXISTING DIRECTORY** 
    :type output_dir: str
    :param verbose:  whether or not to print the output of calling meme to the screen, default is True. 
    :type verbose: bool 
    """
    # define the calling string backbone
    calling_string=""
    # add the input dataset:
    calling_string += str(input_fasta_file)+" "
    ## adding the parameters to the calling string based upon their values 
    calling_string=append_to_calling_string('objfunc', 'classic', objfunc,calling_string,False)
    calling_string=append_to_calling_string('test', 'mhg', test, calling_string,False)
    calling_string=append_to_calling_string('use_llr', False, use_llr, calling_string,True)
    calling_string=append_to_calling_string('shuf', 2, shuf, calling_string,False)
    calling_string=append_to_calling_string('hsfrac', 0.5, hsfrac, calling_string,False)
    calling_string=append_to_calling_string('cefrac', 0.25, cefrac, calling_string,False)
    calling_string=append_to_calling_string('searchsize', -1, searchsize, calling_string,False)
    calling_string=append_to_calling_string('maxsize', -1, maxsize, calling_string,False)
    calling_string=append_to_calling_string('norand', False, norand, calling_string,True)
    calling_string=append_to_calling_string('csites', -1, csites, calling_string,True)
    calling_string=append_to_calling_string('seed', -1, seed, calling_string,False)
    calling_string=append_to_calling_string('mod', 'oops', mod, calling_string, False)
    calling_string=append_to_calling_string('nmotifs', -1, nmotifs, calling_string,False)
    calling_string=append_to_calling_string('evt', -1, evt, calling_string,False)
    calling_string=append_to_calling_string('time', -1, time, calling_string,False)
    calling_string=append_to_calling_string('nsite', -1, nsite, calling_string,False)
    calling_string=append_to_calling_string('minsites', -1, minsites, calling_string,False)
    calling_string=append_to_calling_string('maxsite', -1, maxsite, calling_string,False)
    calling_string=append_to_calling_string('nsites', -1, nsites, calling_string,False)
    calling_string=append_to_calling_string('minw', -1, minw, calling_string,False)
    calling_string=append_to_calling_string('maxw', -1, maxw, calling_string,False)
    calling_string=append_to_calling_string('nomatrim', False, nomatrim, calling_string,True)
    calling_string=append_to_calling_string('wg', -1, wg, calling_string,False)
    calling_string=append_to_calling_string('ws', -1, ws, calling_string,False)
    calling_string=append_to_calling_string('maxiter', -1, maxiter, calling_string,False)
    calling_string=append_to_calling_string('prior', 'dirichlet', prior, calling_string,False)
    calling_string=append_to_calling_string('noendgaps', False, noendgaps, calling_string,True)
    calling_string=append_to_calling_string('b', -1, b, calling_string,False)
    calling_string=append_to_calling_string('p', -1, p, calling_string,False)
    # add the output path
    calling_string+=" -oc "+str(output_dir)
    ## call the program 
    # generate the calling command 
    calling_command=['meme']
    calling_command.extend(calling_string.split(' ')) 
    # remove empty spaces from the calling string 
    calling_command=[elem for elem in calling_command if elem != '']
    # call the porgram 
    try:
        if verbose:
            sp.run(calling_command, check=True)
        else:
            sp.run(calling_command, check=True)
    except sp.CalledProcessError as exp: 
        raise RuntimeError(f'While calling meme to compute the motif, the following was encountered: {exp}')
    




    

    

    
    
    
    
    
    
    
    
    
    
    
    
    

    


    