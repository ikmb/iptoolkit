#!/usr/bin/env Python 
"""
@author: Hesham ElAbd
@brief: The module contains functions to parse & call meme software from the command line  
@version: 0.0.1
@date: 25.08.2020  
"""
# load the module 
from Bio import motifs 
import subprocess as sp 
import os
from IPTK.Utils.UtilityFunction import append_to_calling_string
# DEFINE SOME CONSTANTS 
MEME_IS_NOT_INSTALLED_ERROR=f"meme is either not installed or not part of the PATH ==> {os.environ['PATH']}"

## define the interface functions 
def is_meme_callable()->bool:
    """
    @brief return wether or not meme is callable via a system call or not 
    """
    try: 
        res=sp.run(['meme'],stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        return True
    except FileNotFoundError as exp: 
        return False

def get_meme_help()->None:
    """
    @brief: print the command line help interface for the meme tool
    """
    if is_meme_callable():
        sp.run(['meme','-h'])
    else: 
        raise FileNotFoundError(MEME_IS_NOT_INSTALLED_ERROR)



def call_meme(input_fasta_file:str, output_dir:str, verbose: bool=True, objfunc: str='classic', test: str = 'mhg',
              use_llr: bool = False, shuf: int = 2, hsfrac: float = 0.5, cefrac: float = 0.25, 
              searchsize: int= -1, maxsize: int = -1 , norad: bool = False, csites: int = -1, seed: int = -1,
              mod: str='oops', nmotifs: int= -1, evt: float = -1.0, time: int = -1, nsite: int = -1, 
              minsites: int = -1, maxsite: int = -1, wnsites: int = -1, w: int = -1, minw: int = -1, 
              maxw: int = -1, nomatrim: bool = False, wg: int = -1, ws: int = -1, noendgaps: bool = False,
              maxiter: int = -1, prior: str = 'dirichlet', b: int = -1, p: int = -1 )->None:
    """
    @brief: a wraper for making a system call to meme software for sequence motif finding
    @param: input_fasta_file: The path to input FASTA files 
    @param: output_dir: the output dir to write the results, **IT WILL OVERWRITE EXISTING DIRECTORY** 
    @param: verbose: whether or not to print the output of calling meme to the screen, default is True. 
    @note: for the reset of the function parameters use the function **get_meme_help** defined in the module 
    IO, submodule MEMEInterface 
    """
    calling_string=""
    ## adding the parameters to the calling string based upon their values 
    calling_string=append_to_calling_string(str(objfunc), 'classic', objfunc,calling_string)
    calling_string=append_to_calling_string(str(test), 'mhg', test, calling_string)
    calling_string=append_to_calling_string(str(use_llr), False, use_llr, calling_string)
    calling_string=append_to_calling_string(str(shuf), 2, shuf, calling_string)
    calling_string=append_to_calling_string(str(hsfrac), 0.5, hsfrac, calling_string)
    calling_string=append_to_calling_string(str(cefrac), 0.25, cefrac, calling_string)
    calling_string=append_to_calling_string(str(searchsize), -1, searchsize, calling_string)
    calling_string=append_to_calling_string(str(maxsize), -1, maxsize, calling_string)
    calling_string=append_to_calling_string(str(norad), False, norad, calling_string)
    calling_string=append_to_calling_string(str(csites), -1, csites, calling_string)
    calling_string=append_to_calling_string(str(seed), -1, seed, calling_string)
    calling_string=append_to_calling_string(str(mod), 'oops', mod, calling_string)
    calling_string=append_to_calling_string(str(nmotifs), -1, nmotifs, calling_string)
    calling_string=append_to_calling_string(str(evt), -1, evt, calling_string)
    calling_string=append_to_calling_string(str(time), -1, time, calling_string)
    calling_string=append_to_calling_string(str(nsite), -1, nsite, calling_string)
    calling_string=append_to_calling_string(str(minsites), -1, minsites, calling_string)
    calling_string=append_to_calling_string(str(maxsite), -1, maxsite, calling_string)
    calling_string=append_to_calling_string(str(wnsites), -1, wnsites, calling_string)
    calling_string=append_to_calling_string(str(minw), -1, minw, calling_string)
    calling_string=append_to_calling_string(str(maxw), -1, maxw, calling_string)
    calling_string=append_to_calling_string(str(nomatrim), False, nomatrim, calling_string)
    calling_string=append_to_calling_string(str(wg), -1, wg, calling_string)
    calling_string=append_to_calling_string(str(ws), -1, ws, calling_string)
    calling_string=append_to_calling_string(str(maxiter), -1, maxiter, calling_string)
    calling_string=append_to_calling_string(str(prior), 'dirichlet', prior, calling_string)
    calling_string=append_to_calling_string(str(b), -1, b, calling_string)
    calling_string=append_to_calling_string(str(p), -1, p, calling_string)
    # add the output path
    calling_string+=" -oc "+str(output_dir)
    ## call the program 
    # generate the calling command 
    calling_command=['meme']
    calling_command.extend(calling_string.split(' ')) 
    # call the porgram 
    try:
        if verbose:
            sp.run(calling_command, check=True)
        else:
            sp.run(calling_command, check=True, stderr=sp.DEVNULL, stdout=sp.DEVNULL)
    except sp.CalledProcessError as exp: 
        raise RuntimeError(f'While calling meme to compute the motif, the following was discovered: {exp}')
    
def parse_meme_results(path_to_meme_file):
    """
    @brief: Use BioPython to parse the XML file generate by biopython 
    @param: path_to_results: the path to the XML file
    """
    try: 
        with open(path_to_meme_file,'r') as input_buf: 
            res=motifs.parse(input_buf,'meme')
    except Exception as exp: 
        raise RuntimeError(f'While parsing the file: {path_to_meme_file} The following error was encountered: {exp}')
    # return the results
    return  res 
    



    

    

    
    
    
    
    
    
    
    
    
    
    
    
    

    


    