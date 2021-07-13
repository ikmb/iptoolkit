""""
A configuration wrapper most of the library parameters 
"""
from __future__ import __annotations__
from typing import Dict,Set,List,Union

class ParamConfig:
    def __init__(self,output_directory:str=None,
        matplot_lib_param:Dict[str,Union[int,float,str]]=None, 
        plotly_lib_param:Dict[str,Union[int,float,str]]=None)-> ParamConfig:
        """Create a new parameter configuration, which library global parameters 

        Args:
            output_directory (str, optional): A global output directory where all results will be written to. Defaults to None.
            matplot_lib_param (Dict[str,Union[int,float,str]], optional): a dict of parameters that controls global matplotlib parameters. Defaults to None.
            plotly_lib_param (Dict[str,Union[int,float,str]], optional): [description]. Defaults to None.

        Returns:
            ParamConfig: [description]
        """
        self.params["output_directory"]= output_directory,
        self.params["matplot_lib_param"]=matplot_lib_param,
        self.params["plotly_lib_param"]=plotly_lib_param
        return

    def set(self)->None:
        pass 

    def update(self,param_dict:Dict[str,Union[Dict[str,Union[int,float,str]],int,float,str]])->None:
        pass 
