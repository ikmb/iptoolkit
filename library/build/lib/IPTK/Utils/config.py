""""
A configuration wrapper most of the library parameters 
"""
from __future__ import annotations
from typing import Dict,Set,List,Union
import matplotlib as mpl 
class ParamConfig:
    def __init__(self,output_directory:str=None,
        matplot_lib_param:Dict[str,Union[int,float,str]]=None, 
        plotly_lib_param:Dict[str,Union[int,float,str]]=None)-> ParamConfig:
        """Create a new parameter configuration, which library global parameters 

        Args:
            matplot_lib_param (Dict[str,Union[int,float,str]], optional): a dict of parameters that controls global matplotlib parameters. Defaults to None.
        Returns:
            ParamConfig: [description]
        """
        self.params["matplot_lib_param"]=matplot_lib_param,
        return

    def set(self)->None:
        if self.params["matplot_lib_param"]  is not None:
            mpl.rcParams.update(self.params["matplot_lib_param"])


    def update(self,param_dict:Dict[str,Union[Dict[str,Union[int,float,str]],int,float,str]])->None:
        """Takes a dict to update the parameters of the config dict 

        Args:
            param_dict (Dict[str,Union[Dict[str,Union[int,float,str]],int,float,str]]): the new configuration parameters 
        """
        self.params.update(param_dict) 
