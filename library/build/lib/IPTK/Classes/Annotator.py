#!/usr/bin/env python 
r"""The class provides methods for visualizing different aspects of the protein \
    biology. This is achieved through three main methods: 
        1- add_segmented_track: which visualize information about \
        non-overlapping protein substructures, for example, protein domains.

        2- add_stacked_track: which visualize information about \
        overlapping protein substructures, for example, splice variants. 
        
        3- add_marked_positions_track: which visualize or highlight \
        positions in the protein, for example, sequence variants, or PTM. 
    
    The class also provides functions for visualizing the relationship between 
    a protein and its eluted peptide/peptides in an analogous manner to the way 
    NGS reads are aligned to genomic regions. This can be useful to identify  
    regions in the protein with high/low number of eluted peptides, i.e.,Coverage.
    Also, to link it with other facests of the protein like domain organization,\
    PTM, sequence/splice variants. 
    
    Notes
    -----
    each figure should have a base track this can be done explicitly by calling
    the function add_base_track or by implicitly by calling the function 
    add_coverage_plot with the parameter coverage_as_base=True.

"""
# load the modules 
from __future__ import annotations
import numpy as np 
import matplotlib.pyplot as plt 
import os 
from typing import Tuple, List, Dict, Union
# define the class 
class Annotator:
    r""" A high level API to plot information about the protein, for example, PTM, Splice variant etc, using matplotlib library 
    """
    def __init__(self,protein_length: int,figure_size:Tuple[int,int],figure_dpi:int,
                 face_color="white")->Annotator:
        """
        Parameters
        ----------
        protein_length : int 
            the length of the protein.
        figure_size : tuple/list 
            The figure size in inches. The first element in the tuple/list is 
            the width and the second is the height.
        figure_dpi : int
            the dots per inch, which control the figure resolution.
        face_color : string, optional
            the bachground colors . The default is "white".
        Returns
        -------
        Annotator.

        """
        #check the user input
        if protein_length<0:
            raise ValueError("Protein length must be bigger than zero")

        if len(figure_size)!=2: 
            raise ValueError("figure_size must be a tuple or list with two elements")
        # create the figure 
        self.fig=plt.figure(figsize=figure_size,dpi=figure_dpi,
                            facecolor=face_color) 
        # add the protein length to self
        self.protein_length=protein_length
        self.figure_dpi=figure_dpi
        self.figure_size=figure_size
        return
    
    def add_base_track(self, space_fraction:float=0.3,
                       protein_name_position:float=0.5,
                       track_label:str="base_track",
                       track_label_dict: Dict[str,Union[int,str]]={"fontsize":8,"color":"black"},
                       protein_name: str="A protein", 
                       protein_name_dict:Dict[str,Union[int,str]] ={"fontsize":10,"color":"black"},
                       rect_dict: Dict[str,Union[int,str]]={"color":"olive","capstyle":"butt"},
                       number_ticks: int =10,
                       xticks_font_size: int =4,
                          ):
        """
        Description
        -----------
        Adds a base track to the figure. 

        Parameters
        ----------
        space_fraction : float, optional
            A float between 0 and 1 that represent the fraction of space left 
            below and above the track. The default is 0.3 which means that the
            track will be drown on a 40% while 60% are left as an empty space
            below and above the track.
            
        protein_name_position : float, optional
            A float between 0 and 1 which control the relative position of the 
            protein name on the y-axis. The default is 0.5.
            
        track_label : string, optional
            The name on the track, which will be shown on the y-axis. 
            The default is "base_track".
            
        track_label_dict : Dict[str,Union[int,str]], optional
            The parameters that control the printing of the track_label,
            for example, the font size and the color. These parameters should 
            be provided as dict that will be fed to the function 
            ``axes.set_ylabel``.The default is {"fontsize":8,"color":"black"}.
            
        protein_name : string, optional
            The name of the protein to be printed to the track.
            The default is "A protein".
            
        protein_name_dict : Dict[str,Union[int,str]], optional
            the parameters that control the printing of the protein name,
            for example, the font size and the color. These parameters should 
            be provided as dict that will be fed to the function 
            ``axes.text()``. The default is {"fontsize":10,"color":"black"}.
            
        rect_dict : Dict[str,Union[int,str]], optional
            a dictionary that control the character of the track itself, for example,
            the color and the transparency. this dict will be fed to the function
            ``plt.Rectangle()``. 
            The default is {"color":"olive","capstyle":"butt"}.
        
        number_ticks:int
            The number of ticks on the x-axis. The default is 10.
        
        xticks_font_size: int
            The font size of the x-axis ticks. The default is 4.
            
        Returns
        -------
        None.
        
        Examples
        --------
        >>> example_1=VisTool(250,(3,5),300) 
            # create a graph of size 3 inches by 5 inches with a 300 dots per
            # inch (DPI) as a resolution metric for a protein of length 250 amino acids
        
        >>> example_1.add_base_track() 
            # adds a basic track using the default parameters. 
        
        >>> example_1.add_base_track(space_fraction=0.1,
                                    track_label="example_1",
                                    track_label_dict={"fontsize":5,"color":"blue"}
                                    number_ticks=5,
                                    xticks_font_size=6)
            # generate a base track with 10% empty space above and below 
            #  the track. Track will have the name example_1 and it will be 
            # shown in font 5 instead of 8 and in blue color instead of black.
            # five ticks will be shown on the x-axis using a font of size 6.
        
        Notes
        -----
            calling the function more than once will result in an overriding of
            the previously added base track, for example, in the examples 
            section calling add_base_track for the second time will overrides 
            the graph build by the previous call. 

        """
        # add a subplot to the figure
        ax=self.fig.add_subplot(111)
        # set the axis limits for the axes
        ax.axis(xmin=0,xmax=self.protein_length,ymin=0,
                ymax=1)
        # build the protein representation
        rectangle=plt.Rectangle((0,space_fraction),
            width=self.protein_length,height=1-(2*space_fraction),**rect_dict)
        # add it to the axes
        ax.add_patch(rectangle)
        # add the protein name 
        ax.text(x=self.protein_length/3,
                y=protein_name_position,s=protein_name,**protein_name_dict)
        # remove the spines of the axes
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        # flip the x axes
        ax.xaxis.tick_top()
        # adjust the ticks of the x-axis
        step_size=int(np.floor(self.protein_length/number_ticks))
        ax.set_xticks(ticks=np.arange(0,self.protein_length,step_size),minor=False)
        ax.tick_params(axis='x', which='major', labelsize=xticks_font_size)
        # remove the y-label ticks
        ax.get_yaxis().set_ticks([])
        # set the y-label
        ax.set_ylabel(track_label,**track_label_dict)
    
    def add_segmented_track(self,track_dict:Dict[str,Dict[str,Union[int,str]]],
                            track_label: str ="A segmented Track",
                            track_label_dict: Dict[str,Union[int,str]]={"fontsize":8,"color":"black"},
                            track_element_names_dict:Dict[str,Union[int,str]]={"fontsize":8,"color":"black"},
                            center_line_dict:Dict[str,Union[int,str,float]]={"alpha":0.5,"linewidth":0.5},
                            track_elements_dict: Dict[str,Union[int,str]]={"color":"brown","capstyle":"butt"},
                            show_names: bool=True)->None:
        """Adds a segmentation track which show non-overlapping features of the protein. 
        
        Parameters
        ----------
        track_dict : Dict[str,Dict[str,Union[int,str]]]
            A dict that contain the non-overlapping features of the protein. 
            The dict is assumed to have the following structure: a dict with the
            feature_index as a key and associated features as values. The associated
            features is a dict with the following three keys:
                
                1- Name: which contain the feature name

                2- startIdx: which contain the start position of the protein
                
                3- endIdx: which contain the end position of the protein
                
        track_label : str, optional
            The name of the track, which will be shown on the y-axis. 
            The default is "A segmented Track". 
            
        track_label_dict : Dict[str,Union[int,str]], optional
            the parameters that control the printing of the track_label,
            for example, the font size and the color. These parameters should 
            be provided as dict that will be fed to the function 
            ``axes.set_ylabel``. The default is {"fontsize":8,"color":"black"}.
            
        track_element_names_dict : Dict[str,Union[int,str]], optional
            the parameters that control the printing of the feature names on the track,
            for example, the font size and the color. These parameters should 
            be provided as a dict that will be fed to the function ``axes.text``. 
            The default is {"fontsize":8,"color":"black"}.

        center_line_dict: Dict[str,Union[int,str, float]], optional 
            The parameters that control the printing of the center line of a segmented track object.
            The default is {"fontsize":8,"color":"black"}.
            
        track_elements_dict : Dict[str,Union[int,str]], optional
            the parameters that control the printing of the feature 
            rectangluar representation for example the color, the dict will be fed 
            to the function ``plt.Rectangle``. 
            The default is {"color":"brown","capstyle":"butt"}.
        
        show_names : bool, optional
            whether or not to show the name of the features. The default is True.

        Returns
        -------
        None.
        
        Examples
        --------
        >>> test_dict={"domain1":{"Name":"domain_one","startIdx":55,"endIdx":150},
                       "domain2":{"Name":"domain_Two","startIdx":190,"endIdx":225}}
        # first define a dict object that define some protein features.

        
        >>> example_1=Annotator(protein_length=250, figure_size=(5,3), figure_dpi=200)
        # creating a Annotator instance
        
        >>> example_1.add_base_track()
        # add a base_track

        
        >>> example_1.add_segmented_track(test_dict) # build a segmented track using the default parameters
        # add the segmented track

       
        >>> example_1.add_segmented_track(track_dict=test_dict,
                                          track_label="Domains",
                                          track_elements_dict={"color":"brown"})
        # add a second segmented track with track name set to Domains and elements 
        # of the track shown as brown rectangles.
        
        Notes
        -----
        Any panel can have one or more segmented-tracks. Thus, in the above examples
        calling the method ``add_segmented_track`` for the second time does NOT
        override the previous segmented track it create a new one and added to the
        figure.
        
        """
        # change the geometry of the graph to adapt to the added axes
        number_axes=len(self.fig.axes)
        for i in range(number_axes):
            self.fig.axes[i].change_geometry(number_axes+1,1,i+1)
        ax=self.fig.add_subplot(number_axes+1,1,number_axes+1)
        # set the axes axis limits
        height=1
        ax.axis(xmin=0,xmax=self.protein_length,ymin=0,ymax=height)
        ax.axhline(y=0.45*height,xmin=0,xmax=self.protein_length,color="black", **center_line_dict)
        # get the position of the dictionary objects
        for key in track_dict.keys():
            segment_length=track_dict[key]["endIdx"]-track_dict[key]["startIdx"]
            ax.add_patch(plt.Rectangle((track_dict[key]["startIdx"],0.33*height), 
                width=segment_length, height=0.33*height,**track_elements_dict))
            if show_names:
                ax.text(
                    x=track_dict[key]["startIdx"],
                    y=0.5*height,s=track_dict[key]["Name"],
                    **track_element_names_dict)
        
        # adjust the figure of the track
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        # remove the ticks from the x and y axis
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        # add the track name to the libray
        ax.set_ylabel(track_label,**track_label_dict)
        # set the position of the trach

    def add_stacked_track(self,track_dict: Dict[str,Dict[str,Union[int,str]]],
                            track_label: str ="A stacked Track",
                            track_label_dict: Dict[str,Union[int,str]] ={"fontsize":8,"color":"black"},
                            track_element_names_dict: Dict[str,Union[int,str]]={"fontsize":8,"color":"black"},
                            track_elements_dict: Dict[str,Union[int,str]]={"color":"magenta","capstyle":"butt"},
                            base_line_dict: Dict[str,Union[int,str]]={"color":"black","linewidth":1},
                            show_names:bool =True):
        """
        Description
        -----------
        The function adds a stacked_track to a visualization panel. 
        The stacked track is used to show overlapping protein features, for example,
        different splice variants. 

        Parameters
        ----------
        track_dict : Dict[str,Dict[str,Union[int,str]]]
            A dict that contain the overlapping features of the protein. 
            The dict is assumed to have the following structure, a dict with the
            feature_index as a key and associated features as values. The associated
            features is a dict with the following three keys:

                1- Name: which contain the feature's name
                
                2- startIdx: which contain the start position of the feature.
                
                3- endIdx: which contain the end position of the feature.
                
       
        track_label : str, optional
            The name of the track, which will be shown on the y-axis. 
            The default is "A stacked Track".
            
        track_label_dict : Dict[str,Union[int,str]], optional
            the parameters that control the printing of the track_label,
            for example, the font size and the color. These parameters should 
            be provided as dict that will be fed to the function 
            ``axes.set_ylabel``.The default is {"fontsize":8,"color":"black"}.
            
        track_element_names_dict : Dict[str,Union[int,str]], optional
            the parameters that control the printing of the feature names on the track,
            for example, the font size and the color. These parameters should 
            be provided as a dict that will be fed to the function ``axes.text``. 
            The default is {"fontsize":8,"color":"black"}.
            
        track_elements_dict : Dict[str,Union[int,str]], optional
            the parameters that control the printing of the feature 
            rectangluar representation for example the color, the dict will be fed 
            to the function ``plt.Rectangle``.
            The default is {"color":"magenta","capstyle":"butt"}.
            
        base_line_dict : Dict[str,Union[int,str]], optional
            the parameters that control the shape of the base line,
            for example, color and/or line width. These parameters are going to be fed
            to the function ``axes.hlines``. 
            The default is {"color":"black","linewidth":1}.
            
        show_names : bool, optional
            whether or not to show the name of the features. The default is True.

        Returns
        -------
        None.
        
        Examples
        --------
        >>> test_dict={"feature_1":{"Name":"X","startIdx":55,"endIdx":150},
                       "feature_2":{"Name":"Y","startIdx":85,"endIdx":225},
                       "feature_3":{"Name":"Z","startIdx":160,"endIdx":240}}
         # first define a dict object that define some protein features.
        
        
        >>> example_1=Annotator(protein_length=250, figure_size=(5,3), figure_dpi=200)
        # creating a Annotator instance

        
        >>> example_1.add_base_track()
        # add a base_track

        
        >>> example_1.add_segmented_track(test_dict) # build a stacked track using the default parameters.
        # add the stacked track
       
        >>> example_1.add_segmented_track(track_dict=test_dict,
                                          track_label="OverLappingFeat",
                                          track_elements_dict={"color":"red"})
        # add a second segmented track with track name set to OverLappingFeat and elements 
        # of the track shown as red rectangles.
        
        Notes
        -----
        Any panel can have zero, one or more than one stacked-track.
        Thus, in the above examples calling the method ``add_stacked_track`` 
        for the second time does NOT override the previous stacked track 
        it creates a new one and added to the figure.
        """
        # change the geometry to accommodate a new figure
        number_axes=len(self.fig.axes)
        for i in range(number_axes):
            self.fig.axes[i].change_geometry(number_axes+1,1,i+1)
        # set the base height
        height=1
        # create a new axes:
        ax=self.fig.add_subplot(number_axes+1,1,number_axes+1)
        # adjust the x and y axis limits
        ax.axis(xmin=0,xmax=self.protein_length,ymin=0,ymax=height)
        # add a base line
        ax.hlines(y=0,xmin=0,xmax=self.protein_length,**base_line_dict)
        # adjust the space between the stacked and the formating characters
        number_of_stacked_tracked=len(track_dict)
        hight_per_track=height/number_of_stacked_tracked
        space_margine=0.2*hight_per_track # 0.2 below and above the graph will 
        track_hight=0.6*hight_per_track
        text_hight=0.5*hight_per_track
        track_elements=list(track_dict.keys())
        # plotting the tracks
        for idx in range(len(track_elements)):
            start_index=track_dict[track_elements[idx]]["startIdx"]
            end_index=track_dict[track_elements[idx]]["endIdx"]
            width=end_index-start_index
            name=track_dict[track_elements[idx]]["Name"]
            ax.add_patch(plt.Rectangle(
                (start_index,(idx*hight_per_track)+space_margine),
                width=width,height=track_hight,**track_elements_dict
                ))
            if show_names:
                ax.text(x=start_index,y=(idx*hight_per_track)+text_hight,
                    s=name,**track_element_names_dict)
        # processes the axes:
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.get_yaxis().set_ticks([])
        ax.get_xaxis().set_ticks([])
        ax.set_ylabel(track_label,**track_label_dict)
    
    def add_marked_positions_track(self,positions:List[int],
                    height_frac: float =0.5,
                    marker_bar_dict: Dict[str,Union[int,str]] ={"color":"black","linestyles":"solid"},
                    marker_dict: Dict[str,Union[int,str]] ={"color":"red","s":3},
                    track_label: str ="A marked positions Track",
                    track_label_dict: Dict[str,Union[int,str]] ={"color":"black","fontsize":8},
                    base_line_dict: Dict[str,Union[int,str]] ={"color":"black","linewidth":1}):
        """
        Description
        -----------
        The function adds a marked position to the track which is shown to highlight
        certain amino acid position within the protein, for example, a sequence
        variant position, or PTM position.

        Parameters
        ----------
        positions : List[int]
            a list that contain the position/positions that should be heighlighted
            in the protein sequence.
        
        height_frac : float
            the relative hight of the marked positions. The default is 0.5 which
            means that the hight of the marker will be 50% of the y-axis height.
            
        marker_bar_dict : Dict[str,Union[int,str]], optional
            The parameters of the marker position bar, for example, line width
            or color. These parameters are going to be fed to the function
            ``plt.hlines``. The default is {"color":"black","linestyles":"solid"}.
            
        marker_dict : Dict[str,Union[int,str]], optional
            These are the parameters for the marker points which sits on top of the marker bar,
            for example, the color, the shape or the size. 
            The default is {"color":"red","s":3}.
            
        track_label : str, optional
             The name of the track, which will be shown on the y-axis.
             The default is "A marked positions Track".
            
        track_label_dict : Dict[str,Union[int,str]], optional
            The parameters that control the printing of the track_label,
            for example, the font size and the color. These parameters should 
            be provided as dict that will be fed to the function 
            ``axes.set_ylabel``.The default is {"fontsize":8,"color":"black"}.
            
       base_line_dict : Dict[str,Union[int,str]], optional
            The parameters that control the shape of the base line, for example, color and/or line width. These parameters are going to be fed
            to the function ``axes.hlines``. 
            The default is {"color":"black","linewidth":1}.

        Returns
        -------
        None.
        
        Examples
        --------
        
        >>> test_list=[24,26,75,124,220]
        # first define a dict object that define some protein features.

       
        >>> example_1=Annotator(protein_length=250, figure_size=(5,3), figure_dpi=200)
        # creating a VisTool instance

       
        >>> example_1.add_base_track()
        # add a base_track
        
        >>> example_1.add_marked_positions_track(test_list) # build a marked position track using the default parameters
        # marked positions track

        
        >>> example_1.add_marked_positions_track(positions=test_list,height_frac=0.75,
                                          track_label="Post_translational_modifications",
                                          marker_bar_dict={"color":"blue"})
        # add a second marked position track with the following parameters:
        #track name:  Post_translational_modifications
        #hight of the maker bar = 75%
        #color of the markerbar= blue

        Notes
        ------
        Any panel can have zero, one or more than one marked-position track.
        Thus, in the above examples calling the method ``add_marked_positions_track`` 
        for the second time does NOT override the previous marked-position track 
        it create a new one and added to the figure.
        """
        # adjust the current Geometry of the track.
        number_axes=len(self.fig.axes)
        for i in range(number_axes):
            self.fig.axes[i].change_geometry(number_axes+1,1,i+1)
        # create a new axes
        ax=self.fig.add_subplot(number_axes+1,1,number_axes+1)
        # adjust the axis size
        height=1
        ax.axis(xmin=0,xmax=self.protein_length,ymin=0,ymax=height)
        # draw a base line
        ax.hlines(y=0,xmin=0,xmax=self.protein_length,**base_line_dict)
        # plotting the figure        
        for pos in positions:
            ax.vlines(x=pos,ymin=0,ymax=height_frac,
                      **marker_bar_dict)
            ax.scatter(x=pos,y=height_frac,**marker_dict)
        # adjust the lay out of the figure
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)      
        # adjust the edges
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        # add the label:
        ax.set_ylabel(track_label,track_label_dict)
    
    def add_coverage_track(self, coverage_matrix: np.ndarray,
                          coverage_as_base: bool =False,
                          coverage_dict: Dict[str,Union[int,str]]={"color":"blue","width":1.2},
                          xlabel:str="positions",
                          xlabel_dict:Dict[str,Union[int,str]]={"fontsize":6,"color":"black"},
                          ylabel:str ="coverage",
                          ylabel_dict:Dict[str,Union[int,str]]={"fontsize":6,"color":"black"},
                          number_ticks:int =10,
                          xticks_font_size: int =4,
                          yticks_font_size:int =4):
        """
        Description
        -----------
        Adds a coverage plot to the panel. The coverage plot shows the 
        relationship between a peptide and its experimentally detected eluted 
        peptide/peptides.

        Parameters
        ----------
        coverage_matrix : np.ndarray
            A protein length by one array which summarize information about the protein and 
            the eluted peptides.
            
        coverage_as_base : bool, optional
            Whether or not to plot the coverage as a base track for the figure.
            The default is False which means that the track appended to a figure 
            that have a default base track which can be constructed using the 
            method ``add_base_track``. However, if coverage_as_base is set to 
            True, the function will draw the base track using the coverage matrix and
            calling the function add_base_track should be avoided. 
            
        coverage_dict : Dict[str,Union[int,str]], optional
            The parameters that control the printing of the coverage matrix, 
            for example, the color. These parameters are fed to the function
            ``axes.bar``. The default is {"color":"blue","width":1.2}.
            
        xlabel : str, optional
            The label of the x-axis of the coverage track. The default is "positions".
            
        xlabel_dict : Dict[str,Union[int,str]], optional
            The parameters that control the x-label printing, for example,
            the color and/ the font size. these parameters are fed to the function
            ``axes.set_xlabel``. The default is {"fontsize":6,"color":"black"}.
            
        ylabel : str, optional
            The label of the y-axis of the coverage track. The default is "coverage".
            
        ylabel_dict : Dict[str,Union[int,str]], optional
            The parameters that control the x-label printing, for example,
            the color and/ the font size. these parameters are fed to the function
            ``axes.set_ylabel``. The default is {"fontsize":6,"color":"black"}.
            
        number_ticks : int, optional
            The number of ticks on the x-axis. The default is 10.
            
        xticks_font_size : float, optional
            The font size of the x-axis ticks. The default is 4.
            
        yticks_font_size : float, optional
            The font size of the y-axis ticks. The default is 4.

        """
        # adjust the figure axes:
        if coverage_as_base:
            number_axes=len(self.fig.axes)
            for i in range(number_axes):
                self.fig.axes[i].change_geometry(number_axes+1,1,i+1)
            # create a new axes
            ax=self.fig.add_subplot(number_axes+1,1,number_axes+1)
        else:
            ax=self.fig.add_subplot(111) # coverage is the base track
            
        # Adjust the limits
        heightest_coverage=max(coverage_matrix.reshape(-1))
        ax.axis(xmin=0,xmax=self.protein_length,ymin=0,
                ymax=heightest_coverage+0.1*heightest_coverage)
        
        # Plotting the matrix
        ax.bar(np.arange(coverage_matrix.shape[0]),
                          coverage_matrix.reshape(-1),**coverage_dict)
        # remove the spines of the figure
        ax.spines["left"].set_visible(False)   
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        # adjust the x and y labels
        if not coverage_as_base:
            ax.get_xaxis().set_ticks([])
        else:
            ax.xaxis.tick_top()
        ax.set_ylabel(ylabel,**ylabel_dict)
        ax.set_xlabel(xlabel,xlabel_dict)
        step_size=int(np.floor(self.protein_length/number_ticks))
        ax.set_xticks(ticks=np.arange(0,self.protein_length,step_size),minor=False)
        # adjust the tick size
        ax.tick_params(axis='x', which='major', labelsize=xticks_font_size)
        ax.tick_params(axis='y', which='major', labelsize=yticks_font_size)
    
    def get_figure(self)->plt.Figure:
        """
        Returns
        -------
        matplotlib.figure.Figure
            The figure with all the tracks that have been added to it.
        """
        return self.fig
    
    def save_fig(self,name: str,
                output_path:str=".", 
                format_:str="png",
                figure_dpi:str="same",
                figure_saving_dict:Dict[str,Union[int,str]]={"facecolor":"white"})->None:
        """
        Description
        -----------
        Write the constructed figure to the disk.

        Parameters
        ----------
        name : str
            The name of the figure to save the file.
        
        output_path : str , optional
            The path to write the output, by default the function write to the
            current working directory.
        
        format_ : str, optional
            The output format, this parameter will be fed to the method
            ``plt.savefig``. The default is "png".
        
        figure_dpi: int, optional 
            The dpi of the saved figure. The deafult is same which means the figure 
            will be saved using the same dpi used for creating the figure. 
        
        figure_saving_dict: Dict[str,Union[int,str]],optional
            The parameters that should be fed to the function ``plt.savefig``.
            The default is figure_saving_dict={"facecolor":"white"}
        
        Returns
        -------
        None.

        """
        # Check the user input
        if figure_dpi=="same":
            figure_dpi=self.figure_dpi
        else:
            assert isinstance(figure_dpi,int), """
            figure_dpi must be an int object, or the string same for the figure current dpi"""
        assert isinstance(figure_saving_dict,dict), """ figure_saving_dict must be a dict object"""
        # figure name
        filename=os.path.join(output_path,name)
        # save the figure
        self.fig.savefig(fname=filename,dpi=figure_dpi,format=format_,**figure_saving_dict)

