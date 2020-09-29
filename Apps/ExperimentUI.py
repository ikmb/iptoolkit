#!/usr/bin/env python 
"""
@author: Hesham ElAbd
@contact: h.elabd@ikmb.uni-kiel.de
@brief: A dashboard that utilize the IPTK Library to build to provide an easy-to-use interface to visualize and filter the results of an IP Experiment. 
"""
# load the module 
import dash 
import dash_core_components as dcc 
import dash_html_components as html 
from dash.dependencies import Input, Output, State 
import dash_bootstrap_components as dbc
import pandas as pd 
import os
import base64
import dash_uploader as du
import IPTK.IO.InFunctions as inFunc 
from typing import Dict
import numpy as np
from IPTK.Classes.Database import SeqDB, GeneExpressionDB, CellularLocationDB, OrganismDB
from IPTK.Classes.Experiment import Experiment 
from IPTK.Classes.Proband import Proband 
from IPTK.Classes.Tissue import Tissue 
from IPTK.Classes.HLASet import HLASet
from IPTK.Visualization.vizTools import plotly_num_peptides_per_organism
from IPTK.Visualization.vizTools import plotly_paired_representation 
from IPTK.Visualization.vizTools import plotly_num_peptides_per_parent
from IPTK.Visualization.vizTools import plotly_parent_protein_expression_in_tissue 
from IPTK.Visualization.vizTools import plotly_gene_expression_vs_num_peptides  
from IPTK.Visualization.vizTools import plotly_num_protein_per_location
from IPTK.Visualization.vizTools import plotly_peptide_length_dist 
from IPTK.Visualization.vizTools import plotly_protein_coverage
# define the app and its external style sheets 
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets,suppress_callback_exceptions=True)
# Define The path to upload the data 
UPLOAD_DIRECTORY='./UPLOADED_FILES'
if not os.path.exists(UPLOAD_DIRECTORY):
    os.makedirs(UPLOAD_DIRECTORY)
du.configure_upload(app, UPLOAD_DIRECTORY)
def save_file(name, content):
    """
    @brief: Decode and store a file uploaded with Plotly Dash.
    @obtained from: https://docs.faculty.ai/user-guide/apps/examples/dash_file_upload_download.html
    """
    data = content.encode("utf8").split(b";base64,")[1]
    with open(os.path.join(UPLOAD_DIRECTORY, name), "wb") as fp:
        fp.write(base64.decodebytes(data))
## define some global data 
PEPTIDE_TABLE_PATH:str =None
PEPTIDE_TABLE_FORMAT: str=None
FASTA_DATABASE_PATH: str=None
GENE_EXPRESSION_TABLE: str=None
PROTEIN_LOC_TABLE:str =None
NEW_LINE:str ='\n'
experiment: Experiment = None
MAPPED_PROTEINS: Dict[str, np.ndarray]=None 
# define the layout of the app 
app.layout=html.Div(
    [
        dbc.Row(html.H1('Experimental information')), 
        
        dbc.Row([
                html.Div([
                    html.H3('Identification File '), 
                    html.Div(
                        [
                        du.Upload(
                        id='upload_ident_table', 
                        filetypes=['csv', 'pepXML','idXML','mzTab']
                        
                        )],
                        style={  
                            'textAlign': 'center',
                            'width': '600px',
                            'padding': '10px',
                            'display': 'inline-block'
                            }
                    )
                   ],
                   style={'width': '50%', 'display': 'inline-block'}    
                ),
                html.Div(id='identification_table_hidden_output',style={'display':'none'} )
                ], ),
            
            dbc.Row(html.Div([
                    html.H5("Format: "),
                    dcc.Dropdown(id='ident_table_format',
                    options=[
                        {
                            'label':'pepXML',
                            'value':'pepXML',
                        },
                        {
                            'label':'idXML',
                            'value':'idXML'
                        },
                        {
                            'label':'mzTab',
                            'value':'mzTab'
                        },
                        {
                            'label':'csv',
                            'value':'csv'
                        }
                    ],
                    value='csv' ,
                    style={'width': '40%', 'display': 'inline-block'}
                    )
                ])),
            html.Br(), 
            
            dbc.Row(
                [
                    html.H3('Fasta Database '), 
                    html.Div(
                    [
                    
                    du.Upload(
                    id='upload_seq_table', 
                    filetypes=['fasta'])
                    ],
                    style={  
                            'textAlign': 'center',
                            'width': '600px',
                            'padding': '10px',
                            'display': 'inline-block'
                        }
                    ),
                    html.Div(id='fasta_database_hidden_output',style={'display':'none'} )
                ]
                ),    
        html.Br(),
        html.Div(
            [
                html.H3("HLA alleles"),
                "Enter the alleles, seperated by a semi colons",
                dcc.Input(id='hla_alleles',value='HLA-DRB1*15:01;HLA-DRB1*13:01',
                style={
                            'width': '50%',
                            'height': '30px',
                            'lineHeight': '30px',
                            'borderWidth': '2px',
                            'textAlign': 'left',
                            'margin': '10px'
                            }
                )
            ]
        ),
        html.Br(),
        html.Div(
            [
                html.H3("Tissue Name"),
                "Enter the tissue as defined in Gene expression table: " ,
                dcc.Input(id='tissue_name', value='total PBMC',
                style={
                            'width': '50%',
                            'height': '30px',
                            'lineHeight': '30px',
                            'borderWidth': '2px',
                            'textAlign': 'left',
                            'margin': '10px'
                        }
                )
            ]
        ),
        html.Br(),
        html.Div([
            html.H3('Gene Expression Table: '), 
            dcc.RadioItems(
                id='use_human_protein_atlas', 
                options=[
                    {
                        'label': 'Use Human Protein Atlas reference',
                        'value':'HPA'
                    }, 
                    {
                        'label':'I will upload my expression data',
                        'value':'USER_DATA'
                    }
                ],
                value='HPA',
            )
        ]),
        html.Br(),
        html.Div(id='conditional_expression_table_upload'),
        html.Div([
            html.H3('Protein Localization Table:'), 
            dcc.RadioItems(
                id='use_human_protein_atlas_location', 
                options=[
                    {
                        'label': 'Use Human Protein Atlas reference',
                        'value':'HPA'
                    }, 
                    {
                        'label':'I will upload my localization data',
                        'value':'USER_DATA'
                    }
                ],
                value='HPA',
            )
        ]),
        html.Div(id='conditional_protein_localization_table'),
        html.Br(),
        html.Div(
            [html.Button(id='submit', n_clicks=0, children='Create Experiment')]
        ),
        html.Br(),
        html.Div(id='Experiment_output'),
        html.Br(), 
        html.Div(
            [
                html.H2("Visualization Panel"),
                dcc.Dropdown(
                    id='visualization_functions',
                    options=[
                        {
                            'label':'Number of peptides per organism', 
                            'value':'num_per_per_org'
                        },
                        {
                            'label':'Protein Length Distribution', 
                            'value': 'pep_len_dist'
                        },
                        {
                            'label': 'Number of peptides per parent',
                            'value':'num_pep_per_par'   
                        },
                        {
                            'label':'Parent proteins expression',
                            'value':'par_prot_expression'
                        },
                        {
                            'label':'Number of peptides and gene expression',
                            'value':'num_pep_gene_expr'
                        },
                        {
                            'label': 'Number of proteins per sub-cellular compartment',
                            'value':'num_prot_sub_cell_loc'
                        },
                    ], 
                    placeholder="Select a function",
                    style={
                            'width': '100%',
                            'height': '30px',
                            'lineHeight': '30px',
                            'borderWidth': '2px',
                            'textAlign': 'left',
                            'margin': '10px'
                        }
                ),
                html.Div(id='visualization_figure')
            ]
        ),
        html.Br(), 
        html.Div(
            [
                html.Button(id='remove_org', n_clicks=0, children='Remove Organisms'),
                html.Div(id='organisms_filter')
            ]
        ),
        html.Br(),
        html.Div(
            [
                html.Button(id='protein_coverage', n_clicks=0, children='Look at Coverage'),
                html.Div(id='protein_coverage_fig')
            ]
        )
    ]
)
#*******************************************************************************#
# call back for the uploaded identification table 
@du.callback(
    output=Output('identification_table_hidden_output','children'),
    id='upload_ident_table'
)
def upload_identification_table(file_name):
    if file_name is not None: 
        global PEPTIDE_TABLE_PATH
        PEPTIDE_TABLE_PATH=file_name[0]
        return         
        
#*******************************************************************************#
# call back for the uploaded FASTA DataBase 
@du.callback(Output('fasta_database_hidden_output','children'), 
            id='upload_seq_table')
def upload_database_table(file_name):
    if file_name is not None: 
        global FASTA_DATABASE_PATH
        FASTA_DATABASE_PATH=file_name[0]
        return 
#*******************************************************************************#
# call back to upload expression table
@app.callback(
    Output('conditional_expression_table_upload','children'), 
    [Input('use_human_protein_atlas','value')]
    )
def control_output_state(ref_table):
    if ref_table=='USER_DATA': 
        return html.Div([
                    dcc.Upload(
                        id='gene_expression_table', 
                        children=html.Div(
                            [
                            'Drag and Drop or ', 
                            html.A('Select Files')
                            ]
                        ),
                        style={
                            'width': '100%',
                            'height': '60px',
                            'lineHeight': '60px',
                            'borderWidth': '1px',
                            'borderStyle': 'dashed',
                            'borderRadius': '5px',
                            'textAlign': 'center',
                            'margin': '8px'
                            }
                        ), 
                        html.Div(id="uploaded_gene_expression")
                    ], style={'width': '50%', 'display': 'inline-block'})
    else: 
        GENE_EXPRESSION_TABLE="USE_HPA"
        return 
#*******************************************************************************#
@app.callback(
    Output('uploaded_gene_expression', 'children'),
    [Input('gene_expression_table', 'filename')],
    [State('gene_expression_table','contents')],
)
def upload_expression_table(file_name, content):
    if file_name is not None: 
        save_file(file_name,content)
        global GENE_EXPRESSION_TABLE
        GENE_EXPRESSION_TABLE=os.path.join(UPLOAD_DIRECTORY,file_name)
        return f'File: {file_name} has been uploaded'
#*******************************************************************************#
@app.callback(
    Output('conditional_protein_localization_table','children'), 
    [Input('use_human_protein_atlas_location','value')]
    )
def control_localization_output_state(ref_table):
    if ref_table=='USER_DATA': 
       return html.Div([
                dcc.Upload(
                    id='protein_localization_table', 
                    children=html.Div(
                        [
                        'Drag and Drop or ', 
                        html.A('Select Files')
                        ]
                    ),
                    style={
                        'width': '100%',
                        'height': '60px',
                        'lineHeight': '60px',
                        'borderWidth': '1px',
                        'borderStyle': 'dashed',
                        'borderRadius': '5px',
                        'textAlign': 'center',
                        'margin': '8px'
                        }
                    ),
                    html.Div(id="uploaded_protein_localization_table")
                    ], style={'width': '50%', 'display': 'inline-block'}
            )
    else: 
        PROTEIN_LOC_TABLE="USE_HPA"
        return 
#*******************************************************************************#
@app.callback(
    Output('uploaded_protein_localization_table', 'children'),
    [Input('protein_localization_table', 'filename')],
    [State('protein_localization_table','contents')],
)
def upload_protein_localization_table(file_name, content):
    if file_name is not None: 
        save_file(file_name,content)
        global PROTEIN_LOC_TABLE
        PROTEIN_LOC_TABLE=os.path.join(UPLOAD_DIRECTORY,file_name)
        return f'File: {file_name} has been uploaded'
#*******************************************************************************#
@app.callback(Output('Experiment_output','children'), 
            [Input('ident_table_format','value'),
            Input('tissue_name','value'), 
            Input('hla_alleles','value'), 
            (Input('submit', 'n_clicks'))]
            )
def create_experiment(table_format, tissue_name, hla_alleles, n_clicks):
    if n_clicks > 0: 
        if PEPTIDE_TABLE_PATH is None:
            return "ERROR: The peptide identification file has not been uploaded"
        if FASTA_DATABASE_PATH is None: 
            return "ERROR: The sequence database identification has not been uploaded"
        # try to load the peptide table 
        try: 
            if table_format=='pepXML': 
                table_pep: pd.DataFrame = inFunc.parse_xml_based_format_to_identification_table(
                    path2XML_file=PEPTIDE_TABLE_PATH, 
                    path2fastaDB=FASTA_DATABASE_PATH,
                    is_idXML=False
                )
            elif table_format=='idXML':
                table_pep: pd.DataFrame = inFunc.parse_xml_based_format_to_identification_table(
                    path2XML_file=PEPTIDE_TABLE_PATH, 
                    path2fastaDB=FASTA_DATABASE_PATH,
                    is_idXML=True
                )
            elif table_format=='mzTab':
                table_pep: pd.DataFrame = inFunc.parse_mzTab_to_identification_table(
                    path2mz_tab=PEPTIDE_TABLE_PATH, 
                    path2fastaDB=FASTA_DATABASE_PATH,) 
            else: 
                table_pep: pd.DataFrame = inFunc.parse_text_table(
                    path2file=PEPTIDE_TABLE_PATH, 
                    path2fastaDB=FASTA_DATABASE_PATH,
                    sep=',') 
        except Exception as exp: 
            return f'ERROR:: While parsing the identification table, the following error was encountered: {exp} '
        # create a proband 
        proband: Proband = Proband(name='UI_PROBAND')
        # create an sequence database 
        try: 
            seqs: SeqDB = SeqDB(path2fasta=FASTA_DATABASE_PATH)
        except Exception as exp:
            return f'ERROR:: While creating the sequence database: the following error was encountered; {exp}'
        # create the OrgDB
        try: 
            org_db: OrganismDB =OrganismDB(FASTA_DATABASE_PATH)
        except Exception as exp: 
            return f'ERROR:: While creating the sequence database: the following error was encountered; {exp}'
        # create the expression profile 
        if GENE_EXPRESSION_TABLE is None:
            try: 
                expresson_profile: GeneExpressionDB= GeneExpressionDB(
                path2data='https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip',sep='\t')
            except Exception as exp: 
                return f'While Downloading the online table the following error was encountered: {exp}'
        else:
            try: 
                expresson_profile: GeneExpressionDB= GeneExpressionDB(
                path2data=GENE_EXPRESSION_TABLE,sep=',')
            except Exception as exp: 
                return f'While parsing the expression table, the following error was encountered: {exp}'
        # create the location table 
        if PROTEIN_LOC_TABLE is None: 
            try: 
                protein_locations: CellularLocationDB=CellularLocationDB(
                path2data='https://www.proteinatlas.org/download/subcellular_location.tsv.zip', 
                sep='\t')
            except Exception as exp: 
                return f'While Downloading the protein sub-cellular location table, the following error was encountered: {exp}'
        else: 
            try: 
                protein_locations: CellularLocationDB=CellularLocationDB(
                path2data=PROTEIN_LOC_TABLE,sep=',')
            except Exception as exp: 
                return f'While parsing the location table, the following error was encountered: {exp}'
        # create the tissue instance 
        tissue: Tissue = Tissue(name=tissue_name,
                        main_exp_value=expresson_profile, 
                        main_location=protein_locations)
        # create the hla_set
        hlas: HLASet = HLASet(hlas=hla_alleles.split(';'))
        # create the experiment object 
        global experiment 
        try: 
            experiment = Experiment(proband=proband,hla_set=hlas,tissue=tissue,database=seqs, 
                 ident_table=table_pep)
        except Exception as exp: 
            return f'while creating an experimental object,the following error was encounter {exp}'
        # annoatate the experiment 
        experiment.annotate_proteins(org_db)
        # return the experiment output 
        return str(experiment)
#*******************************************************************************#
@app.callback(
    Output('visualization_figure','children'),
    [Input('visualization_functions','value')]
)
def visualization_function(vis_func): 
    if experiment is None: 
        return [html.A('The Experiment object has not been created yet!') ]      
    else: 
        if vis_func == "num_per_per_org": 
            fig=plotly_num_peptides_per_organism(experiment.get_peptides_per_organism())
            return [dcc.Graph(figure=fig)]
        
        elif vis_func == 'pep_len_dist': 
            fig= plotly_peptide_length_dist(experiment.get_peptides_length())
            return [dcc.Graph(figure=fig)]
        
        elif vis_func == 'num_pep_per_par': 
            fig = plotly_num_peptides_per_parent(experiment.get_peptides_per_protein(), 25, 
            title = 'Number of peptides per protein, Top 25 protein')
            return [dcc.Graph(figure=fig)]
        
        elif vis_func == 'par_prot_expression':
            parent_protein_exp: pd.DataFrame = experiment.get_expression_of_parent_proteins()
            reference_gene_expression: pd.DataFrame = experiment.get_experiment_reference_tissue_expression()
            fig=plotly_parent_protein_expression_in_tissue(parent_protein_exp,reference_gene_expression, 
                                             tissue_name=experiment.get_tissue_name()) 
            return [dcc.Graph(figure=fig)]
        
        elif vis_func == 'num_pep_gene_expr':
            num_peptides_parent: pd.DataFrame = experiment.get_num_peptide_expression_table() 
            fig=plotly_gene_expression_vs_num_peptides(num_peptides_parent,
                    tissue_name=experiment.get_tissue_name())
            return [dcc.Graph(figure=fig)]

        elif vis_func == 'num_prot_sub_cell_loc':
            protein_counts: pd.DataFrame = experiment.get_number_of_proteins_per_compartment() 
            fig=plotly_num_protein_per_location(protein_counts,
                    drop_unknown=False) 
            return [dcc.Graph(figure=fig)]
        else:
            return [html.A('Not done yet, working on it ....')]
#*******************************************************************************#
@app.callback(
    Output('organisms_filter','children'),
    [Input('remove_org','n_clicks')]
)
def return_org_filter_interface(n_clicks):
    if n_clicks>0:
        if experiment is None: 
            return html.A('The Experiment object has not been created yet!')
        else:
            return html.Div(
                [
                    dcc.Dropdown(
                        id='organism_filter', 
                        options=[{'label':org,'value':org} for org in experiment.get_orgs()],
                        placeholder='select on organism to remove'
                    ), 
                    html.Div(id='exp_after_org_rm'), 
                    "WARNING:: REMOVING AN ORGANISM FROM THE EXPERIMENT OBJECT IS IRREVERSIBLE"
                ]
            )
@app.callback(
    Output('exp_after_org_rm','children'), 
    [Input('organism_filter','value')]
)
def filter_by_org_and_show_exp(org2rm):
    experiment.drop_peptide_belong_to_org(org2rm)
    # return the figure shown the number of peptide per org after the update
    fig=plotly_num_peptides_per_organism(experiment.get_peptides_per_organism())
    return [dcc.Graph(figure=fig)]
#*******************************************************************************#
@app.callback(
    Output('protein_coverage_fig','children'),
    [Input('protein_coverage','n_clicks')]
)
def return_protein_coverage_interface(n_clicks):
    if n_clicks>0:
        if experiment is None: 
            return html.A('The Experiment object has not been created yet!')
        else:
            return html.Div(
                [
                    dcc.Dropdown(
                        id='protein_id', 
                        options=[{'label':org,'value':org} for org in experiment.get_proteins()],
                        placeholder='select a protein to visualize its coverage'
                    ), 
                    html.Div(id='pep_per_proteins'), 
                ]
            )
@app.callback(
    Output('pep_per_proteins','children'), 
    [Input('protein_id','value')]
)
def plot_protein_coverage(prot_id):
    global MAPPED_PROTEINS
    if MAPPED_PROTEINS is None:
        MAPPED_PROTEINS=experiment.get_mapped_proteins()
    # return the figure shown the number of peptide per org after the update
    fig=plotly_protein_coverage(MAPPED_PROTEINS[prot_id], prot_id)
    return [dcc.Graph(figure=fig)]
#*******************************************************************************#
# run the app 
if __name__=='__main__':
    app.run_server(debug=True)
