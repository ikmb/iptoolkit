#!/usr/bin/env Python 
"""An interface to query and retrive data from AFND
"""
## LOAD THE MODULES 
#------------------
from __future__ import  annotations
import urllib3
import mhcnames
import pandas as pd 
from bs4 import BeautifulSoup
##DEFINE THE CLASS
#-----------------
class AFND:
    def __init__(self)->AFND:
        """class constructor  

        :return: an AFND instance that can be used to query data from the website
        :rtype: AFND
        """
        self._base='http://www.allelefrequencies.net/hla6002a.asp?all_name='
        self._http=urllib3.PoolManager()
        return

    def get_frequency_of(self,allele_name:str)->pd.DataFrame:
        """The worker function of the class, Queries AFND for the frequency of the provided allele in different populations

        :param allele_name: the allele name in standard notation, for example, A*01:01 or DRB1*15:01
        :type allele_name: string
        :return: a tuples of two columns, the first contain the country name and the second contains the frequency
        :rtype: pd.DataFrame
        """
        try:
            allele_name=mhcnames.parse_allele_name(allele_name)
        except mhcnames.AlleleParseError as exp: 
            raise ValueError(f'While parsing the provided name: {allele_name} the following error was encountered: {str(exp)}')
        q_name = allele_name.gene+'*'+allele_name.allele_family+':'+allele_name.allele_code
        try: 
            print(f'sending the query to the server ....')
            response=self._http.request('GET',self._base+q_name) 
            print(f'response received')
        except Exception as exp: 
            raise RuntimeError(f'While querying the server the following error was encountered: {str(exp)}')
        # check that the response is valid 
        if response.status==404:
            raise RuntimeError(f'The provided allele list is not defined in the database, server returned 404 ')        
        # get the frequency table from the HTML page 
        html_soup=BeautifulSoup(response.data.decode('utf-8',errors='ignore'),'html.parser')
        html_table=html_soup.body.find_all("table")[3].find_all("tr")[1:]
        ## extract the data from all elements of the table 
        rows=[]
        for row in html_table:
            columns=[]
            for elem in row:
                try:
                    columns.append(elem.get_text())
                except AttributeError:
                    continue
            rows.append(columns)
        ## create the table 
        temp_table=pd.DataFrame(data=rows)
        freq_table=pd.concat([temp_table.iloc[:,1],temp_table.iloc[:,3]],axis=1)
        freq_table.columns=['Country','Frequency']
        return freq_table
    
    
   
    