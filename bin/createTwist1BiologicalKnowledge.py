"""
Reads ./database/supp_jcb.201306088_JCB_201306088_TableS1.xlsx
Changes column name
Writes ./networks/Twist1/biologicalknowledge.csv
"""

#!/usr/bin/env python

import sys
import os.path
import logging
import warnings
import pandas as pd

warnings.filterwarnings('ignore', category=UserWarning, module='openpyxl')

logging.basicConfig(format='%(name)s %(levelname)s: %(message)s')
logger = logging.getLogger('NetResponse')
logger.setLevel(logging.INFO)

def readData():
    """
    Reads ./database/supp_jcb.201306088_JCB_201306088_TableS1.xlsx
    This file is Table 1 from Shamir et al, 2014.
    Returns pandas database for sheet All Sequenced Genes.
    The column name log2(Fold Change) has been changed to Quantity.
    """

    #check filepath exists
    filepath = './database/supp_jcb.201306088_JCB_201306088_TableS1.xlsx'

    if not os.path.isfile(filepath):
        logger.error('%s does not exist.',filepath)
        sys.exit()

    #check sheet exists
    sheetname = 'All Sequenced Genes'
    xl = pd.ExcelFile(filepath)
    if sheetname not in xl.sheet_names:
        logger.error('%s does not have sheet %s.',filepath,sheetname)
        sys.exit()

    #read data
    ret = pd.read_excel(filepath,sheet_name=sheetname)

    #check column exists
    cols = ['Gene','log2(Fold Change)']
    for col in cols:
        if col not in ret.columns:
            logger.error('%s sheet %s does not have column %s.',filepath,sheetname,col)
            sys.exit()

    #change column name
    ret = ret.rename(columns={'log2(Fold Change)': 'Quantity'})
    return(ret)

def writeData(data):
    """
    Writes a pandas dataframe to ./networks/Twist1/biologicalknowledge.csv
    """
    filepath = './networks/Twist1/biologicalknowledge.csv'
    data.to_csv(filepath,sep=',',index=False)
    
    return None

def main():
    #read ./database/supp_jcb.201306088_JCB_201306088_TableS1.xlsx
    data = readData()

    #write data to ./networks/Twist1/biologicalknowledge.csv
    writeData(data)

if __name__=='__main__':
    main()
