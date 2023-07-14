#!/usr/bin/env python

import sys
import os.path
import argparse
import logging
import pandas as pd
import numpy as np
import warnings
from math import log, e
from scipy import linalg,integrate,stats
from collections import defaultdict

warnings.filterwarnings('ignore', category=UserWarning, module='openpyxl')

logging.basicConfig(format='%(name)s %(levelname)s: %(message)s')
logger = logging.getLogger('NetResponse')
logger.setLevel(logging.INFO)

class CMDParser():
    #Command line parser

    def __init__(self):
        parser = argparse.ArgumentParser()
        subparser = parser.add_subparsers()

        #NetResponse analysis
        analysis_parser = subparser.add_parser('analysis', help='Run NetResponse.')
        analysis_parser.add_argument('network_dir',help='Network directory path.')
        analysis_parser.add_argument('driver',help='Network driver gene symbol.')
        analysis_parser.add_argument('-s',choices=['human','mouse'],default='human',help='species. Default human')
        analysis_parser.set_defaults(func = runAnalysis)
    
        #Compare NetResponse rankings to biological knowledge and dissemination assay results from Georgess et al. (2020)
        comparison_parser = subparser.add_parser('comparison', help='Compare NetResponse rankings to biological knowledge and dissemination assay results from Georgess et al. (2020).')
        comparison_parser.add_argument('network_dir',help='Network directory path.')
        comparison_parser.add_argument('driver',help='Network driver gene symbol.')
        comparison_parser.add_argument('-s',choices=['human','mouse'],default='human',help='species. Default human')
        comparison_parser.set_defaults(func = runComparison)

        args = parser.parse_args()
        self.__dict__.update(args.__dict__)
        return

def runAnalysis(args):
    """
    Runs NetResponse analysis
    """
    #logger.info('*******************************')
    print('*************************************************')
    logger.info('Running NetResponse analysis...')
    #logger.info('*******************************')
    print('*************************************************')
    getNetworkInteractions(args)
    driver,intermediates,responses = getNetworkGenes(args,verbose=True)
    getPrioritization(args,driver,intermediates,responses)
    return None

def getNetworkInteractions(args):
    """
    Compiles interactions into one file and stores it in database dir, if needed.
    Removes unwanted PPI between a protein and itself.
    """
    logger.info('Compiling %s network interactions from data files...',args.s)
    
    #get filepath
    if args.s == 'human':
        interaction_filepath = './database/human_interactions.txt'
    elif args.s == 'mouse':
        interaction_filepath = './database/mouse_interactions.txt'
    else:
        logger.error('%s is not a supported species.' % args.s)
        sys.exit()

    #check if interaction file exists
    if os.path.isfile(interaction_filepath):
        logger.info('Found compiled %s interactions file %s',args.s,interaction_filepath)
    else:
        #read PPI and TF
        if args.s == 'human':
            ppi = readIIDHuman()
            tf = readTF('human')
        elif args.s == 'mouse':
            ppi = readIIDMouse()
            tf = readTF('mouse')
        else:
            logger.error('%s is not a supported species.', args.s)
            sys.exit()
        #combine data and write
        ppi_frame = pd.DataFrame(ppi,columns=['source','target'])
        ppi_frame['type'] = 'PPI'
        tf_frame = pd.DataFrame(tf,columns=['source','target'])
        tf_frame['type'] = 'TF'
        frames = [ppi_frame,tf_frame]
        wholeinteract = pd.concat(frames)
        wholeinteract.to_csv(interaction_filepath,sep='\t',index=False,header=False)
        logger.info('Wrote compiled %s interactions file to %s',args.s,interaction_filepath)

    edges = readEdges(args,verbose=True)
       
    return None

def getNetworkGenes(args,verbose=False):
    """
    Reads user-defined driver and respones. Checks they are in the network.
    Identifies intermediates from interaction file.
    Returns driver, list of intermediates, and list of network respones.
    """
    if verbose:
        logger.info('Compiling network genes...')

    #Read network interactions
    if args.s == 'human':
        gene_filepath = './database/human_genes.txt'
        interaction_filepath = './database/human_interactions.txt'
    elif args.s == 'mouse':
        gene_filepath = './database/mouse_genes.txt'
        interaction_filepath = './database/mouse_interactions.txt'
    else:
        logger.error('%s is not a supported species.' % args.s)
        sys.exit()

    if not os.path.isfile(gene_filepath):

        #read interaction file
        wholeinteract = readInteractions(interaction_filepath)

        #get whole gene list
        wholegenelist = wholeinteract[0].tolist() + wholeinteract[1].tolist()
        wholegenelist = list(set(wholegenelist))

        #write to file
        wholegenedf = pd.DataFrame(wholegenelist)
        wholegenedf.to_csv(gene_filepath,sep='\t',index=False,header=False)
        logger.info('Wrote compiled %s gene list file to %s', args.s, gene_filepath)
        
    #check if gene list is imported
    try:
        wholegenelist
    except NameError:
        wholegenelist = readGenes(gene_filepath)
        if verbose:
            logger.info('Found compiled %s gene list file %s',args.s,gene_filepath)

    if verbose:
        logger.info('Network directory: %s',args.network_dir)
    
    #read driver and check it is in network
    driver = checkDriver(args,wholegenelist,verbose)
    
    #get responses in network
    responses = getNetworkResponses(args,driver,wholegenelist,verbose)

    #get intermediates in network
    intermediatesFilepath = args.network_dir + '/networkIntermediates.txt' 
    if os.path.isfile(intermediatesFilepath):
        if verbose:
            logger.info('Found network intermediates file %s. Delete this file if you want the intermediates to be compiled again',intermediatesFilepath)
    else:
        #get intermediates
        intermediates = wholegenelist

        #remove driver from intermediates
        intermediates.remove(driver)

        #remove responses from intermediates
        for r in responses:
            intermediates.remove(r)
    
        #output intermediates
        fp = open(intermediatesFilepath, 'w')
        for intermediate in intermediates:
            fp.write(intermediate + '\n')
        fp.close()
        if verbose:
            logger.info('Found %s intermediate genes in the network. Wrote to %s.' % (str(len(intermediates)),intermediatesFilepath))

    #check if intermediates is imported
    try:
        intermediates
    except NameError:
        intermediates = readGenes(intermediatesFilepath)

    #output network size
    if verbose:
        network_size = 1 + len(intermediates) + len(responses)
        logger.info('Network consists of %d genes: the driver %s, %d intermediates, and %d responses.',network_size,driver,len(intermediates),len(responses))

    return(driver,intermediates,responses)

def getPrioritization(args,driver,intermediates,responses):
    """
    Calculates NetResponse weights.
    Ranks network genes.
    Writes results
    """

    #get whole gene list
    wholegenelist = list()
    wholegenelist.append(driver)
    wholegenelist.extend(intermediates)
    wholegenelist.extend(responses)
    
    #get network edges 
    edge_dict = readEdges(args)
    edge_ab, edge_ba = edgeTuplesToDicts(edge_dict)

    #get genename (vertex) to index dict and vice versa
    v2i, i2v = nameToIndex(wholegenelist)

    #get genename (vertex) to indegree dict and outdegree dict
    v2indeg, v2outdeg = nameToDegree(wholegenelist,edge_ba,edge_ab)

    #get genename (vertex) to type (driver, intermediate, response)
    v2type = nameToType(wholegenelist,driver,intermediates,responses)

    #get laplacian matrix
    (adj_mat, lap_mat) = getLaplacian(edge_dict, v2i)

    #get r, nt, dt
    (dt,nt) = getTimeSteps(driver,v2outdeg)

    #get G(t) list
    g_list = getGlist(lap_mat,dt,nt)

    #get weights
    v2weight = getWeights(g_list,v2i,driver,responses,nt,dt)

    #get rankings
    (vertices_ranked,v2rank) = getRankings(v2weight,reverse_order=True)

    #write to file
    writeWeights(args,v2weight,v2rank,v2indeg,v2outdeg,v2type,vertices_ranked)

    return None

def runComparison(args):
    """
    Compares NetResponse rankings to biological knowledge and dissemination assay results from Georgess et al. (2020).
    Writes Tables 1 and 2 in the manuscript.
    """
    print('***************************************************')
    logger.info('Running NetResponse comparison...')
    print('***************************************************')

    logger.info('Comparing NetResponse rankings to biological knowledge and dissemination assay results from Georgess et al. (2020).')
    logger.info('Network directory: %s',args.network_dir)

    #get network genes
    driver, intermediates, networkResponses = getNetworkGenes(args,verbose=False)
    responses = readResponses(args,driver)


    #map human gene to Dan results
    hgene2dan = getDanResults()

    #map human gene to repurposing hub drug
    hgene2repurp = getRepurpDrugs()

    #get repurp drugs and dan results for driver, intermediates, shortest-path intermediates, and responses
    d2dan,d2repurp,i2dan,i2repurp,r2dan,r2repurp,spi2dan,spi2repurp = getGenes2DrugsAndResults(args,driver,intermediates,responses,hgene2repurp,hgene2dan)

    #write protein count table
    writeTable1(args,driver,intermediates,responses,d2repurp,i2repurp,r2repurp,spi2repurp)

    #genes to NetResponse weights
    v2weight = readWeights(args)

    #genes to biological knowledge 
    v2knowledge = readBiologicalKnowledge(args)

    #calculate spearman correlation for intermediates and dan results
    spearman_intermediates,spearman_responses = getSpearman(i2dan,r2dan,v2weight,v2knowledge)

    #write protein test table
    writeTable2(args,driver,intermediates,responses,d2repurp,d2dan,i2repurp,i2dan,r2repurp,r2dan,spearman_intermediates,spearman_responses)
    
    return None

def writeTable1(args,driver,intermediates,responses,d2repurp,i2repurp,r2repurp,spi2repurp):
    """
    Writes the manuscript Table 1
    """

    #get edges
    edge_dict = readEdges(args)
    edge_ab,edge_ba = edgeTuplesToDicts(edge_dict)

    #get shortestpath intermediates
    sp_intermediates = getShortestPathIntermediates(args,intermediates,driver,responses,edge_ab,edge_ba)

    #Driver + Response row
    dr_tot = 1 + len(responses)
    dr_tar = len(d2repurp) + len(r2repurp)
    
    #All intermediates row
    i_tot = len(intermediates)
    i_tar = len(i2repurp)

    #Shortest path intermediates row
    spi_tot = len(sp_intermediates)
    spi_tar = len(spi2repurp)

    #Write to file
    filepath = args.network_dir + '/Table1.tsv'
    f = open(filepath,'w')
    f.write(""+"\ttotal\ttargetable\n")
    f.write("Driver + Response\t"+str(dr_tot)+"\t"+str(dr_tar)+"\n")
    f.write("Shortest-path intermediates\t"+str(spi_tot)+"\t"+str(spi_tar)+"\n")
    f.write("All intermediates\t"+str(i_tot)+"\t"+str(i_tar)+"\n")
    f.close()

    logger.info('Wrote Table 1 to %s',filepath)

    return None

def writeTable2(args,driver,intermediates,responses,d2repurp,d2dan,i2repurp,i2dan,r2repurp,r2dan,spearman_intermediates,spearman_responses):
    """
    Writes the manuscript Table 2
    """

    #Driver + Response column
    dr_tot = 1 + len(responses)
    dr_tar = len(d2repurp) + len(r2repurp)
    dr_dan = len(d2dan) + len(r2dan)
    dr_r = spearman_responses[0]
    dr_p = spearman_responses[1]
    #Tested small molecules
    dr_dan_drug_dict = dict()
    for driver in d2dan:
        for drug in d2dan[driver][0]:
            dr_dan_drug_dict[drug] = True
    for response in r2dan:
        for drug in r2dan[response][0]:
            dr_dan_drug_dict[drug] = True
    dr_dan_drug = len(dr_dan_drug_dict)
    #Repurp drugs
    dr_repurp_drug_dict = dict()
    for driver in d2repurp:
        for drug in d2repurp[driver]:
            dr_repurp_drug_dict[drug] = True
    for response in r2repurp:
        for drug in r2repurp[response]:
            dr_repurp_drug_dict[drug] = True
    dr_repurp_drug = len(dr_repurp_drug_dict)

    #Intermediates column
    i_tot = len(intermediates)
    i_tar = len(i2repurp)
    i_dan = len(i2dan)
    i_r = spearman_intermediates[0]
    i_p = spearman_intermediates[1]
    #Tested small molecules
    i_dan_drug_dict = dict()
    for intermediate in i2dan:
        for drug in i2dan[intermediate][0]:
            i_dan_drug_dict[drug] = True
    i_dan_drug = len(i_dan_drug_dict)
    #Repurp drugs
    i_repurp_drug_dict = dict()
    for intermediate in i2repurp:
        for drug in i2repurp[intermediate]:
            i_repurp_drug_dict[drug] = True
    i_repurp_drug = len(i_repurp_drug_dict)

    #Write to file
    filepath = args.network_dir + '/Table2.tsv'
    f = open(filepath,'w')
    f.write(""+"\tD+R\tIntermediates\n")
    f.write("Number of proteins\t"+str(dr_tot)+"\t"+str(i_tot)+"\n")
    f.write("Targetable proteins\t"+str(dr_tar)+"\t"+str(i_tar)+"\n")
    f.write("Tested proteins\t"+str(dr_dan)+"\t"+str(i_dan)+"\n")
    f.write("Rank correlation (p-value)\t"+"{:.2f}".format(dr_r)+" ("+"{:.1e}".format(dr_p)+")"+"\t"+"{:.2f}".format(i_r)+" ("+"{:.1e}".format(i_p)+")"+"\n")
    f.write("Tested small molecules\t"+str(dr_dan_drug)+"\t"+str(i_dan_drug)+"\n")
    f.write("Small molecules in RepHub\t"+str(dr_repurp_drug)+"\t"+str(i_repurp_drug)+"\n")
    f.close()

    logger.info('Wrote Table 2 to %s',filepath)

    return None

def getShortestPathIntermediates(args,intermediates,driver,responses,edge_ab,edge_ba):
    """
    Returns a list of intermediates with an edge from the driver and an edge to a response.    
    """
    sp_intermediates = list()

    for intermediate in intermediates:
        if intermediate in edge_ba and intermediate in edge_ab:
            for response in responses:
                if driver in edge_ba[intermediate] and response in edge_ab[intermediate]:
                    sp_intermediates.append(intermediate)

    sp_intermediates = list(set(sp_intermediates))

    return(sp_intermediates)


def readNetworkGenes(args):
    """
    Reads file of genes stored in database.
    Returns a list of genes.
    """
    if args.s == 'human':
        filepath = './database/human_genes.txt'
    elif args.s == 'mouse':
        filepath = './database/mouse_genes.txt'
    else:
        logger.error('%s is not a supported species.', args.s)
        sys.exit()

    genes = readGenes(filepath)

    return(genes)

def getSpearman(i2dan,r2dan,v2weight,v2knowledge):
    """
    Returns spearman rank correlation and p-value for intermediates and responses.
    """

    weights = list()
    results = list()

    for i in i2dan:
        weights.append(v2weight[i])
        results.append(100-i2dan[i][1][0])

    spearman_intermediates = stats.spearmanr(weights,results)

    knowledge = list()
    results = list()

    for r in r2dan:
        knowledge.append(v2knowledge[r])
        results.append(100-r2dan[r][1][0])

    spearman_responses = stats.spearmanr(knowledge,results)


    return( (spearman_intermediates,spearman_responses) )

def readBiologicalKnowledge(args):
    """
    Reads file called biologicalknowledge.
    File may have extension .xlsx, .csv, or .tsv.
    File must have columns named Gene and Quanity.
    Gene is a column of gene symbols.
    Quantity is a column of numbers used to rank genes.
    Returns dict of gene to quantity.
    """
    v2knowledge = dict()

    #find biological knowledge filepath
    knowledgeFilepath = args.network_dir + '/biologicalknowledge.xlsx'
    if os.path.isfile(knowledgeFilepath):
        #logger.info('Found biological knowledge file %s. Delete this file if you want to use a file with a .csv or .tsv extension.',knowledgeFilepath)
        pass
    else:
        knowledgeFilepath = args.network_dir + '/biologicalknowledge.csv'
        if os.path.isfile(knowledgeFilepath): 
            #logger.info('Found biological knowledge file %s. Delete this file if you want to use a file with a .tsv extension.',knowledgeFilepath)
            pass
        else:
            knowledgeFilepath = args.network_dir + '/biologicalknowledge.tsv'
            if os.path.isfile(knowledgeFilepath): 
                #logger.info('Found biological knowledge file %s.',knowledgeFilepath)
                pass
            else:
                logger.error('No biological knowledge file found. Must provide a file named biologicalknowledge with .xlsx, .csv, or .tsv extension in directory %s.', args.network_dir)
                sys.exit()

    #check if file is empty
    if not os.stat(knowledgeFilepath).st_size:
        logger.error('File %s is empty.',knowledgeFilepath)
        sys.exit()

    #read data
    if knowledgeFilepath.split(".")[-1] == 'xlsx':
        xl = pd.ExcelFile(knowledgeFilepath)
        res = len(xl.sheet_names)
        if res != 1:
            logger.error('Biological knowledge file %s must contain only 1 sheet.', knowledgeFilepath)
            sys.exit()
        df = pd.read_excel(knowledgeFilepath)
    elif knowledgeFilepath.split(".")[-1] == 'csv':
        df = pd.read_csv(knowledgeFilepath, sep=',', engine='python')
    else:
        df = pd.read_csv(knowledgeFilepath, sep='\t', engine='python')

    #check for columns
    header = list(df.columns)
    header_reqd = ['Gene','Quantity']
    if 'Gene' not in header:
        logger.error('Biologial knowledge file %s must have a column named Gene containing gene symbols.',knowledgeFilepath)
        sys.exit()
    if 'Quantity' not in header:
        logger.error('Biological knowldege file %s must contain a column named Quantity containing the quantity used to rank genes by biological knowledge.',knowledgeFilepath)

    #store data into dict
    for index,row in df.iterrows():
        gene = str(row['Gene'])
        v2knowledge[gene] = float(row['Quantity'])

    logger.info('Read biological knowledge for %d genes from %s',len(v2knowledge),knowledgeFilepath)

    return(v2knowledge)

def readWeights(args):
    """
    Reads weights written to file and stores in dict.
    """
    v2weight = dict()
    weightsFilepath = args.network_dir + '/geneRankings.tsv'
    df = pd.read_csv(weightsFilepath, sep='\t')

    for index,row in df.iterrows():
        gene = row['gene']
        weight = row['weight']
        v2weight[gene] = weight
    
    return(v2weight)

def getGenes2DrugsAndResults(args,driver,intermediates,responses,hgene2repurp,hgene2dan):
    """
    Returns eight dict's.
    1. a mapping of driver to list A and list B.
    List A is a list of drugs tested by dan targeting the protein of the human gene, in order of effectiveness.
    List B is the results of those drugs in list A.
    2. a dict of dict's mapping intermediate to drug repurposing compounds targeting the protein of the human gene.
    3. a mapping of intermediates to list A and list B.
    List A is a list of drugs tested by dan targeting the protein of the human gene, in order of effectiveness.
    List B is the results of those drugs in list A.
    4. a dict of dict's mapping intermediate to drug repurposing compounds targeting the protein of the human gene.
    5. a mapping of responses to list A and list B.
    List A is a list of drugs tested by dan targeting the protein of the human gene, in order of effectiveness.
    List B is the results of those drugs in list A.
    6. a dict of dict's mapping responses to drug repurposing compounds targeting the protein of the human gene.
    7. a mapping of shortest path intermediates to list A and list B.
    List A is a list of drugs tested by dan targeting the protein of the human gene, in order of effectiveness.
    List B is the results of those drugs in list A.
    8. a dict of dict's mapping shortest path intermediates to drug repurposing compounds targeting the protein of the human gene.
    """

    d2dan = defaultdict(list)
    d2repurp = defaultdict(dict)
    i2dan = defaultdict(list)
    i2repurp = defaultdict(dict)
    r2dan = defaultdict(list)
    r2repurp = defaultdict(dict)
    spi2dan = defaultdict(list)
    spi2repurp = defaultdict(dict)

    #get mouse to human conversion, if needed
    if args.s == 'mouse':
        mouse2human = getMouse2Human()
    else:
        mouse2human = dict()
    
    #filepaths
    intermediatesFilepath = args.network_dir + '/networkIntermediates.txt'
    responsesFilepath = args.network_dir + '/responses.txt'

    #check intermediate filepath is there
    if not os.path.isfile(intermediatesFilepath):
        logger.error('Intermediates file is missing in %s. Run NetResponse analysis.',args.network_dir)
        sys.exit()

    #driver
    if driver in mouse2human:
        driver_h = mouse2human[driver]
    else:
        driver_h = driver
    if driver_h in hgene2repurp:
        d2repurp[driver] = hgene2repurp[driver_h]
    dan_results = list()
    if driver_h in hgene2dan:
        dan_results.extend(hgene2dan[driver_h])
    dan_results.sort(key=lambda x: x[1])
    if dan_results:
        d2dan[driver] = [ list(), list() ]
        for pair in dan_results:
            d2dan[driver][0].append(pair[0])
            d2dan[driver][1].append(pair[1])

    #intermediates
    for gene in intermediates:
        if gene in mouse2human:
            gene_h = mouse2human[gene]
        else:
            gene_h = gene
        if gene_h in hgene2repurp:
            i2repurp[gene] = hgene2repurp[gene_h]
        dan_results = list()
        if gene_h in hgene2dan:
            dan_results.extend(hgene2dan[gene_h])
        dan_results.sort(key=lambda x: x[1])
        if dan_results:
            i2dan[gene] = [ list() , list() ]
            for pair in dan_results:
                i2dan[gene][0].append(pair[0])
                i2dan[gene][1].append(pair[1])

    #responses
    for gene in responses:
        if gene in mouse2human:
            gene_h = mouse2human[gene]
        else:
            gene_h = gene
        if gene_h in hgene2repurp:
            r2repurp[gene] = hgene2repurp[gene_h]
        dan_results = list()
        if gene_h in hgene2dan:
            dan_results.extend(hgene2dan[gene_h])
        dan_results.sort(key=lambda x: x[1])
        if dan_results:
            r2dan[gene] = [ list() , list() ]
            for pair in dan_results:
                r2dan[gene][0].append(pair[0])
                r2dan[gene][1].append(pair[1])

    #shortest path intermediates
    edge_dict = readEdges(args)
    edge_ab,edge_ba = edgeTuplesToDicts(edge_dict)
    sp_intermediates = getShortestPathIntermediates(args,intermediates,driver,responses,edge_ab,edge_ba)
    for gene in sp_intermediates:
        if gene in mouse2human:
            gene_h = mouse2human[gene]
        else:
            gene_h = gene
        if gene_h in hgene2repurp:
            spi2repurp[gene] = hgene2repurp[gene_h]
        dan_results = list()
        if gene_h in hgene2dan:
            dan_results.extend(hgene2dan[gene_h])
        dan_results.sort(key=lambda x: x[1])
        if dan_results:
            spi2dan[gene] = [ list() , list() ]
            for pair in dan_results:
                spi2dan[gene][0].append(pair[0])
                spi2dan[gene][1].append(pair[1])


    return(d2dan,d2repurp,i2dan,i2repurp,r2dan,r2repurp,spi2dan,spi2repurp)


def getRepurpDrugs():
    """
    Returns dict mapping of human gene symbol to Broad Repurposing Hub drug.
    """
    hgene2repurp = defaultdict(dict)
    drugs = dict()
    filepath = './database/repurposing_drugs_20200324.txt'
    df = pd.read_csv(filepath, sep='\t', engine='python',comment='!')

    for index,row in df.iterrows():
        drug = str(row['pert_iname'])
        genes = str(row['target']).split('|')
        for gene in genes:
            if gene != 'nan':
                hgene2repurp[gene][drug] = True
                drugs[drug] = True

    logger.info('Read %d drugs and %d targets from the Broad Drug Repurposing Hub.', len(drugs),len(hgene2repurp))

    return(hgene2repurp)

def getMouse2Human():
    """
    Returns a dict mapping of mouse gene symbol to homolog human gene symbol. 
    """

    #Check if mouse to human gene mapping file exists
    filepath = './database/mouse_to_human_genes.tsv'
    if not os.path.isfile(filepath):
        #read data
        rawfilepath = './database/HOM_MouseHumanSequence.rpt'
        #logger.info('Compiling mouse to human gene mapping from %s...',rawfilepath)
        data = pd.read_csv(rawfilepath, sep='\t')

        num2species = dict()
        m2h = dict()
        mouse2human = dict()

        for index,row in data.iterrows():
            num = row['DB Class Key']
            if num not in num2species.keys():
                num2species[num] = dict()
                num2species[num]['mouse'] = dict()
                num2species[num]['human'] = dict()
            species = row['Common Organism Name']
            symbol = row['Symbol']
            if species == 'mouse, laboratory':
                num2species[num]['mouse'][symbol] = True
            if species == 'human':
                num2species[num]['human'][symbol] = True

        for num in num2species.keys():
            for m in num2species[num]['mouse'].keys():
                if m not in m2h.keys():
                    if num2species[num]['human']:
                        m2h[m] = dict()
                for h in num2species[num]['human'].keys():
                    m2h[m][h] = True
            
        for m in m2h.keys():
            if len(m2h[m]) > 1:
                if m.upper() in m2h[m].keys():
                    mouse2human[m]= m.upper()
                else:
                    choice = str(list(m2h[m].keys())[0])
                    mouse2human[m] = str(choice)
            else:
                mouse2human[m] = str(list(m2h[m].keys())[0])

        hgenes = dict()

        f = open(filepath,'w')
        f.write("mouse\thuman\n")
        for m in mouse2human.keys():
            h = mouse2human[m]
            f.write("%s\t%s\n" % (m,h))
            hgenes[h] = True
        f.close()
        logger.info('Mapped %d mouse genes to %d human genes. Wrote to %s.',len(mouse2human),len(hgenes),filepath)
    else:
        #read the existing file
        data = pd.read_csv(filepath, sep='\t')

        mouse2human = dict()

        for index,row in data.iterrows():
            m = row['mouse']
            h = row['human']
            mouse2human[m] = h

    return(mouse2human)

def getDanResults():
    """
    Returns dict mapping of human gene targets to compound names and dissemination assay results from Georgess et al. (2020)
    """

    hgene2dan = defaultdict(list)
    drug2dan = dict()

    #raw file from Georgess et al. (2020)
    rawfilepath = './database/Georgess_Prkd1_2019_Data_Final.xlsx'

    #check raw file exists
    if not os.path.isfile(rawfilepath):
        logger.error('Missing database file %s.', rawfilepath)
        sys.exit()
    
    #read raw file
    dfraw = pd.read_excel(rawfilepath,sheet_name='Fig. 1',usecols='F:H',skiprows=3,nrows=24,header=None)
    for index,row in dfraw.iterrows():
        drug = row[5]
        result = row[7]
        if np.isnan(result):
            result = 100.0
        #correcting a compound name error in Georgess et al. (2020)
        if drug == 'CID-755636':
            drug = 'CID-755673'
        drug2dan[drug] = result
   
    #mapping drugs to targets file
    mappingfilepath = './database/DisseminationAssayTargets.tsv'

    #check mapping file exists
    if not os.path.isfile(mappingfilepath):
        logger.error('Missing database file %s.', mappingfilepath)
        sys.exit() 

    #read mapping file
    dfmap = pd.read_csv(mappingfilepath, sep='\t', engine='python')
    for index,row in dfmap.iterrows():
        drug = row['Inhibitor']
        if drug not in drug2dan:
            logger.error('%s was not found in %s',drug,mappingfilepath)
            sys.exit()
        genes = str(row['Targets']).split(',')
        result = drug2dan[drug]
        for gene in genes:
            hgene2dan[gene].append( (drug,result) )

    logger.info('Read %d compounds and %d targets for the dissemination assay from Georgess et al. (2020).',len(drug2dan),len(hgene2dan))

    return(hgene2dan)

def writeWeights(args,v2weight,v2rank,v2indeg,v2outdeg,v2type,vertices_ranked):
    """
    Writes NetResponse weights to network directory.
    """
    filepath = args.network_dir + '/geneRankings.tsv' 
    fp = open(filepath,'w')
    fields = ['gene','rank','type','weight','indeg','outdeg']
    fp.write('\t'.join(fields) + '\n')
    for v in vertices_ranked:
        toks = [v,str(v2rank[v]),v2type[v],"{:.2e}".format(v2weight[v]),str(v2indeg[v]),str(v2outdeg[v])]
        fp.write('\t'.join(toks) + '\n')
    fp.close()
    logger.info('Wrote ranked genes to %s.',filepath)
    return None

def getRankings(v2weight,ties='high',reverse_order=True):
    """
    Returns the following:
    vertices_ranked: list of genes ranked by v2weight. Default order is descending (reverse), and
    v2rank: a dict of gene names to rank. Default order is descending (reverse).
    """
    pairs = sorted([ (v2weight[v], v) for v in v2weight ], reverse = reverse_order)
    vertices_ranked = [ v for (d,v) in pairs ]
    v2rank = dict()
    if ties == 'high':
        rank = 1
        (d0,v0) = pairs[0]
        v2rank[v0] = rank
        previous_d =d0
        counter = 0
        for (d,v) in pairs:
            if d == previous_d and v != v0:
                #print('Tie',rank)
                v2rank[v] = rank
                counter += 1
            elif d < previous_d:
                rank = rank + counter + 1
                counter = 0
                v2rank[v] = rank
                previous_d = d
            elif d > previous_d:
                logger.error('Pairs are not ranked')
                sys.exit()
    return(vertices_ranked,v2rank)


def getWeights(g_list,v2i,driver,responses,nt,dt):
    """
    Returns the weight for each gene in the network in a dict.
    """
    logger.info('Calculating NetResponse weights for network genes...')
    v2w = dict()
    vertices = sorted(v2i.keys())

    for v in vertices:
        v2w[v] = 0.0
    d = driver
    a = v2i[d]
    for r in responses:
        b = v2i[r]
        for v in vertices:
            x = v2i[v]
            y = [(g_list[nt-i][b,x])*(g_list[i][x,a]) for i in range(nt+1)]
            v2w[v] = v2w[v] + integrate.romb(y,dx=dt)
    norm = 0
    for v in vertices:
        norm += v2w[v]
    if norm == 0.:
        norm = 1.
    for v in vertices:
        v2w[v] = v2w[v] / norm
    return(v2w)

def getGlist(lap_mat,dt,nt):
    """
    Returns the G(t) as a list of matrices.
    """
    logger.info('Calculating the Green\'s function...')
    ret = [None] * (nt + 1)
    (nr,nc) = lap_mat.shape
    assert(nr == nc), 'Bad laplacian shape: %s' % str(lap_mat.shape)
    my_eye = np.eye(nr)
    epsilon = - dt * lap_mat
    g1 = linalg.expm(epsilon)
    curval = my_eye
    ret[0] = curval
    for i in range(nt):
        curval = np.matmul(curval, g1)
        ret[i + 1] = curval
    return(ret)

def getTimeSteps(driver,v2outdeg):
    """
    Returns the time step size and number of time steps.
    """
    #relaxation parameter (the reduction in driver activity)
    r = 0.1

    #number of time steps to reach r (must be a non-negative power of 2)
    nt = 8

    #calculate dt (time step size)
    den = 1 - r
    dt = -log(den,e) / (v2outdeg[driver] * nt)

    #logger.info('dt calculated to be ~%f. Driver density will reach %g in %d steps. Driver outdeg is %d.', round(dt,5),den,nt,v2outdeg[driver])
    return(dt,nt)

def getLaplacian(edge_dict,v2i):
    """
    Returns the adjacency matrix and Laplacian.
    """
    nvert = len(v2i)
    adj_mat = np.zeros( (nvert,nvert) )
    indeg_arr = np.zeros(nvert)
    outdeg_arr = np.zeros(nvert)
    lap_mat = np.zeros( (nvert,nvert) )
    for (src,dest) in edge_dict.keys():
        (a,b) = (v2i[src],v2i[dest])
        adj_mat[b,a] += 1.0
        indeg_arr[b] += 1.0
        outdeg_arr[a] += 1.0
    for a in range(nvert):
        lap_mat[a,a] = outdeg_arr[a]
    lap_mat = lap_mat - adj_mat
    return(adj_mat, lap_mat)


def nameToType(genelist,driver,intermediates,responses):
    """
    Returns a dict of gene names to type.
    D: driver.
    I: intermediate.
    R: response.
    """
    v2type = dict()
    for v in genelist:
        if v == driver:
            v2type[v] = 'D'
        elif v in intermediates:
            v2type[v] = 'I'
        elif v in responses:
            v2type[v] = 'R'
        else:
            logger.warning('Gene %s is not in the driver, intermediate, or response lists.',v)
    return(v2type)

def nameToDegree(genelist,edge_ba,edge_ab):
    """
    Returns indegree and outdegree dicts.
    Gene name is the key.
    """
    v2indeg = dict()
    v2outdeg = dict()
    for v in genelist:
        v2indeg[v] = 0
        if v in edge_ba:
            v2indeg[v] = len(edge_ba[v].keys())
        v2outdeg[v] = 0
        if v in edge_ab:
            v2outdeg[v] = len(edge_ab[v].keys())
    return(v2indeg,v2outdeg)


def nameToIndex(genelist):
    """
    Reads a list of gene names.
    Returns dicts of names to indices and vice versa
    """
    n2i = dict()
    i2n = dict()
    names = sorted(genelist)
    for (n,i) in zip( names, list(range(len(names))) ):
        n2i[n] = i
        i2n[i] = n
    return(n2i,i2n)

def edgeTuplesToDicts(edge_dict):
    """
    Reads dict of edges.
    Returns two dict-of-dicts.
    a is source, be is destination.
    ab uses sources as first key.
    ba uses dests as first key.
    """
    edge_ab = dict()
    edge_ba = dict()
    for (a,b) in edge_dict.keys():
        if a not in edge_ab:
            edge_ab[a] = dict()
        edge_ab[a][b] = True
        if b not in edge_ba:
            edge_ba[b] = dict()
        edge_ba[b][a] = True

    return(edge_ab, edge_ba)

def readEdges(args,verbose=False):
    """
    Reads interactions.
    Returns a dictionary of directed edges.
    Undirected edges are recorded as two directed edges in the dictionary.
    """
    if args.s == 'human':
        filepath = './database/human_interactions.txt'
    elif args.s == 'mouse':
        filepath = './database/mouse_interactions.txt'

    #check filepath exists
    if not os.path.isfile(filepath):
        logger.error('Missing database file %s. Run NetResponse analysis.',filepath)
        sys.exit()

    ret = dict()
    ndirected = 0
    nundirected = 0
    fp = open(filepath, 'r')
    for line in fp:
        my_list = [ ]
        toks = line.strip().split()
        if (len(toks) < 2) or (len(toks) >3):
            if verbose:
                logger.warning('Bad line, wrong token count: %s', line)
            continue
        (a, b) = (toks[0], toks[1])
        if len(toks)==2:
            my_list.append( (a,b) )
            ndirected += 1
        elif len(toks)==3 and toks[2]=='TF':
            my_list.append( (a,b) )
            ndirected += 1
        elif len(toks)==3 and toks[2]=='PPI':
            my_list.append( (a,b) )
            my_list.append( (b,a) )
            nundirected += 1
        else:
            if verbose:
                logger.warning('Bad line, unknown edge type: %s', line)
            continue
        for key in my_list:
            ret[key] = True

    nedge = len(ret)
    if verbose:
        logger.info('Network consists of %d edges: %d directed edges and %d undirected edges from %s.', nedge, ndirected, nundirected, filepath)
    return(ret)

def checkDriver(args,wholegenelist,verbose=False):
    """
    Reads the driver from the commandline.
    Checks if driver is in the network.
    Outputs driver as a string.
    """
    #read driver
    ret = args.driver

    #check if driver is in the network
    if ret in wholegenelist:
        if verbose:
            logger.info('Found driver %s in the network.' % ret)
    else:
        logger.error('Driver %s was not found in the network. Pick another driver.' % ret)
        sys.exit()

    return(ret)

def readDriver(args):
    """
    Reads the driver from the commandline.
    Outputs driver as a string.
    """
    #read driver
    ret = args.driver

    return(ret)

def readResponses(args,driver):
    """
    Reads a file provided by the user containing the gene symbols of the response genes.
    Removes driver from the response genes.
    Returns a list of response genes.
    """

    responseFilepath = args.network_dir + '/responses.xlsx'
    if os.path.isfile(responseFilepath):
        #logger.info('Found response file %s. Delete this file if you want to use a file with a .csv or .tsv extension.',responseFilepath)
        pass
    else:
        responseFilepath = args.network_dir + '/responses.csv'
        if os.path.isfile(responseFilepath):
            #logger.info('Found response file %s. Delete this file if you want to use a file with a .tsv extension.',responseFilepath)
            pass
        else:
            responseFilepath = args.network_dir + '/responses.tsv'
            if os.path.isfile(responseFilepath):
                #logger.info('Found response file %s.',responseFilepath)
                pass
            else:
                logger.erro('No response file found. Must provide a file named responses with .xlsx, .csv, or .tsv extension in directory %s.',args.network_dir)
                sys.exit()

    #check if file is empty
    if not os.stat(responseFilepath).st_size:
        logger.error('The response file %s is empty. It must contain at least one gene symbol.',responseFilepath)
        sys.exit()

    #read responses
    if responseFilepath.split(".")[-1] == 'xlsx':
        xl = pd.ExcelFile(responseFilepath)
        res = len(xl.sheet_names)
        if res != 1:
            logger.error('Response file %s must contain only 1 sheet.',responseFilepath)
            sys.exit()
        responses = pd.read_excel(responseFilepath,header=None)
    elif responseFilepath.split(".")[-1] == 'csv':
        responses = pd.read_csv(responseFilepath,sep=',',header=None)
    else:
        responses = pd.read_csv(responseFilepath,sep='\t',header=None)

    responses = list(set(responses[0].tolist()))

    #remove driver from responses
    if driver in responses:
        responses.remove(driver)

    return(responses)



def getNetworkResponses(args,driver,wholegenelist,verbose=False):
    """
    Reads a file provided by the user containing the gene symbols of the response genes.
    Removes driver from the response genes.
    Checks if each response gene is in the network.
    Outputs list of network response genes called networkResponse.
    Returns a list of network responses.
    """

    #check if network response file exists
    networkResponsesFilepath = args.network_dir + '/networkResponses.txt' 
    if os.path.isfile(networkResponsesFilepath):
        if verbose:
            logger.info('Found network responses file %s. Delete this file if you want the network respones to be compiled again.',networkResponsesFilepath)
        
        #check if file is empty
        if not os.stat(networkResponsesFilepath).st_size:
            logger.error('The network responses file is empty. It must contain at least one gene symbol. Please delete this file and run NetResponse again.')
            sys.exit()

        responses = pd.read_csv(networkResponsesFilepath,sep='\t',header=None)
        responses = list(set(responses[0].tolist()))

        #remove driver from network responses
        if driver in responses:
            responses.remove(driver)
            if verbose:
                logger.warning('Found driver %s in the network response file. It was removed.' % (driver))
    
        #check if responses are in network and output list
        ret = list()
        for response in responses:
            if response in wholegenelist:
                ret.append(response)
            else:
                logger.error('Not all network response genes in %s are in the network. Please delete this file and run NetResponse analysis, again.', networkResponsesFilepath)
                sys.exit()

    else:
        #check if user response file exists
        responseFilepath = args.network_dir + '/responses.xlsx'
        if os.path.isfile(responseFilepath):
            #logger.info('Found response file %s. Delete this file if you want to use a file with a .csv or .tsv extension.',responseFilepath)
            pass
        else:
            responseFilepath = args.network_dir + '/responses.csv'
            if os.path.isfile(responseFilepath):
                #logger.info('Found response file %s. Delete this file if you want to use a file with a .tsv extension.',responseFilepath)
                pass
            else:
                responseFilepath = args.network_dir + '/responses.tsv'
                if os.path.isfile(responseFilepath):
                    #logger.info('Found response file %s.',responseFilepath)
                    pass
                else:
                    logger.erro('No response file found. Must provide a file named responses with .xlsx, .csv, or .tsv extension in directory %s.',args.network_dir)
                    sys.exit()

        #check if file is empty
        if not os.stat(responseFilepath).st_size:
            logger.error('Response file %s is empty. Must contain at least one gene symbol.' % responseFilepath)
            sys.exit()
        
        #read responses and convert to list
        if responseFilepath.split(".")[-1] == 'xlsx':
            xl = pd.ExcelFile(responseFilepath)
            res = len(xl.sheet_names)
            if res != 1:
                logger.error('Response file %s must contain only 1 sheet.',responseFilepath)
                sys.exit()
            responses = pd.read_excel(responseFilepath,header=None)
        elif responseFilepath.split(".")[-1] == 'csv':
            responses = pd.read_csv(responseFilepath,sep=',',header=None)
        else:
            responses = pd.read_csv(responseFilepath,sep='\t',header=None)
        
        responses = list(set(responses[0].tolist()))
        if verbose:
            logger.info('Read %s response genes from %s.' % (str(len(responses)),responseFilepath))
    
        #remove driver from responses
        if driver in responses:
            responses.remove(driver)
            if verbose:
                logger.info('Found driver %s in the response list. It was removed. %s response genes remain.' % (driver,str(len(responses))))
        else:
            if verbose:
                logger.info('Did not find driver %s in the response list.' % driver)
    
        #check if responses are in network and output list
        ret = list()
        outfilepath = args.network_dir + '/networkResponses.txt' 
        fp = open(outfilepath, 'w')
        for response in responses:
            if response in wholegenelist:
                ret.append(response)
                fp.write(response + '\n')
        fp.close()
        if verbose:
            logger.info('Found %s of %s response genes in the network. Wrote to %s.' % (str(len(ret)),str(len(responses)),outfilepath))
    return(ret)

def readInteractions(filepath):
    """
    Reads a NetResponse Interaction file.
    Returns a pandas dataframe
    """
    ret = pd.read_csv(filepath, sep='\t',header=None)
    return(ret)

def readGenes(filepath):
    """
    Read a NetResponse Gene file.
    Returns a list
    """
    ret = pd.read_csv(filepath,sep='\t',header=None)
    ret = list(set(ret[0].tolist()))
    return(ret)

def readIIDMouse():
    """
    Reads mouse ppi file database/mouse_annotated_PPIs.txt.
    Converts ppi from uniprot id to mgi symbol.
    Returns a list of unique ppi interactions as tuple (a,b).
    a and b are the genes of the two interacting proteins.
    Removes unwanted ppi between a protein and itself.
    """
    filename = './database/mouse_annotated_PPIs.txt'

    #get mapping from mouse uniprot id to mgi symbol
    uniprot2mgisymbol = get_uniprot2mgisymbol()

    #read the ppi file
    #read_iid_ppi
    logger.info('Reading mouse PPI from %s. Removing PPI between a protein and itself...', filename)
    reftable = pd.read_csv(filename,sep='\t',low_memory=False)
    ret = dict()
    for index,row in reftable.iterrows():
        (a,b) = ( row['uniprot1'], row['uniprot2'] )
        if a == b:
            continue
        if a in uniprot2mgisymbol.keys() and b in uniprot2mgisymbol.keys():
            for genea in uniprot2mgisymbol[a]:
                for geneb in uniprot2mgisymbol[b]:
                    if (genea,geneb) not in ret and (geneb,genea) not in ret:
                        ret[(genea,geneb)] = True
    logger.info('Mapped %d mouse PPI interactions from uniprot id to mgi symbol.', len(ret))
    return(list(ret.keys()))

def get_uniprot2mgisymbol():
    """
    Returns a dict mapping from mouse uniprot id to mgi symbol.
    """
    #get uniprot to mgi num mapping
    uniprot2mginum = get_uniprot2mginum()
    #get mgi num to mgi symbol mapping
    mginum2symbol = get_mginum2symbol()
    ret = dict()
    counter = 0
    for prot in uniprot2mginum.keys():
        ret[prot] = dict()
        for num in uniprot2mginum[prot]:
            if num in mginum2symbol:
                ret[prot][mginum2symbol[num]] = True
                counter += 1
    logger.info('Mapped %d mouse uniprot ids to %d mgi symbols.', len(ret), counter)
    
    return(ret)

def get_mginum2symbol():
    """
    Reads database/MRK_List2.rpt'
    Returns a dict mapping from mouse mgi num to symbol.
    """
    mginum2symbol_file = './database/MRK_List2.rpt'
    logger.info('Reading mgi marker accession number to mgi marker symbol from %s...', mginum2symbol_file)
    reftable = pd.read_csv(mginum2symbol_file,sep='\t')
    ret = dict()
    for index,row in reftable.iterrows():
        num = row['MGI Accession ID']
        symbol = row['Marker Symbol']
        ret[num] = symbol
    logger.info('Mapped %d mgi marker accession numbers to symbols.',len(ret))
    return(ret)

def get_uniprot2mginum():
    """
    Reads database/MOUSE_10090_idmapping.dat.
    Returns a dict mapping from mouse uniprot id to MGI number.
    """
    uniprot2mginum_file = './database/MOUSE_10090_idmapping.dat'
    logger.info('Reading uniprot to mgi marker accession number from %s...', uniprot2mginum_file)
    ret = dict()
    ret2 = dict()
    reftable = pd.read_csv(uniprot2mginum_file, sep='\t')
    for index,row in reftable.iterrows():
        name = row[0]
        if name not in ret:
            ret[name] = dict()
        id_type = row[1]
        if id_type == 'MGI':
            mgi = row[2]
            ret[name][mgi] = True
    for prot in ret.keys():
        if ret[prot]:
            ret2[prot] = ret[prot]
    logger.info('Mapped %d of %d uniprot ids to a MGI number.', len(ret2), len(ret))
    return(ret2)

def readIIDHuman():
    """
    Reads ppi from Human IID file in database dir.
    Returns list of unique ppi interactions as tuple (src,des).
    Removes unwanted PPI between a protein and itself.
    """
    filename = './database/human_annotated_PPIs.txt'
    logger.info('Reading human PPI from %s. Removing PPI between a protein and itself...', filename)
    ref = pd.read_csv(filename, sep='\t',low_memory=False)

    ret = dict()

    for index,row in ref.iterrows():
        (a,b) = (row['symbol1'],row['symbol2'])
        if len(a.strip().split()) != 1:
            continue
        if len(b.strip().split()) != 1:
            continue
        if a == b:
            continue
        if (a,b) not in ret and (b,a) not in ret:
            ret[(a,b)] = True

    logger.info('Read %d unique human PPI.', len(ret))
    return(list(ret.keys()))

def readTF(species):
    """
    Reads regulator interactions from files with format: GENE1  GENE2   Mode    PMID.
    Files stored in database dir. 
    File should not have a header.
    """
    if species == 'human': 
        filename = './database/trrust_rawdata.human.tsv'
    elif species == 'mouse':
        filename = './database/trrust_rawdata.mouse.tsv'
    else:
        logger.error('%s is not a supported species.',species)
        sys.exit()

    ret = dict()
    fp = open(filename,'r')
    for line in fp:
        toks = line.strip().split('\t')
        if len(toks)!=4:
            logger.warning('Bad line: wrong token count %s',line)
            continue
        (genea,geneb) = (toks[0],toks[1])
        ret[(genea,geneb)] = True
    logger.info('Read %d %s TF regulatory interactions.',len(ret),species)
    return(list(ret.keys()))

def main():
    args = CMDParser()
    args.func(args)
    return None

if __name__=='__main__':
    main()
