#!/usr/bin/env python3
from argparse import ArgumentParser

#parsing input
parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="Input",default=[], nargs='+', required=True, help='One or multiple cooler files')
parser.add_argument("-l", "--list", dest="List", required=True, help='Sorted TAB list with no header and  Chr Start End Type_1_clusters Type_2_clusters colums' )
parser.add_argument("-p", "--processors",dest="processors", type=int, default=True,help='number of cores to use, if not specified it uses all the processors available')
parser.add_argument("-o", "--output", dest="outputDir",default="Results", help='output directory')
parser.add_argument("-f", "--feature_name", dest="element",default="Element", help='Name of the genomic feature')
parser.add_argument("-b", "--boundaries", dest="boundaries",action='store_true',default= False, help='calculates interaction numbers over consecutive elements (the List provided must be sorted)')
parser.add_argument("-pl", "--plot", dest="plot",action='store_true',default= False, help='activates exploratory plotting options')
parser.add_argument("-s", "--sex_Chr", dest="sex",action='store_true',default= False, help='Keep sexual chromosomes')
parser.add_argument("-X", "--ChrX", dest="chrX",action='store_true',default= False, help='Keep X chromosome')
parser.add_argument("-coo", "--Coordinates", dest="coo",action='store_true',default= False, help='Recover Feature number with genomic positions')
parser.add_argument("-trans", "--incl_trans", dest="trans",action='store_true',default= False, help='Include in trans interactions, False by default')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
options= parser.parse_args()


import pandas as pd
import cooler
import numpy as np
from multiprocessing import Pool, cpu_count
from functools import partial
from itertools import repeat
import os
if options.plot:
    import seaborn as sns
    import matplotlib.pyplot as plt
    plt.switch_backend('agg')
# recovers data from the cooler files
def run (Data,List,Name,Element_Name,Trans,RANGE):
    c=[]
    # format coordinates for cooler files
    List["Start"]=round(List["Start"]/Data.info['bin-size']).astype(int)
    List["Stop"]=round(List["Stop"]/Data.info['bin-size']).astype(int)
    
    for i in range(RANGE[0],RANGE[1]):
        chr1=Data.chromnames[i]
        if Trans:
            for n in range(0,RANGE[2]):
                chr2=Data.chromnames[n]
                #load the matrix
                matrix = Data.matrix(sparse=True).fetch(chr1,chr2).tocsr()
                if chr1==chr2:
                    #remouves the diagonal
                    matrix.setdiag(0)
                    mem='Cis'
                    CisTransIntra='Cis'
                else:
                    mem='Trans'
                    CisTransIntra='Trans'
                
                l1 = List[List["Chr"]==chr1]
                l2 = List[List["Chr"]==chr2]
                for i1 in range(l1.index[0],l1.index[-1]+1):
                    sl1 = l1.loc[i1]
                    for i2 in range(l2.index[0],l2.index[-1]+1):
                        sl2 = l2.loc[i2]
                        v = np.nan_to_num(matrix[sl1[1]:sl1[2], sl2[1]:sl2[2]].toarray()).sum()
                        if all(sl1==sl2):
                            CisTransIntra='Intra'
                        
                        c.append([Name,'{}-{}'.format(Element_Name,i1),'{}-{}'.format(Element_Name,i2),v,('{}to{}'.format(sl1['Type_1'],sl2['Type_1'])),('{}to{}'.format(sl1['Type_2'],sl2['Type_2'])),CisTransIntra])
                        CisTransIntra=mem
        else:
            #load the matrix
            matrix = Data.matrix(sparse=True).fetch(chr1).tocsr()
            #remouves the diagonal
            matrix.setdiag(0)
            mem='Cis'
            CisTransIntra='Cis'
            l1 = List[List["Chr"]==chr1]
            for i1 in range(l1.index[0],l1.index[-1]+1):
                sl1 = l1.loc[i1]
                for i2 in range(l1.index[0],l1.index[-1]+1):
                    sl2 = l1.loc[i2]
                    v = np.nan_to_num(matrix[sl1[1]:sl1[2], sl2[1]:sl2[2]].toarray()).sum()
                    if all(sl1==sl2):
                        CisTransIntra='Intra'
                        
                    c.append([Name,'{}-{}'.format(Element_Name,i1),'{}-{}'.format(Element_Name,i2),v,('{}to{}'.format(sl1['Type_1'],sl2['Type_1'])),('{}to{}'.format(sl1['Type_2'],sl2['Type_2'])),CisTransIntra])
                    CisTransIntra=mem

    return pd.DataFrame(c,columns=['Name','Anchored_Element','Interactiong_Element','Interactions_number','Type_1','Type_2','CisTransIntra'])

#checks directory and creates it if not existing
def check_dir(path):
    directory = os.path.dirname(path)
    if not os.path.exists(directory):
        os.makedirs(directory)
def boundaries_count(Result):
    Boundaries=Result[Result.apply(lambda x: (int(x['Anchored_Element'].split('-')[-1])-int(x['Interactiong_Element'].split('-')[-1]))== -1, axis=1)]
    return Boundaries
#plot function 
def boundaries_plot(Boundaries,outputDir,Lname):
    g=sns.factorplot(
        data=Boundaries,
        x='Name',
        y='Interactions_number',
        palette = 'Set1',
        sharey = False
        )
    g.set_xticklabels(rotation=30)
    g.savefig("{}Boundaries_Strengh_{}.pdf".format(outputDir,Lname))
    g=sns.factorplot(
        data=Boundaries,
        x='Name',
        y='Interactions_number',
        col='Type_1',
        kind='box',
        palette = 'Set1',
        sharey = False
    )
    g.set_xticklabels(rotation=30)
    g.savefig("{}Boundaries_Strengh_{}_split_Type_1_feature.pdf".format(outputDir,Lname))
    g=sns.factorplot(
        data=Boundaries,
        x='Name',
        y='Interactions_number',
        col='Type_2',
        kind='box',
        size=5,
        palette = 'Set1',
        sharey = False
    )
    g.set_xticklabels(rotation=30)
    g.savefig("{}Boundaries_Strengh_{}_split_Type_2_feature.pdf".format(outputDir,Lname))

#plot function 
def Feat_Cis_Trans_intra_plot(Feat_Cis_Trans_intra,outputDir,Lname):
    g=sns.factorplot(
        data=Feat_Cis_Trans_intra,
        x='Feature_genotype',
        y='Interactions_number',
        hue='Name',
        row='CisTransIntra',
        col='Type',
        kind='box',
        palette = 'Set1',
        sharey = False
        )
    g.set_xticklabels(rotation=30)
    g.savefig("{}results_{}_split_by_feature.pdf".format(outputDir,Lname))

#plot function 
def Feat_plot(Feat,outputDir,Lname):
    g=sns.factorplot(
        data=Feat,
        x='Feature_genotype',
        y='Interactions_number',
        hue='Name',
        col='Type',
        kind='box',
        palette = 'Set1',
        sharey = False
        )
    g.set_xticklabels(rotation=30)
    g.savefig("{}results_{}_split_by_feature_no_CIS_TRANS_INTRA.pdf".format(outputDir,Lname))

#plot function 
def Cis_Trans_intra_plot(Cis_Trans_intra,outputDir,Lname):
    g=sns.factorplot(
        data=Cis_Trans_intra,
        x='Name',
        y='Interactions_number',
        col='CisTransIntra',
        kind='box',
        palette = 'Set1',
        sharey = False
    )
    g.set_xticklabels(rotation=30)
    g.savefig("{}results_{}_Cis_Trans_intra.pdf".format(outputDir,Lname))
    
def main(options):

    Input = [os.path.join(os.getcwd(), path) for path in options.Input]
    listin = os.path.join(os.getcwd(),options.List)
    outputDir = '{}/'.format(os.path.join(os.getcwd(),options.outputDir))

    if options.processors:
        processors = cpu_count()-1
    else:
        processors = options.processors

    #create directory fore sults
    check_dir(outputDir)
            
    List=pd.read_table(listin, sep="\t",names=("Chr","Start","Stop","Type_1","Type_2"))
    List=List.sort_values(["Chr","Start","Stop"])
    Lname=listin.split('/')[-1].split(".")[0]
    if options.coo==True:
        List2=List
        List2['element']=options.element+'-'+ List2.index.astype(str)
        List2.to_csv('{}featuren_number_{}.txt'.format(outputDir,Lname), sep='\t', index=False)

    # Names of the samples based on cooler files
    Name = [(i.split('/')[-1].split(".")[0].split("_")[0])for i in Input]
    # Import cooler data
    DATA=[(cooler.Cooler(i))for i in Input]
            
    # split the job in parallel according to the number of processors and the lenght of the list
    if options.sex:
        Chrs = len(DATA[0].chromnames)
    else:
        if options.chrX:
            Chrs = len(DATA[0].chromnames)-1
        else:
            Chrs = len(DATA[0].chromnames)-2
        
    he = [int(round(Chrs*i/processors)) for i in range(1, (processors+1))] 
    hs = [int(round(Chrs*i/processors)) for i in range(0, processors)]
    p = Pool(processors)
    
    # first result
    Result=pd.concat([(pd.concat(p.map(partial(run, DATA[f], List, Name[f],options.element,options.trans),zip(hs,he,repeat(Chrs)))))for f in range(len(DATA))])

    
    Result.to_csv('{}results_{}.txt'.format(outputDir,Lname), sep='\t', index=False)
    
    # create table for boundaries strengh, list file must be ordered by genomic position
    if options.boundaries:
        
        Boundaries=pd.concat(p.map(boundaries_count,[(Result[(Result['Name']==i) & (Result['CisTransIntra']=='Cis')])for i in Name]))
        Boundaries.to_csv('{}Boundaries_Strengh_{}.txt'.format(outputDir,Lname), sep='\t', index=False)
        
        #plot
        if options.plot:
            boundaries_plot(Boundaries,outputDir,Lname)
        
        del Boundaries
        
    #calculate interactions divided for list features (reshape dataframe)
    Type_1=Result[['Name','Anchored_Element','Interactions_number','Type_1','CisTransIntra']]
    Type_2=Result[['Name','Anchored_Element','Interactions_number','Type_2','CisTransIntra']]
    Type_1=Type_1.rename(columns = {'Type_1':'Type'})
    Type_2=Type_2.rename(columns = {'Type_2':'Type'})
    Type_1.loc[:,'Feature_genotype']='Type_1'
    Type_2.loc[:,'Feature_genotype']='Type_2'

    # no need for the origina result, reuse variable 
    Result=pd.concat([Type_1,Type_2])
    Result=Result.groupby(['Name','Anchored_Element','Type','CisTransIntra','Feature_genotype']).aggregate(sum).reset_index()
    
    #print the table devided by Feature and CisTransIntra
    Result.to_csv('{}results_{}_split_by_feature.txt'.format(outputDir,Lname), sep='\t', index=False)
    
    #plot
    if options.plot:
        Feat_Cis_Trans_intra_plot(Result,outputDir,Lname)
    #no need for the Type_2 feature later
    del Type_2
    
    # delete cis trans intra information and sum the interactions devided by the list features
    del Result['CisTransIntra']
    Result=Result.groupby(['Name','Anchored_Element','Type','Feature_genotype']).aggregate(sum).reset_index()
    
    #print the table devided by Feature
    Result.to_csv('{}results_{}_split_by_feature_no_CIS_TRANS_INTRA.txt'.format(outputDir,Lname), sep='\t', index=False)
    #plot
    if options.plot:
        Feat_plot(Result,outputDir,Lname)
    del Result
    
    #caluclate Interactions on the base of the CisTransIntra
    del Type_1['Type']
    del Type_1['Feature_genotype']
    Type_1=Type_1.groupby(['Name','Anchored_Element','CisTransIntra']).aggregate(sum).reset_index()

    # write file 
    Type_1.to_csv('{}results_{}_Cis_Trans_intra.txt'.format(outputDir,Lname), sep='\t', index=False)
        #plot
    if options.plot:
        Cis_Trans_intra_plot(Type_1,outputDir,Lname)

if __name__ == "__main__":
    main(options)

