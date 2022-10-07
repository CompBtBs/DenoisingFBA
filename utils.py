#load libraries
import scanpy as sc
import pandas as pd 
from sklearn import metrics 
import numpy as np
import pandas as pd
import re
from progressbar import ProgressBar,Bar,Percentage
from scanpy import AnnData
import itertools
import cobra as cb
import matplotlib.pyplot as plt
import scipy

cb.Configuration.solver="glpk"
cb.Configuration().tolerance=1E-07


"""
Class to compute the RAS values

"""

class RAS_computation:

    def __init__(self,adata,model):
                                                       
        self._logic_operators = ['and', 'or', '(', ')']
        self.val_nan = np.nan

        # Build the dictionary for the GPRs
        df_reactions = pd.DataFrame(index=[reaction.id for reaction in model.reactions])
        gene_rules=[reaction.gene_reaction_rule for reaction in model.reactions]
        
        
        gene_rules=[el.replace("OR","or").replace("AND","and").replace("(","( ").replace(")"," )") for el in gene_rules]        
        df_reactions['rule'] = gene_rules
        df_reactions = df_reactions.reset_index()
        df_reactions = df_reactions.groupby('rule').agg(lambda x: sorted(list(x)))
        
        self.dict_rule_reactions = df_reactions.to_dict()['index']

        # build useful structures for RAS computation
        self.model = model
        self.count_adata = adata.copy()
        self.genes = self.count_adata.var.index.intersection([gene.id for gene in model.genes])
        
        
        self.cell_ids = list(self.count_adata.obs.index.values)
        self.count_df_filtered = self.count_adata.to_df().T.loc[self.genes]
 
    def compute(self):

        self.or_function = np.nansum
        self.and_function = np.nanmin
        regexp=re.compile(r"\([a-zA-Z0-9-.\s]+\)")  # regular expression inside a parenthesis
        
        ras_df = pd.DataFrame(index=range(len(self.dict_rule_reactions)), columns=self.cell_ids)
        ras_df[:][:] = self.val_nan

        pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=len(self.dict_rule_reactions)).start()
        i = 0
        
        # for loop on reactions
        ind = 0       
        for rule, reaction_ids in self.dict_rule_reactions.items():
            if len(rule) != 0:
                # there is one gene at least in the formula
                rule_split = rule.split()
                rule_split_elements = list(filter(lambda x: x not in self._logic_operators, rule_split))  # remove of all logical operators
                rule_split_elements = list(np.unique(rule_split_elements))                                # genes in formula
                
                # which genes are in the count matrix?                
                genes_in_count_matrix = list(set([el for el in rule_split_elements if el in self.genes]))
                genes_notin_count_matrix = list(set([el for el in rule_split_elements if el not in self.genes]))


                if len(genes_in_count_matrix) > 0: #there is at least one gene in the count matrix
                     if len(rule_split) == 1:
                         #one gene --> one reaction
                         ras_df.iloc[ind] = self.count_df_filtered.loc[genes_in_count_matrix]
                     else:                        
                        # more genes in the formula
                        lista = re.findall(regexp, rule)
                        if len(lista) == 0:
                             #or/and sequence
                             matrix = self.count_df_filtered.loc[genes_in_count_matrix].values
                             if len(genes_notin_count_matrix) > 0:
                                matrix = np.vstack([matrix, [self.val_nan for el in self.cell_ids]])

                             if 'or' in rule_split: 
                                ras_df.iloc[ind] = self.or_function(matrix, axis=0)
                             else:
                                ras_df.iloc[ind] = self.and_function(matrix, axis=0)
                        else:
                            # ho almeno una tonda
                            data = self.count_df_filtered.loc[genes_in_count_matrix]  # dataframe of genes in the GPRs
                            genes = data.index
                            j = 0
                             
                            for cellid in self.cell_ids:    #for loop on the cells
                                lista_cell = lista.copy()
                                rule_cell = rule
                                 
                                while len(lista_cell) > 0:
                                    #
                                    for el in lista_cell:
                                        #print(el[1:-1])
                                        value = self._evaluate_expression(el[1:-1].split(), data[cellid], genes)
                                        rule_cell = rule_cell.replace(el, str(value))   
                                    lista_cell = re.findall(regexp, rule_cell)      
         
                                ras_df.iloc[ind, j] = self._evaluate_expression(rule_cell.split(), data[cellid], genes)
                                j=j+1
      
            ind = ind+1
            #update percentage
            pbar.update(i+1)
            i = i+1

        pbar.finish()
        
        ras_df=ras_df.astype("float")    
        ras_df['REACTIONS'] = [reaction_ids for rule,reaction_ids in self.dict_rule_reactions.items()]
        
        reactions_common = pd.DataFrame()
        reactions_common["REACTIONS"] = ras_df['REACTIONS']
        reactions_common["proof2"] = ras_df['REACTIONS']
        reactions_common = reactions_common.explode('REACTIONS')
        reactions_common = reactions_common.set_index("REACTIONS")

        ras_df = ras_df.explode("REACTIONS")
        ras_df = ras_df.set_index("REACTIONS")

        # Drop na rules
        ras_df=ras_df.dropna()        
        
        #create AnnData structure for RAS
        ras_adata = AnnData(ras_df.T)

        #add metadata
        for el in self.count_adata.obs.columns:
            ras_adata.obs["countmatrix_"+el]=self.count_adata.obs[el]


        return ras_adata


    def _check_number(self,value):
      try:
        float(value)
        return True
      except ValueError:
        return False

    def _evaluate_expression(self, rule_split, values_cellid, genes):
        
        #ci sono per forza solo or
        rule_split2 = list(filter(lambda x: x != "or" and x!="and", rule_split))   

        values = list()
        i=0
        for el in rule_split2:
             if self._check_number(el):
                 values.append(float(el))
             elif el in genes:
                 values.append(values_cellid[el])
             else:
                 values.append(self.val_nan)
                 i=i+1
                 
        if i==len(rule_split2):
            return self.val_nan
        if "or" in rule_split:
            #or sequence
            return self.or_function(values)
        else:
            #and sequence
            return self.and_function(values)


""" Funcion to perform single cell FBA"""

def scFBA(model,ras_adata,dfFVA,eps=0,verbose=False):


    each_count=20
    reactions=[reaction.id for reaction in model.reactions]
    
    #normalize ras matrix
    ras_matrix=ras_adata.to_df().T
    ras_matrix=ras_matrix.T.div(ras_matrix.T.max(axis=0)).T
    ras_matrix.fillna(0,inplace=True)

    indexes=ras_matrix.index
    cells=list(ras_matrix.columns)    

    i=0
    valori=[]

    dfOpt=pd.DataFrame(index=reactions,columns=cells)

    
    for cell in cells:
        model2=model.copy()
        for reaction in reactions:
            if reaction in indexes:
                lower_bound=dfFVA.loc[reaction,"minimum"]
                upper_bound=dfFVA.loc[reaction,"maximum"]

                if dfFVA.loc[reaction,"maxABS"]>0:

                    valMax=eps+(dfFVA.loc[reaction,"maximum"]-eps)*ras_matrix.loc[reaction,cell]
                    valMin=-eps+(dfFVA.loc[reaction,"minimum"]+eps)*ras_matrix.loc[reaction,cell]                

                    if upper_bound>0 and lower_bound==0:
                        model2.reactions.get_by_id(reaction).upper_bound=valMax #

                    if upper_bound==0 and lower_bound<0:
                        model2.reactions.get_by_id(reaction).lower_bound=valMin

                    if upper_bound!=0 and lower_bound!=0: 
                        model2.reactions.get_by_id(reaction).lower_bound=valMin
                        model2.reactions.get_by_id(reaction).upper_bound=valMax

            model2.solver.configuration.timeout=5
            opt=model2.optimize()
            dfOpt[cell]=opt.fluxes
            if verbose:
                print(opt.objective_value)

        if not verbose and i % each_count==0:
            print(i,cell)   
        i=i+1

            
    flux_adata=AnnData(dfOpt.T.round(10))
    flux_adata.obs=ras_adata.obs
    
    return flux_adata


def find_essential(model):
    list_essential_reactions=cb.flux_analysis.find_essential_reactions(model)
    list_essential_reactions=[reaction.id for reaction in list_essential_reactions]
    return list_essential_reactions


def plot_distributions(adatasets):
    
    i=0
    for adata in adatasets:
        if i==1: 
            groupby="countmatrix_Factor Value[disease]"
        else:
            groupby="countmatrix_Type"

        len_reactions=adata.to_df().shape[1]
        adata.obs["perc_zeros"]=(adata.to_df()==0).sum(1)/len_reactions*100
        axes=sc.pl.violin(adata,keys=["perc_zeros"],groupby=groupby,stripplot=False,show=False,scale="area")
        axes.set_ylim([0,100])
        axes.grid()
        axes.set_xlabel("")
        axes.set_ylabel("% of zero RAS values")
        i=i+1
        
def plot_correlations(adatasets):
    
    i=0
    for adata in adatasets:
      
        df=adata.to_df()
        dfCorr=df.corr(method='spearman')
        dfCorr=dfCorr.where(np.triu(np.ones(dfCorr.shape),k=1).astype(np.bool))
        
        plt.figure()
        plt.hist(dfCorr.values.ravel()[dfCorr.values.ravel()>-100000000])
        plt.xlim([-1,1])
        plt.grid()
        plt.xlabel("")
        plt.ylabel("N° of reaction pairs")
        
        
def plot_two_reaction_correlation(adatasets,reactions):
    i=0
    for adata in adatasets:
        if i==1: 
            groupby="countmatrix_Factor Value[disease]"
        else:
            groupby="countmatrix_Type" 

        df=adata.to_df()
        df["type"]=adata.obs[groupby].values

        names=list(set(adata.obs[groupby].values))

        df1=df.loc[:,reactions[0]]
        df2=df.loc[:,reactions[1]]

        names1=df["type"][df["type"]==names[0]].index
        names2=df["type"][df["type"]==names[1]].index

        r1, p = scipy.stats.pearsonr(df1.values,df2.values)

        plt.figure()
        plt.scatter(df1.loc[names1].values,df2.loc[names1].values)
        plt.scatter(df1.loc[names2].values,df2.loc[names2].values)
        plt.xlabel(reactions[0])
        plt.ylabel(reactions[1])
        plt.grid()
        plt.title("".join(['r: ',str(np.round(r1,2)),", ",'p: ',str(np.round(p,4)) ]))

        i=i+1
        
      
def table_sparsity(adatasets_countmatrix,adatasets_rasmatrix,names_datasets):
    
    dfTable=pd.DataFrame(index=names_datasets,
                   columns=["count_matrix","metabolic_countmatrix",
                            "ras_matrix","essential ras_matrix"])
    
    df_genes=pd.read_csv("data/genes_ENGRO2.csv",index_col=0)
    dfFVA=pd.read_csv("data/FVA.csv",index_col=0)["essential"]
    essential=dfFVA[dfFVA==True].index
    
    for count_adata,ras_adata,name in zip(adatasets_countmatrix,
                                          adatasets_rasmatrix,
                                          names_datasets):
        
         
        
        cells,len_reactions=ras_adata.to_df().shape
        cells,genes=count_adata.to_df().shape
        
        #count matrix
        dfTable.loc[name,"count_matrix"]=str(np.round((count_adata.to_df()==0).sum(0).sum(0)/(genes*cells)*100,2))+"%"
        
        #metabolic count matrix
        if name=="datasetE-GEOD-86618":
            met_genes=list(df_genes["metabolic_genes_ensg"].values)
        else:
            met_genes=list(df_genes["metabolic_genes_genesymbol"].values) 
        met_genes=[el for el in met_genes if el in count_adata.var.index]
        
        dfTable.loc[name,"metabolic_countmatrix"]=str(np.round((count_adata[:,met_genes].to_df()==0).sum(0).sum(0)/(len(met_genes)*cells)*100,2))+"%"

        #ras matrix
        dfTable.loc[name,"ras_matrix"]=str(np.round((ras_adata.to_df()==0).sum(0).sum(0)/(len_reactions*cells)*100,2))+"%"

        #essential ras values
        reactions=list(ras_adata.var.index)
        essential_specific=[el for el in essential if el in reactions]
        val=(ras_adata[:,essential_specific].to_df()==0).sum(0).sum(0)

        dfTable.loc[name,"essential ras_matrix"]=str(np.round(val/(len(essential_specific)*cells)*100,2))+"%"

        
    return dfTable



def table_sparsity_denoised(adatasets_magic,adatasets_enhance,adatasets_saver,names_datasets):
    
    dfTable=pd.DataFrame(index=names_datasets,
                   columns=["magic","enhance","saver"])
    

    dfFVA=pd.read_csv("data/FVA.csv",index_col=0)["essential"]
    essential=dfFVA[dfFVA==True].index
    
    for ras_adata_magic,ras_adata_enhance,ras_adata_saver,name in zip(adatasets_magic,
                                          adatasets_enhance,
                                                      adatasets_saver,
                                          names_datasets):
        

        cells,len_reactions=ras_adata_magic.to_df().shape

        reactions=list(ras_adata_magic.var.index)
        essential_specific=[el for el in essential if el in reactions]        
        
        #ras matrix magic
        val=np.round((ras_adata_magic.to_df()==0).sum().sum()/(len_reactions*cells)*100,2)
        val2=np.round((ras_adata_magic[:,essential_specific].to_df()==0).sum().sum()/(len(essential_specific)*cells)*100,2)
        
        dfTable.loc[name,"magic"]=str(val)+"%("+str(val2)+"%)"
        
        #ras matrix enhance
        val=np.round((ras_adata_enhance.to_df()==0).sum().sum()/(len_reactions*cells)*100,2)
        val2=np.round((ras_adata_enhance[:,essential_specific].to_df()==0).sum().sum()/(len(essential_specific)*cells)*100,2)
        dfTable.loc[name,"enhance"]=str(val)+"%("+str(val2)+"%)"

        #ras matrix enhance
        val=np.round((ras_adata_saver.to_df()==0).sum().sum()/(len_reactions*cells)*100,2)
        val2=np.round((ras_adata_saver[:,essential_specific].to_df()==0).sum().sum()/(len(essential_specific)*cells)*100,2)
        dfTable.loc[name,"saver"]=str(val)+"%("+str(val2)+"%)"
        
    return dfTable


def find_bh(
        ras_adata,                        
        n_pcs=[10],                        
        n_neighbors=[5],                    
        ):

    ras_adata_clustering=ras_adata.copy()

    
    #%%pca analysis
    sc.tl.pca(ras_adata_clustering, svd_solver='arpack',n_comps=max(n_pcs))

    cluster_values_sil=list()

    num_cluster=list()
    

    pcs_values=list()
    neigh_values=list()

    for npc in n_pcs:          
        for neigh in n_neighbors:
            adata=ras_adata_clustering.copy()
            sc.pp.neighbors(adata, n_neighbors=neigh, n_pcs=npc)   
            sc.tl.leiden(adata,key_added = "leiden")

            pcs_values.append(npc)
            neigh_values.append(neigh)
            if len(set(adata.obs["leiden"].values))>1:
                cluster_values_sil.append(metrics.silhouette_score(adata.obsm['X_pca'][:,0:npc],adata.obs["leiden"],metric='euclidean'))
            else:
                cluster_values_sil.append(None)

            num_cluster.append(len(set(adata.obs["leiden"].values)))                   

    df=pd.DataFrame()

    df["pcs_values"]=pcs_values
    df["neigh_values"]=neigh_values
    df["num_cluster"]=num_cluster
    df["cluster_values_sil"]=cluster_values_sil

    return df
