"""
Created on Wed Jan 23 17:57:51 2019

@author: Haohe Liu
"""
import networkx as nx
import csv_read as cr
import numpy as np
import scipy.sparse as sp
import pickle

class data_cache():
    def __init__(self):
        print("Initializing training data cache...")
        self.drug2num_dict = {}
        self.num2drug_dict = {}
        self.gene2num_dict = {}
        self.num2gene_dict = {}

        self.gene_list = []
        self.gene_drug_list = []

        self.drug_drug_dict = {}
        self.drug_drug_adj_list = []

        self.n_genes = 0
        self.n_drugs = 0

        self.gene_drug_adj_mat = None
        self.gene_graph = nx.Graph()
        self.lost_drug_count = 0

        #initialize gene2num & num2gene
        self.load_gene_dict()
        # initialize drug2num & num2drug
        self.load_drug_dict()
        self.update_length()
        #load adjacent matrixes
        self.gene_net_load()
        self.drug_drug_adj_load()
        self.gene_drug_adj_load()



    def drug_update(self,drug):
        if (drug not in self.drug2num_dict.keys()):
            self.drug2num_dict[drug] = len(self.drug2num_dict.keys())
            self.num2drug_dict[self.drug2num_dict[drug]] = drug
    def gene_update(self,gene):
        if (gene not in self.gene2num_dict.keys()):
            self.gene2num_dict[gene] = len(self.gene2num_dict.keys())
            self.num2gene_dict[self.gene2num_dict[gene]] = gene

    def load_drug_dict(self,fname = './bio-data/bio-decagon-mono.csv'):
        print("loading drug data")
        drug_name,_,_ = cr.csv_read(fname,line = 0)
        for each in drug_name:
            self.drug_update(each[0])
        drug_rel, _, _ = cr.csv_read("./bio-data/bio-decagon-combo.csv", line=0)
        for each in drug_rel:
            self.drug_update(each[0])
            self.drug_update(each[1])
        '''
        drug_targets, _, _ = cr.csv_read("./bio-data/bio-decagon-targets-all.csv", line=0)
        for each in drug_targets:
            self.drug_update(each[0])
        drug_targets, _, _ = cr.csv_read("./bio-data/bio-decagon-targets.csv", line=0)
        for each in drug_targets:
            self.drug_update(each[0])
        '''

    def load_gene_dict(self,fname = "./bio-data/bio-decagon-ppi.csv",gene_line = 0):
        print("loading gene data")
        if(gene_line == 0):
            gene_rel,_,_ = cr.csv_read(fname,line = 0 ,trans_int = True)
        else:
            gene_rel,_,_ = cr.csv_read(fname,line = gene_line ,trans_int = True)

        for each in gene_rel:
            self.gene_update(each[0])
            self.gene_update(each[1])
            self.gene_list.append(tuple(each))

        rel,l,c = cr.csv_read('./bio-data/bio-decagon-targets.csv')
        for each in rel:
            self.gene_update(int(each[1]))

        rel,l,c = cr.csv_read('./bio-data/bio-decagon-targets-all.csv')
        for each in rel:
            self.gene_update(int(each[1]))

    def gene_net_load(self):
        print("constructing gene graph")
        self.gene_graph.add_edges_from(self.gene_list)
        for each in self.gene2num_dict.keys():
            self.gene_graph.add_node(each)

    def drug_drug_adj_load(self, fname="./bio-data/bio-decagon-combo.csv", drug_line=0):
        print("constructing drugs adj-matrix from file:"+fname)
        if (drug_line == 0):
            drug_rel, l, c = cr.csv_read(fname, line=0)
        else:
            drug_rel, l, c = cr.csv_read(fname, line=drug_line)

        # Record drug pairs for each kind of side effect
        for each in drug_rel:
            if (each[2] not in self.drug_drug_dict.keys()):
                self.drug_drug_dict[each[2]] = [((each[0], each[1]), each[3])]
            else:
                self.drug_drug_dict[each[2]].append(((each[0], each[1]), each[3]))

        # Construct drug2num dict and form the adjacent matrix for drug-drug relations
        for each in set(list(self.drug_drug_dict.keys())):
            mat = np.zeros((self.n_drugs, self.n_drugs))
            for rel in self.drug_drug_dict[each]:
                mat[self.drug2num_dict[rel[0][1]], self.drug2num_dict[rel[0][0]]] = mat[
                    self.drug2num_dict[rel[0][0]], self.drug2num_dict[rel[0][1]]] = 1

            self.drug_drug_adj_list.append(sp.csr_matrix(mat))

    def gene_drug_adj_load(self,fname='./bio-data/bio-decagon-targets-all.csv', rel_line=0):
        print("constructing gene-drug adj-matrix from file:"+fname)
        self.gene_drug_adj_mat = np.zeros((self.n_genes, self.n_drugs))
        if (rel_line == 0):
            rel, l, c = cr.csv_read(fname, line=0)
        else:
            rel, l, c = cr.csv_read(fname, line=rel_line)
        for each in rel:
            if(each[0] not in self.drug2num_dict.keys()):
                self.lost_drug_count += 1
                continue
            self.gene_drug_adj_mat[self.gene2num_dict[int(each[1])], self.drug2num_dict[each[0]]] = 1
        self.gene_drug_adj_mat = sp.csr_matrix(self.gene_drug_adj_mat)

    def update_length(self):
        self.n_genes = len(self.gene2num_dict.keys())
        self.n_drugs = len(self.drug2num_dict.keys())
        print("Totally:\n "+str(self.n_drugs)+" drugs")
        print(" " + str(self.n_genes) + " genes")

    def save(self,path = "./"):
        print("Saving dics...")
        print("lost drugs: ",str(self.lost_drug_count))
        f0 = open(path+"drug2num.pkl",'wb')
        f1 = open(path + "num2drug.pkl", 'wb')
        f2 = open(path + "gene2num.pkl", 'wb')
        f3 = open(path + "num2gene.pkl", 'wb')
        pickle.dump(self.drug2num_dict,f0)
        pickle.dump(self.gene2num_dict,f2)
        pickle.dump(self.num2gene_dict,f3)
        pickle.dump(self.num2drug_dict,f1)
        f0.close()
        f1.close()
        f2.close()
        f3.close()
