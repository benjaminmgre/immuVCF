import pandas
import config
import gene_identifiers as gid


class Cdg:
    """
    Class for generating the Closest Disease Causing Gene (CDG) list and annotations.

    Methods:
    --------
    __init__:
        Take the file handle for the CDG list and generates a CDG gene list in python.
    annotate_cdg:
        Takes a list of genes and returns a list of the same size with 1 for gene in the CDG list and 0 for not.

    Attributes:
    -----------
    cdg_df: pandas dataframe
        Loads the CDG file into python for data analysis
    cdg_list: list
        All CDG Genes (Does not contain genes on the known IEI list unless it is also related with another known
        gene through CDG.)
    cdg_dict: dictionary
        Stores CDG genes as keys, and every respective gene's known IEI interactors as list for values.
    """

    def __init__(self):
        print('\nLoading CDG...')
        # Generate CDG list from file
        with open(config.cdg_filename) as cdg:
            self.cdg_df = pandas.read_table(cdg)
        self.cdg_list = []
        self.cdg_dict = {}
        for i, row in self.cdg_df.iterrows():
            hgmd_gene = row['HGMD_Gene']
            interactor = row['Gene']
            if hgmd_gene == 'is_known_HGMD_Gene':
                continue
            elif hgmd_gene in self.cdg_list:
                self.cdg_dict[hgmd_gene].append(interactor)
            else:
                self.cdg_dict[hgmd_gene] = [interactor]
                self.cdg_list.append(hgmd_gene)
        print(str(len(self.cdg_list)) + ' CDG Genes Generated')

        # Initialize Annotation Attribute
        self.indicator_annotation = []
        self.interaction_annotation = []
        print('CDG Initialization Complete')
        return

    def annotate_cdg(self, gene_list, id_type):
        """
        Generates an annotation list attribute that indicates whether a given gene is on the Closest Disease Gene (CDG)
        list.

        Parameters:
        -----------
        gene_list: list
            Genes to be annotated by the CDG module.
        id_type: str
            Gene ID type of the ID list
            Available types: 'NCBI': Entrez ID, 'Ensembl': Ensembl ID, 'RefSeq': RefSeq ID,
            'HGNC_symbol': Approved Gene Symbol, 'HGNC_id': HGNC ID

        Returns:
        --------
        self.indicator_annotation: list
            Respective list of True (1) or False (0) that indicates whether the given gene is on the CDG list.
        self.interaction_annotation: list
            Respective list of known IEI genes that interact with the annotated gene on the CDG list.
        """
        # Initialize Gene Conversion Structure
        gid_structure = gid.GeneID(gene_list, id_type)
        hgnc_list = gid_structure.id_conversion(gene_list, id_type, 'HGNC_symbol')
        # DEBUG:  print(f'Gene List: {gene_list}, hgnc List: {hgnc_list}, CDG List: {self.cdg_list}')
        for gene in hgnc_list:
            if gene in self.cdg_list:
                self.indicator_annotation.append(1)
                self.interaction_annotation.append(self.cdg_dict[gene])
            else:
                self.indicator_annotation.append(0)
                self.interaction_annotation.append('')
        print('CDG Annotation Complete')
        return
