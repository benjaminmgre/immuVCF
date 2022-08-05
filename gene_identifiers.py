import pandas
import config


def id_tab(id_format):
    """
    Converts the given ID format string into the string which can be used to index the gene
    identity dataframes.

    Parameters:
    ----------
    id_format: str
        A string representing the identification forms of the initial gene list.

    Returns:
    --------
    id_t: str
        A string that can be used to index the gene identification dataframes. Matches id_format
    """
    if id_format == 'HGNC_symbol':
        id_t = 'Approved symbol'
    elif id_format == 'HGNC_id':
        id_t = 'HGNC ID'
    elif id_format == 'RefSeq':
        id_t = 'RefSeq IDs'
    elif id_format == 'NCBI':
        id_t = 'NCBI Gene ID'
    elif id_format == 'Ensembl':
        id_t = 'Ensembl gene ID'
    else:
        id_t = ''
        print("""
            Error! id_type was not entered correctly. Valid str: "HGNC_symbol", "HGNC_id", "RefSeq", "NCBI",
            "Ensembl".
              """)
    return id_t


class GeneID:
    """
    A dataframe for referencing various forms of a gene identifications.
    Supports HGNC IDs, HGNC Approved Symbols, RefSeq IDs, NCBI Gene IDs, Ensembl gene IDs

    Attributes:
    -----------
    self.id_format: str
        A dataframe column header which represents the given identification form from the gene list.
    self.gid_df: pandas dataframe
        Gene Identifier Dataframe: Generated from the hgnc identification file
    self.gene_log: pandas dataframe
        A filtered self.gid_df: Only contains genes given in the gene list
    """

    def __init__(self, gene_list, id_form):
        """
        Initiates a gene identity dataframe for all genes in the gene list.

        Parameters:
        -----------
        gene_list: list
            A list of genes in the specified id format
        id_form: str
            String representing the ID form of the gene list.
            'HGNC_symbol': Approved HGNC symbols
            'HGNC_id': HGNC ID
            'RefSeq': RefSeq ID
            'NCBI': NCBI Gene ID
            'Ensembl': Ensembl gene ID
        """
        self.id_format = id_tab(id_form)
        with open(config.hgnc_filename) as gid:
            self.gid_df = pandas.read_table(gid, dtype='str')
        gene_log_list = []
        self.gene_searched_list = []
        for gene in gene_list:
            # DEBUG: print(gene)
            if gene in self.gene_searched_list:
                continue
            else:
                gene_search = self.gid_df[self.gid_df[self.id_format] == gene]
                # print(gene_search)
                gene_log_list.append(gene_search)
                self.gene_searched_list.append(gene)
        self.gene_log = pandas.concat(gene_log_list)
        return

    def id_conversion(self, gene_list, from_id, to_id):
        """
        Converts a list of genes from one identification form into another.

        Parameters:
        -----------
        gene_list: list
            Genes to be returned in a different form
        from_id: str
            The ID format the gene list is in
        to_id: str
            The ID format the gene list will be converted to

        See __init__ for ID forms.

        Returns:
        --------
        converted_list: list
            Genes listed in the to_id format.
        """
        from_form = id_tab(from_id)
        to_form = id_tab(to_id)
        converted_list = []
        for gene in gene_list:
            gene_search = self.gene_log[self.gene_log[from_form] == gene]
            if gene_search.empty:
                # print(f'{gene} can not be converted to {to_id}.')
                # print(gene_search)
                converted_form = gene
            elif len(gene_search) > 1:
                # print(gene_search)
                converted_form = gene_search.iloc[0][to_form]
            else:
                converted_form = gene_search.iloc[0][to_form]
            converted_list.append(converted_form)
        return converted_list

    def single_conversion(self, gene, from_id, to_id):
        """
        Duplicate of id_conversion method, but only takes a single gene and returns a single gene.
        """
        from_form = id_tab(from_id)
        to_form = id_tab(to_id)
        gene_search = self.gene_log[self.gene_log[from_form] == gene]
        if gene_search.empty:
            converted_form = 'GENE ID CONVERSION ERROR!'
        else:
            converted_form = gene_search.iloc[0][to_form]
        return converted_form


class ProteinID:
    """
    A class for finding protein IDs encoded by gene IDs and vise-versa.

    Attributes
    ----------
    self.gp_df: pandas dataframe
        Gene protein dataframe. Contains a gene column, and a protein column. The protein is what the gene encodes.

    Methods
    -------
    self.gene_to_proteins:
        Takes a gene in Ensembl ID, returns a list of proteins it encodes.
    self.protein_to_gene:
        Takes a protein, and returns the gene that encodes it.
    """
    def __init__(self):
        with open(config.gp_filename) as gp:
            self.gp_df = pandas.read_table(gp)
        # Remove all gene columns without protein entries
        self.gp_df = self.gp_df.dropna()
        return

    def gene_to_proteins(self, gene):
        """
        Takes a gene in Ensembl ID format, returns a list of proteins in Ensembl ID format.

        Parameters
        ----------
        gene: str
            Gene in Ensembl ID format

        Returns
        -------
        related_prots: list
            List of proteins in Ensembl ID format that are encoded by the gene.
        """
        search_results = self.gp_df[self.gp_df['Gene stable ID'] == gene]
        if search_results.empty:
            return None
        else:
            related_prots = []
        for i, row in search_results.iterrows():
            prot = row['Protein stable ID']
            related_prots.append(prot)
        return related_prots

    def protein_to_gene(self, prot):
        """
        Takes a protein in Ensembl ID format, returns the gene that encodes the protein in Ensembl ID format.
        """
        # Search should only return one gene.
        search = self.gp_df[self.gp_df['Protein stable ID'] == prot]
        if search.empty:
            return 'missing gene'
        else:
            gene = search.iloc[0]['Gene stable ID']
        return gene
