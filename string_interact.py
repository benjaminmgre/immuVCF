import pandas
import sys
import config
import gene_identifiers as gi
import format_functions as ff


def format_search(i):
    """
    Formats a protein ID to fit the STRING interaction dataframe. 9606 is the identifier for human.
    """
    format_id = '9606.' + i
    return format_id


def clear_format(i):
    """
    Removes the STRING formatting.
    """
    clear_id = i[5:]
    return clear_id


class String:
    """
    A class for generating STRING interaction annotations.

    Attributes
    ----------
    self.s_df: pandas dataframe
        STRING interaction information dataframe
    self.iei_list: list
        List of known IEI genes in HGNC ID form
    self.conversion_structure: immuVCF GeneID object
        Used to convert genes between different ID forms
    self.IEI_ens: list
        Known IEI genes in Ensembl Gene ID form
    self.threshold: int
        The STRING confidence score threshold value. All interactions with confidence value under the threshold are
        removed.
    self.indicator_list: list
        A list of 1 (True) or 0 (False) to indicate whether a gene has a known IEI interaction above the set threshold
    self.interaction_list: list
        Contains all known IEI interactions. Every interaction is recorded with a tuple: (IEI_gene, STRING_score)
    self.iei_indicator_annotation: list
        Indicates if a given gene is a known IEI gene
    """
    def __init__(self):
        print('\nLoading STRING...')
        # Initialize the STRING interaction dataframe
        with open(config.string_filename) as string_fh:
            self.s_df = pandas.read_table(string_fh, delimiter=' ')

        # Initialize the known IEI list
        iei_df = pandas.read_excel(config.iei_gene_filename, sheet_name=0, header=None)
        self.iei_list = iei_df[0].tolist()

        # Convert IEI HGNC symbol list to IEI ensembl list
        self.conversion_structure = gi.GeneID(self.iei_list, 'HGNC_symbol')
        self.IEI_ens = self.conversion_structure.id_conversion(self.iei_list, 'HGNC_symbol', 'Ensembl')

        # Set the STRING search threshold
        self.threshold = config.string_threshold

        # Initialize annotation attributes
        self.indicator_list = []
        self.interaction_list = []
        self.iei_indicator_annotation = []

        print('STRING Initialization Complete')
        return

    def count_interactions(self):
        """
        Generates a STRING interaction df that is saved in a csv file.
        Contains all interactions with known IEI genes.
        This function only needs to be run to create or update the dataframe. It does not have to be run for every VCF.

        Returns
        -------
        interaction_df: csv file
            Contains all known and predicted STRING interactions with known IEI genes above the self.threshold value
        """
        gp_structure = gi.ProteinID()
        interaction_list = []
        index_list = []
        # TODO: See what conversion errors are problems and fix
        for gene_ens in self.IEI_ens:
            # Check what position the gene is in.
            i = self.IEI_ens.index(gene_ens)
            l = len(self.IEI_ens)
            if i % 10 == 0:
                print(f'{i} out of {l} complete.')
            # Convert Gene Ensembl ID to its possible Ensembl Protein IDs
            prot_ens = gp_structure.gene_to_proteins(gene_ens)
            if prot_ens is None:
                # print(f'{gene_ens} has no protein links.')
                continue
            for prot in prot_ens:
                # Format the search term to match the STRING file format
                fprot = format_search(prot)
                search_results = self.s_df[self.s_df['protein1'] == fprot]
                if search_results.empty:
                    continue
                else:
                    for _, row in search_results.iterrows():
                        interactor_protein = clear_format(row['protein2'])
                        interactor_gene = gp_structure.protein_to_gene(interactor_protein)
                        score = row['combined_score']
                        gene_hgnc = self.conversion_structure.single_conversion(gene_ens, 'Ensembl', 'HGNC_symbol')
                        interaction = f"{gene_hgnc}:{score}"
                        if interactor_gene not in index_list:
                            interaction_dict = {'gene': interactor_gene, 'protein': interactor_protein,
                                                'interactions': [interaction]}
                            interaction_list.append(interaction_dict)
                            index_list.append(interactor_gene)
                        else:
                            i = index_list.index(interactor_gene)
                            interaction_list[i]['interactions'].append(interaction)
        # Turn interactions into printable strings
        for dct in interaction_list:
            dct['interactions'] = ff.print_list(dct['interactions'])
        interaction_df = pandas.DataFrame(interaction_list)
        if config.write_to != config.string_df_name:
            print('Caution: In config.py: write_to does not equal string_df_name. Writing to string_df_name')
        if config.string_df_name.endswith('.csv'):
            interaction_df.to_csv(config.string_df_name)
        else:
            print('Error! config.string_df_name must end in .csv!')
            sys.exit()
        print('Interaction counting complete!')
        return

    def generate_annotations(self, ens_list):
        """
        Generate STRING interaction annotation list.
        NOTE: This function requires the generation of an interaction dataframe, which is generated
        by the count_interactions method.

        Parameters:
        -----------
        ens_list: list
            List of genes in Ensembl ID format to be annotated

        Returns:
        --------
        indicator_list: list
            Annotation list for VCF. 1 represents an interaction was found, 0 means no interaction was found within
            the threshold value of the string dataframe.
        interaction_list: list
            Annotation list for VCF. Contains couples of interactor genes and STRING interaction score.
        """
        with open(config.string_df_name) as sl_fh:
            # Initialize the string interaction list dataframe.
            sl_df = pandas.read_csv(sl_fh)
        # Generate a list of all possible interactor genes.
        sl_ens = sl_df['gene'].tolist()
        for ens in ens_list:
            if ens in sl_ens:
                self.indicator_list.append(1)
                # Query for the interaction entry, search the interaction columns,
                # get the first row (there should only be 1 row).
                interaction_string = sl_df[sl_df['gene'] == ens]['interactions'].iloc[0]
                interactions = interaction_string.split(", ")
                interaction_annotation = []
                for interaction in interactions:
                    # Interactions are in format "interactor:score"
                    group = interaction.split(":")
                    interactor = group[0]
                    score = group[1]
                    if int(score) > self.threshold:
                        inter_annotate = f"({interactor}, {score})"
                        interaction_annotation.append(inter_annotate)
                    else:
                        # Does not meet score threshold
                        continue
                self.interaction_list.append(ff.print_list(interaction_annotation))
            else:
                self.indicator_list.append(0)
                self.interaction_list.append('')
        print('STRING Annotation Complete')
        return self.indicator_list, self.interaction_list

    def annotate_iei(self, gene_list, id_type='Ensembl'):
        """
        Generates an annotation list attribute that indicates whether a given gene is a known IEI disease gene.

        Parameters
        ----------
        gene_list: list
            Genes to be annotated for IEI genes
        id_type: str
            The type of gene id the gene list is formatted in. (See vcf.py)

        Returns
        -------
        iei_indicator_annotation: list
            List of indicators of whether a gene is a known IEI disease gene. 1 is True, 0 is False
        """
        gid_structure = gi.GeneID(gene_list, id_type)
        gene_list = gid_structure.id_conversion(gene_list, id_type, 'Ensembl')
        for gene in gene_list:
            if gene in self.IEI_ens:
                self.iei_indicator_annotation.append(1)
            else:
                self.iei_indicator_annotation.append(0)
        return self.iei_indicator_annotation
