import pandas
import sys
import config
import gene_identifiers as gid
import format_functions as ff


class Bio:
    """
    A class for generating BioGRID interaction annotations.

    Attributes
    ----------
    self.bio_df: pandas dataframe
        A dataframe containing the BioGRID interaction data.
    self.gid_structure: immuVCF object
        For converting genes between different identifications
    self.ncbi_iei_list: list
        List of known IEI genes in NCBI Entrez identification
    self.interaction_dict: dict
        For storing interactions with known IEI genes
    self.gene_list: list
        List of known IEI interactors
    self.interactor_list: list
        Stores all IEI interactors relevant to every gene in self.gene_list
    self.indicator_annotation: list
        1 codes an interaction, 0 codes no known interactions with IEI genes
    self.interaction_annotation: list
        Interaction annotation data for the VCF. Contains the IEI interactor
    """
    def __init__(self):
        print('\nLoading BioGRID...')
        with open(config.bio_filename) as bio_fh:
            self.bio_df = pandas.read_table(bio_fh, dtype='str')

        iei_df = pandas.read_excel(config.iei_gene_filename, sheet_name=0, header=None)
        iei_list = iei_df[0].tolist()

        # Create a list of genes that may need to be translated
        gid_gene_list = self.bio_df['Entrez Gene Interactor A'].tolist()

        # Initialize Gene ID and generated NCBI IEI list
        self.gid_structure = gid.GeneID(gid_gene_list, 'NCBI')
        self.ncbi_iei_list = self.gid_structure.id_conversion(iei_list, 'HGNC_symbol', 'NCBI')
        # print(self.ncbi_iei_list)

        # Create interaction dictionary
        self.interaction_dict = {}
        self.gene_list = []
        self.interactor_list = []

        # Annotation Attributes
        self.indicator_annotation = []
        self.interaction_annotation = []

        return

    def generate_interaction_list(self):
        """
        Generates a custom BioGIRD CSV file to store relevant interaction data for generating annotations.
        NOTE: This method should only be used when updating BioGRID data, as it takes a long time to process.
        """
        for ncbi in self.ncbi_iei_list:
            search_results = self.bio_df[self.bio_df['Entrez Gene Interactor A'] == ncbi]
            # print(search_results)
            interactors = search_results['Entrez Gene Interactor B'].tolist()
            # print(interactors)
            for gene in interactors:
                if gene not in self.interaction_dict.keys():
                    self.interaction_dict[gene] = [ncbi]
                    self.gene_list.append(gene)
                else:
                    self.interaction_dict[gene].append(ncbi)

        # Create interaction annotation lists
        for gene in self.gene_list:
            iei_interactors = self.interaction_dict[gene]
            self.interactor_list.append(self.gid_structure.id_conversion(iei_interactors, 'NCBI', 'HGNC_symbol'))
            # self.interactor_list.append(interactors)

        # Convert interaction lists to print strings
        self.interactor_list = ff.list_print_list(self.interactor_list)

        # Convert gene list to HGNC symbol
        self.gene_list = self.gid_structure.id_conversion(self.gene_list, 'NCBI', 'HGNC_symbol')

        # Create dataframe and Excel file
        annotation_dict = {'Gene': self.gene_list, 'IEI_Interactor': self.interactor_list}
        annotation_df = pandas.DataFrame(annotation_dict)
        if config.write_to != config.biogrid_df_name:
            print('Caution: In config.py: write_to is not the same file as biogrid_df_name. '
                  'Writing dataframe to biogrid_df_name')
        if config.biogrid_df_name.endswith('.csv'):
            annotation_df.to_csv(config.biogrid_df_name)
        else:
            print('Error! config.biogrid_df_name must end with .csv!')
            sys.exit()

        return self.gene_list, self.interactor_list

    def generate_annotations(self, gene_list, id_type):
        """
        Generates two lists to annotate the VCF with BioGRID data.

        Parameters
        ----------
        gene_list: list
            Genes to annotate
        id_type: str
            Type of gene ID. Must be compatible with the GeneID class.
        """
        bio_df = pandas.read_csv(config.biogrid_df_name, dtype='str')
        biogrid_genes = bio_df['Gene'].tolist()
        biogrid_interactions = bio_df['IEI_Interactor'].tolist()
        sym_list = self.gid_structure.id_conversion(gene_list, id_type, 'HGNC_symbol')
        for sym in sym_list:
            if sym in biogrid_genes:
                self.indicator_annotation.append(1)
                i = biogrid_genes.index(sym)
                self.interaction_annotation.append(biogrid_interactions[i])
            else:
                self.indicator_annotation.append(0)
                self.interaction_annotation.append('')
        return
