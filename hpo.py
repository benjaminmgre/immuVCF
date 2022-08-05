import pandas
import openpyxl
import config
import gene_identifiers as gid
from format_functions import print_list as pl


class Hpo:
    """
    Object for generating HPO annotations.

    Methods:
    --------
    annotate_genes:
        Takes a list of genes and generates HPO annotations (which can be accessed as object attributes)

    Attributes
    ----------
    id_list: list
        All HPO IDs
    pheno_to_gene_df: pandas dataframe
        Loads the HPO phenotype to gene annotation document into a dataframe
    null_phenotypes: list
        HP IDs for phenotypes with no relevant gene associations
    id_description_dict: dict
        Relates HPO IDs to HPO descriptions
    ncbi_to_ids_dict: dict
        Converts genes to their relevant phenotypes (Gene Name => HP ID)
    gen_description_annotation: list
        General HPO description annotations to be added to the VCF
    gen_id_annotation: list
        General HP ID annotations to be added to the VCF
    immune_description_annotation: list
        Immune HPO description annotations to be added to the VCF. Empty if config.separate_immune_hpo_terms is false
    immune_id_annotation: list
        Immune HPO ID annotations to be added to the VCF. Empty if config.separate_immune_hpo_terms is false
    count_annotation: list
        HPO Phenotype count annotations to be added to the VCF
    """

    def __init__(self, hpo_id_handle: str):
        """
        Initializes phenotype to gene dataframe, initializes hpo ID list,
        loads HPO descriptions from the ID list, generates gene to HP ID dictionary

        Parameters
        ----------
        hpo_id_handle: str
            Handle for file containing HPO IDs
            Format: ID values separated by newline characters
        """

        print('\nLoading HPO...')
        # INITIALIZE PHENO TO GENE DATAFRAME
        labels = ['HPO-id', 'HPO label', 'entrez-gene-id', 'entrez-gene-symbol', 'Additional Info from G-D source',
                  'G-D source', 'disease-ID for link']
        try:
            with open(config.ph2g_filename) as pheno_to_gene:
                self.pheno_to_gene_df = pandas.read_table(pheno_to_gene, dtype='str', comment='#',
                                                          names=labels)
        except IOError:
            print('Error! HPO Phenotype to Disease file could not be found.')

        # INITIALIZE HPO ID LIST
        if config.input_hpo == 'List':
            try:
                id_file = open(hpo_id_handle)
            except IOError:
                print('Error! HPO ID file could not be found.')
                return
            self.id_list = [line.rstrip() for line in id_file]
            id_file.close()
        else:
            # Phenotips
            pt_df = pandas.read_table(config.hpo_id_filename, sep='\t', comment='#', dtype='str')
            # Select HPO IDs column
            hpo_col = pt_df['HPO IDs'].tolist()
            self.id_list = []
            for group in hpo_col:
                ids = group.split('; ')
                for hp_id in ids:
                    if hp_id in self.id_list:
                        continue
                    else:
                        self.id_list.append(hp_id)
            print(f'{len(self.id_list)} HPO IDs found. List: {self.id_list}')

        # LOAD HPO DESCRIPTIONS
        if len(self.id_list) == 0:
            print('Error! No initiated HPO IDs found')
            return
        self.null_phenotypes = []
        self.id_description_dict = {}
        for hpo_id in self.id_list:
            hpo_phenotypes = self.pheno_to_gene_df[self.pheno_to_gene_df['HPO-id'] == hpo_id]
            # If there are no genes associated with a phenotype, add HP ID to null_phenotypes list
            if hpo_phenotypes.empty:
                self.null_phenotypes.append(hpo_id)
                continue
            else:
                hpo_label = hpo_phenotypes.iloc[0]['HPO label']
                self.id_description_dict[hpo_id] = hpo_label

        # GENERATE NCBI ID => HP ID DICTIONARY
        self.ncbi_to_ids_dict = {}
        for hpo_id in self.id_list:
            related_genes = self.pheno_to_gene_df[self.pheno_to_gene_df['HPO-id'] == hpo_id]['entrez-gene-id'].tolist()
            related_genes = list(set(related_genes))
            for ncbi_gene in related_genes:
                if ncbi_gene not in self.ncbi_to_ids_dict:
                    self.ncbi_to_ids_dict[ncbi_gene] = [hpo_id]
                else:
                    self.ncbi_to_ids_dict[ncbi_gene].append(hpo_id)

        # Generate Immune HPO Term List
        immune_hpo_df = pandas.read_excel(io=config.immune_HPO_filename, sheet_name=0)
        self.immune_HPO_list = immune_hpo_df['HPO ID'].tolist()

        # Initialize Annotation Attributes
        self.count_annotation = []

        self.gen_description_annotation = []
        self.gen_id_annotation = []

        self.immune_description_annotation = []
        self.immune_id_annotation = []

        print('HPO Initialization Complete')
        return

    def annotate_genes(self, gene_list: list, id_type: str):
        """
        Generates HPO annotations for a list of genes (in any available format).

        Parameters
        ----------
        gene_list: list of strings
            List of genes to annotate
        id_type: str
            Gene ID type of the ID list
            Available types: 'NCBI': Entrez ID, 'Ensembl': Ensembl ID, 'RefSeq': RefSeq ID,
            'HGNC_symbol': Approved Gene Symbol, 'HGNC_id': HGNC ID

        Returns
        -------
        self.count_annotation: list
            Respective list with an integer representing the number of phenotypes relevant to the gene
        self.id_annotation: list
            Respective list of lists containing all relevant HPO phenotype identification numbers for each gene
        self.description_annotation: list
            Respective list of lists containing all relevant HPO phenotype labels for each gene
        """
        separate_immune = config.separate_immune_hpo_terms

        # Initialize Gene Conversion Structure
        gid_structure = gid.GeneID(gene_list, id_type)
        ncbi_list = gid_structure.id_conversion(gene_list, id_type, 'NCBI')

        # Generate HPO Annotations for genes.
        # Description: HPO Label, ID: HPO ID, Count: Number of related phenotypes per gene
        for gene in ncbi_list:
            if gene in self.ncbi_to_ids_dict.keys():
                hpo_ids = self.ncbi_to_ids_dict[gene]
                num_hpo = len(self.ncbi_to_ids_dict[gene])
                gen_pheno_descs = []
                immune_pheno_descs = []
                gen_hpo_ids = []
                immune_hpo_ids = []
                for hpo_id in hpo_ids:
                    pheno_desc = self.id_description_dict[hpo_id]
                    # If immune HPO, add to the immune HPO annotation
                    if hpo_id in self.immune_HPO_list and separate_immune:
                        immune_hpo_ids.append(hpo_id)
                        immune_pheno_descs.append(pheno_desc)
                    # Else add to the general annotation
                    else:
                        gen_hpo_ids.append(hpo_id)
                        gen_pheno_descs.append(pheno_desc)
            else:
                gen_hpo_ids = []
                gen_pheno_descs = []
                immune_hpo_ids = []
                immune_pheno_descs = []
                num_hpo = ''
            self.count_annotation.append(num_hpo)
            self.gen_id_annotation.append(gen_hpo_ids)
            self.gen_description_annotation.append(gen_pheno_descs)
            self.immune_id_annotation.append(immune_hpo_ids)
            self.immune_description_annotation.append(immune_pheno_descs)
        print('HPO Annotation Complete')
        return

    def hpo_genes_list(self, file_name):
        """
        Generates an annotated gene list in CSV format.

        Parameters
        ----------
        file_name: str
            What to name the generated file
        """
        annotation_list = []
        unique_gene_list = []
        ncbi_list = []
        # Initiate
        # Iterate over every given HP Phenotype ID
        for hpo in self.id_list:
            hpo_genes = self.pheno_to_gene_df[self.pheno_to_gene_df['HPO-id'] == hpo]
            # Iterate over every gene associated with the HP Phenotype ID
            for gene in hpo_genes.iterrows():
                hpo_id = gene[1].loc["HPO-id"]
                hpo_label = gene[1].loc["HPO label"]
                gene_name = gene[1].loc["entrez-gene-symbol"]
                ncbi_id = gene[1].loc["entrez-gene-id"]
                ncbi_list.append(ncbi_id)
                # If the gene is not in the annotation list, add it as a dictionary
                if gene_name not in unique_gene_list:
                    gene_entry = {"Gene-Name": gene_name, "HPO-Count": 1, "HPO-IDs": [hpo_id],
                                  "HPO-Labels": [hpo_label], "Entrez-ID": ncbi_id, "Ensembl-ID": None}
                    annotation_list.append(gene_entry)
                    # Unique gene list makes sure repeated genes are accounted for in future steps
                    unique_gene_list.append(gene_name)
                else:
                    # Determine the index of the repeated gene (Should have the same index as annotation_list
                    i = unique_gene_list.index(gene_name)
                    # If the specific HPO term is already related to the gene, skip. (HPO file has phenotype duplicates)
                    if hpo_id in annotation_list[i]["HPO-IDs"]:
                        continue
                    else:
                        # Make HPO ID, Count, and Label additions
                        annotation_list[i]["HPO-Count"] += 1
                        annotation_list[i]["HPO-IDs"].append(hpo_id)
                        annotation_list[i]["HPO-Labels"].append(hpo_label)
        # Generate Ensembl Gene ID
        gid_structure = gid.GeneID(ncbi_list, 'NCBI')
        for entry in annotation_list:
            ncbi_id = entry["Entrez-ID"]
            ens_id = gid_structure.single_conversion(ncbi_id, "NCBI", "Ensembl")
            entry["Ensembl-ID"] = ens_id
            # For each entry, turn lists into print strings
            entry["HPO-IDs"] = pl(entry["HPO-IDs"])
            entry["HPO-Labels"] = pl(entry["HPO-Labels"])

        # Convert to a pandas dataframe, sort, and export to CSV file format
        df = pandas.DataFrame(annotation_list)
        df.sort_values('HPO-Count', ascending=0)
        df.to_excel(file_name)
        return


class IeiHpo:
    """
    A class for annotating known IEI lists.

    Methods
    -------
    generate_annotations: Generates a CSV of IEI genes and respective HPO phenotype annotations
    annotate_HPO_terms: Generates a CSV of immune HPO term numbers and labels
    """

    def __init__(self):
        # Initialize the gene to phenotype dataframe
        labels = ['entrez-gene-id', 'entrez-gene-symbol', 'HPO-Term-ID', 'HPO-Term-Name', 'Frequency-Raw',
                  'Frequency-HPO', 'Additional Info from G-D source', 'G-D source', 'disease-ID for link']
        with open(config.g2ph_filename) as g2ph:
            self.g2ph_df = pandas.read_table(g2ph, dtype='str', comment='#', names=labels)

        # Generate the IEI lists
        iei_df = pandas.read_excel(config.iei_gene_filename, sheet_name=0, header=None)
        blueprint_df = pandas.read_excel(config.iei_gene_filename, sheet_name=1, header=None)
        veoibd_df = pandas.read_excel(config.iei_gene_filename, sheet_name=2, header=None)

        self.iei_list = iei_df[0].tolist()
        self.blueprint_list = blueprint_df[0].tolist()
        self.veoibd_list = veoibd_df[0].tolist()

        print(f'IEI List: {self.iei_list}')
        print(f'Blueprint list: {self.blueprint_list}')
        print(f'VEO IBD List: {self.veoibd_list}')

        # Initialize the Immune HPO terms list
        immune_hpo_df = pandas.read_excel(io=config.immune_HPO_filename, sheet_name=0)
        self.immune_HPO_list = immune_hpo_df['HPO ID'].tolist()

        # Initialize gene id conversion structure
        self.gid_structure = gid.GeneID(self.iei_list, 'HGNC_symbol')
        # Generate NCBI IEI List
        self.ncbi_iei_list = self.gid_structure.id_conversion(self.iei_list, 'HGNC_symbol', 'NCBI')
        # Initialize attributes:
        self.immu_hpo_id_annotations = []
        self.immu_hpo_desc_annotations = []
        self.gen_hpo_id_annotations = []
        self.gen_hpo_desc_annotations = []

        self.blueprint_annotations = []
        self.veoibd_annotations = []
        self.cdg_annotations = []

        self.immune_HPO_label_annotation = []
        return

    def generate_annotations(self, filename):
        """
        Generates a CSV file containing known IEI genes and their respective HPO phenotypes (separated by immunity)

        Parameters
        ----------
        filename: str
            Where to write the CSV file
        """
        for iei in self.iei_list:
            ncbi = self.gid_structure.single_conversion(iei, 'HGNC_symbol', 'NCBI')
            search_results = self.g2ph_df[self.g2ph_df['entrez-gene-id'] == ncbi]
            phenotype_ids = search_results['HPO-Term-ID'].tolist()
            phenotype_labels = search_results['HPO-Term-Name'].tolist()
            immune_ids = []
            immune_labels = []
            gen_ids = []
            gen_labels = []
            for i in range(0, len(phenotype_ids)):
                hpo_id = phenotype_ids[i]
                label = phenotype_labels[i]
                # Skip all HPO terms that have been annotated already
                if hpo_id in immune_ids or hpo_id in gen_ids:
                    continue
                if hpo_id in self.immune_HPO_list:
                    # Generate Immunity Annotations
                    immune_ids.append(hpo_id)
                    immune_labels.append(label)
                else:
                    # Generate General Annotations
                    gen_ids.append(hpo_id)
                    gen_labels.append(label)
            # pl() formats all lists into printable strings
            self.immu_hpo_id_annotations.append(pl(immune_ids))
            self.immu_hpo_desc_annotations.append(pl(immune_labels))
            self.gen_hpo_id_annotations.append(pl(gen_ids))
            self.gen_hpo_desc_annotations.append(pl(gen_labels))

            # Generate Panel Annotations
            if iei in self.blueprint_list:
                self.blueprint_annotations.append(1)
            else:
                self.blueprint_annotations.append(0)
            if iei in self.veoibd_list:
                self.veoibd_annotations.append(1)
            else:
                self.veoibd_annotations.append(0)

        # Open CDG Database and generate annotations
        with open(config.cdg_filename) as cdg:
            cdg_df = pandas.read_table(cdg)
        for iei in self.iei_list:
            cdg_interactions = cdg_df[cdg_df['Gene'] == iei]['HGMD_Gene'].tolist()
            if 'is_known_HGMD_gene' in cdg_interactions:
                cdg_interactions.remove('is_known_HGMD_gene')
            if 'No_such_gene_in_HGC_database' in cdg_interactions:
                cdg_interactions.remove('No_such_gene_in_HGC_database')
            self.cdg_annotations.append(pl(cdg_interactions))

        # Generate Annotated CSV
        frame = {'gene': self.iei_list, 'Immune-HPO-IDs': self.immu_hpo_id_annotations,
                 'Immune-HPO-Names': self.immu_hpo_desc_annotations, 'General-HPO-IDs': self.gen_hpo_id_annotations,
                 'General-HPO-Names': self.gen_hpo_desc_annotations, 'Blueprint Panel': self.blueprint_annotations,
                 'VEO IBD Panel': self.veoibd_annotations, 'CDG Interactions': self.cdg_annotations}
        export_df = pandas.DataFrame(frame)
        export_df.to_excel(filename)
        return

    def annotate_hpo_terms(self, filename):
        """
        Creates a CSV file storing the contents of self.immune_HPO_list and every related phenotype label

        Parameters
        ----------
        filename: str
            Where to write the annotation file
        """
        for hpo in self.immune_HPO_list:
            search_results = self.g2ph_df[self.g2ph_df['HPO-Term-ID'] == hpo]
            if search_results.empty:
                self.immune_HPO_label_annotation.append('')
                continue
            else:
                label = search_results['HPO-Term-Name'].iloc[0]
                self.immune_HPO_label_annotation.append(label)
        annotation_dict = {'HPO ID': self.immune_HPO_list, 'HPO Name': self.immune_HPO_label_annotation}
        annot_hpo_df = pandas.DataFrame(annotation_dict)
        annot_hpo_df.to_csv(filename)
        return
