import pandas
import config
import gene_identifiers as gid


class FindGene:
    """
    A class for creating gene ID annotations from chromosome and base pair coordinate numbers. Annotation provides
    information to duplicate coordinates which may have multiple protein-coding genes.

    Attributes:
    -----------
    self.ens_annotation: list
        A list of Ensembl gene IDs the size of the VCF + number of duplicates that contains the Ensembl Gene IDs.
        Note: duplicates are placed underneath the original entry.
    self.duplicate_indicator: list
        A list of integers the size of the chrom_cord_list that indicates how many duplicates per row should be
        generated.
    self.ec_df: pandas dataframe
        Ensembl coordinate dataframe; used to find all genes related to a specific chromosome and coordinate.
    """
    def __init__(self, chrom_coord_list):
        """
        Parameters
        ----------
        chrom_coord_list: list
            Chromosome and coordinate values. In format: "chrom#:coord#". Taken from the VCF.
        """
        # Initialize Gene Coordinates Dataframe
        if config.reference_genome == 'GRCh38':
            ensembl_coords_name = config.gencode_GRCh38_filename
        else:
            ensembl_coords_name = config.gencode_GRCh37_filename
        try:
            with open(ensembl_coords_name) as ec:
                # Define column names
                column_labels = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'genomic phase',
                                 'additional information']
                # Read the GENCODE GFT File
                self.ec_df = pandas.read_table(ec, sep="\t", names=column_labels, comment='#', dtype='str')
                # Set base pair coordinates to int data type for future comparisons
                self.ec_df = self.ec_df.astype({'start': 'int', 'end': 'int'})
        except IOError:
            print('Error! Ensembl coordinates file could not be found.')
        self.ens_annotation = []
        self.duplicate_indicator = []
        removed_genes = 0
        for entry in chrom_coord_list:
            vals = entry.split(':')
            chrom = 'chr' + str(vals[0])
            coord = int(vals[1])

            """
            Feature Search: Removes all non-gene entries (ex. RNA or exons or proteins)
            Chrom Search: Removes all entries not in the correct chromosome
            Chrom Search Start: Removes all genes that start after the coordinate
            Chrom Search: Removes all genes that end before the coordinate
            Additional Info: Take the additional info column of every gene, and returns a list
            """
            feature_search = self.ec_df[self.ec_df['feature'] == 'gene']
            chrom_search = feature_search[feature_search['chrom'] == chrom]
            coord_search_start = chrom_search[chrom_search['start'] <= coord]
            coord_search = coord_search_start[coord_search_start['end'] >= coord]
            add_info = coord_search['additional information'].tolist()
            # If no genes were found report no gene found
            if len(add_info) == 0:
                gene = 'No gene found'
                self.ens_annotation.append(gene)
                dups = 0
            # If one gene found, add to Ensembl annotation
            elif len(add_info) == 1:
                gene = add_info[0].split(';')[0].split('"')[1].split('.')[0]
                self.ens_annotation.append(gene)
                dups = 0
                # gene_type = add_info[0].split('; ')[1].split('"')[1]
            # If multiple genes are found, remove all non-protein-coding genes if filter_npc is True, then present
            # every remaining gene. Entries with multiple valid genes should have duplicate rows.
            else:
                # DEBUG: print(f'Multiple genes found for entry: {entry}. Genes: {add_info}.')
                genes = []
                for row in add_info:
                    g = row.split(';')[0].split('"')[1].split('.')[0]
                    g_type = row.split(';')[1].split('"')[1]
                    # Do not append g to genes if filter_npc is True and g is not protein-coding
                    if g_type != 'protein_coding' and config.filter_non_protein_genes is True:
                        removed_genes += 1
                        continue
                    else:
                        genes.append(g)
                # DEBUG:
                if len(genes) == 0:
                    self.ens_annotation.append('No protein-coding genes')
                    dups = 0
                else:
                    dups = len(genes) - 1
                # Create annotation for every valid redundant gene
                for gene in genes:
                    self.ens_annotation.append(gene)
            self.duplicate_indicator.append(dups)

            # Print progress
            i = chrom_coord_list.index(entry)
            if i % 1000 == 0:
                print(f'{i} / {len(chrom_coord_list)} variants complete.')

        # Convert Ensembl Gene IDs to Approved Symbol
        gid_structure = gid.GeneID(self.ens_annotation, 'Ensembl')
        self.hgnc_annotation = gid_structure.id_conversion(self.ens_annotation, 'Ensembl', 'HGNC_symbol')
        return
