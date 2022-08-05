# PARAMETERS - DESCRIPTIONS AVAILABLE IN README.TXT

mode = 'annotate_VCF'
# Available Modes: 'annotate_VCF', 'generate_HPO_list', 'generate_STRING_df', 'generate_BioGRID_df', 'annotate_iei_list'

# Patient Dependent Files
vcf_filename = 'input_files/'
hpo_id_filename = 'input_files/'

write_to = 'generated_files/'

# Toggle Parameters
generate_genes = True
input_hpo = 'Phenotips'
# Available Input HPO: 'Phenotips' or 'List'
separate_immune_hpo_terms = True
filter_non_protein_genes = True
string_threshold = 700
reference_genome = 'GRCh37'
# Available Input: 'GRCh37' or 'GRCh38'

# Ensembl Coordinate File
gencode_GRCh38_filename = 'file_dependencies/GENCODE_GRCh38.gtf'
gencode_GRCh37_filename = 'file_dependencies/GENCODE_GRCh37.gtf'

# Gene ID Files
hgnc_filename = 'file_dependencies/hgnc_gene_id.txt'
gp_filename = 'file_dependencies/gene_proteins.txt'

# HPO Files
ph2g_filename = 'file_dependencies/phenotype_to_genes.txt'
g2ph_filename = 'file_dependencies/genes_to_phenotype.txt'
immune_HPO_filename = 'file_dependencies/immune_HPO_terms_July13.xlsx'

# Interaction Analysis Files
iei_gene_filename = 'file_dependencies/IEI_Genes_July26.xlsx'
cdg_filename = 'file_dependencies/CDG_data.txt'
bio_filename = 'file_dependencies/BIOGRID_data.txt'
string_filename = 'file_dependencies/STRING_data.txt'

# Generated Dataframe
string_df_name = 'generated_files/string_df_Aug3.csv'
biogrid_df_name = 'generated_files/biogrid_df_Aug3.csv'

# Column Names
ensembl_column = 'Ensembl_gene_id'
coord_column = 'Position'

# Generated Column Names
duplicate_column_label = 'Variant Duplicates'
gene_column_label = 'immuVCF_Gene'
