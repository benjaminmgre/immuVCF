import config

import CDG
import hpo
import vcf
import biogrid as bp
import string_interact as si
import find_gene as fg
from format_functions import list_print_list as lpl


def annotate_vcf(vcf_file, hpo_file):
    print('\n\nWelcome to immuVCF!')
    print('Version: 1.0')
    print('Date Last Updated: August 2nd, 2022')

    # Initialize VCF Object, Call the ENSEMBL ID column
    vcf_structure = vcf.Vcf(vcf_file)
    if config.generate_genes is True:
        chrom_coords = vcf_structure.load_column(config.coord_column)
        fg_structure = fg.FindGene(chrom_coords)

        duplicate_column = config.duplicate_column_label
        vcf_structure.add_column(duplicate_column, fg_structure.duplicate_indicator)
        vcf_structure.duplicate_rows(duplicate_column)
        vcf_structure.add_column(config.gene_column_label, fg_structure.hgnc_annotation)
        variant_ens = fg_structure.ens_annotation
    else:
        variant_ens = vcf_structure.load_column(config.ensembl_column)

    # Initialize HPO Object
    hpo_structure = hpo.Hpo(hpo_file)
    # Generate the HPO annotation
    hpo_structure.annotate_genes(variant_ens, 'Ensembl')

    # Initialize CDG Object
    cdg_structure = CDG.Cdg()
    # Generate CDG Annotation
    cdg_structure.annotate_cdg(variant_ens, 'Ensembl')

    # Initialize STRING Object
    string_structure = si.String()
    # Generate STRING Annotation
    string_structure.generate_annotations(variant_ens)
    string_structure.annotate_iei(variant_ens, 'Ensembl')

    # Initialize BIOGRID Object
    bio_structure = bp.Bio()
    bio_structure.generate_annotations(variant_ens, 'Ensembl')

    print('\nAnnotating VCF...')
    # Add columns to the VCF
    # HPO
    vcf_structure.add_column('HPO Count', hpo_structure.count_annotation)
    # Check if separating immune HPO IDs
    if config.separate_immune_hpo_terms:
        vcf_structure.add_column('Immune HPO ID', lpl(hpo_structure.immune_id_annotation))
        vcf_structure.add_column('Immune HPO Description', lpl(hpo_structure.immune_description_annotation))
        vcf_structure.add_column('General HPO ID', lpl(hpo_structure.gen_id_annotation))
        vcf_structure.add_column('General HPO Description', lpl(hpo_structure.gen_description_annotation))
    else:
        vcf_structure.add_column('HPO ID', lpl(hpo_structure.gen_id_annotation))
        vcf_structure.add_column('HPO Description', lpl(hpo_structure.gen_description_annotation))

    # Known IEI Genes
    vcf_structure.add_column('Known IEI Gene', string_structure.iei_indicator_annotation)

    # CDG
    vcf_structure.add_column('CDG Interaction', cdg_structure.indicator_annotation)
    vcf_structure.add_column('CDG Data', lpl(cdg_structure.interaction_annotation))

    # STRING
    vcf_structure.add_column('STRING-Interaction', string_structure.indicator_list)
    vcf_structure.add_column('STRING-Data', string_structure.interaction_list)

    # BIOGRID
    vcf_structure.add_column('BioGRID Interaction', bio_structure.indicator_annotation)
    vcf_structure.add_column('BioGRID Data', bio_structure.interaction_annotation)

    # Export to Excel File
    vcf_structure.df_to_excel(config.write_to)
    print('Process Complete, VCF ready!')


def generate_hpo_list():
    hpo_structure = hpo.Hpo(config.hpo_id_filename)
    hpo_structure.hpo_genes_list(config.write_to)
    return


def generate_string_df():
    string_structure = si.String()
    string_structure.count_interactions()


def generate_biogrid_df():
    biogrid_structure = bp.Bio()
    biogrid_structure.generate_interaction_list()


def annotate_iei_list():
    iei_structure = hpo.IeiHpo()
    iei_structure.generate_annotations(config.write_to)


if config.mode == 'annotate_VCF':
    annotate_vcf(config.vcf_filename, config.hpo_id_filename)
elif config.mode == 'generate_HPO_list':
    generate_hpo_list()
elif config.mode == 'generate_STRING_df':
    generate_string_df()
elif config.mode == 'generate_BioGRID_df':
    generate_biogrid_df()
elif config.mode == 'annotate_iei_list':
    annotate_iei_list()
