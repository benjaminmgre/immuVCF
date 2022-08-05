import pandas
import openpyxl


class Vcf:
    """
    VCF Object for parsing, loading, and annotating a VCF file.

    Methods:
    --------
    load_dataframe:
        Loads the VCF in CSV format into a pandas dataframe
    load_column:
        Returns a column from the VCF in the format of a list
    add_column:
        Adds a column to the VCF dataframe
    df_to_csv:
        Loads the current dataframe into an exportable CSV file

    Attributes:
    -----------
    vcf_df: pandas dataframe
        The annotated VCF dataframe
    vcf_len: int
        Number of rows on the VCF
    """

    def __init__(self, f_name: str):
        """
        Parameters
        ---------
        f_name: str
            The CSV VCF handle
        """
        print('\nLoading VCF...')
        self.vcf_df = pandas.read_excel(io=f_name, sheet_name=0)
        self.vcf_len = len(self.vcf_df.index)
        return

    def load_column(self, column_name: str):
        """
        Loads the indexed column name into a list

        Parameters
        ----------
        column_name: string
            Header name of the column in the VCF

        Returns
        -------
        col_list: list (of strings)
            An ordered list of strings of the column
        """

        col_list = []
        if column_name not in self.vcf_df.columns:
            print('Error! Column name not found.')
            return
        for index, row in self.vcf_df.iterrows():
            col_list.append(row[column_name])
        return col_list

    def add_column(self, column_name: str, column_list: list):
        """
        Adds a new column to the VCF

        Parameters
        ----------
        column_name: str
            The column's title
        column_list: list
            Entries for every row of the column. Must be the length of self.vcf_len
        """

        try:
            self.vcf_df[column_name] = column_list
        except IOError:
            print('Error! Check list format.')
        return

    def duplicate_rows(self, duplicate_column_name: str):
        """
        Duplicates selected rows on the VCF based on the duplicate indicator column.

        duplicate_column_name: str
            Name of the duplicate column label on the VCF
        """
        dup_count = 0
        df_copy = self.vcf_df.iloc[:0, :].copy()
        l = len(self.vcf_df)
        for i, row in self.vcf_df.iterrows():
            r = pandas.DataFrame(row).transpose()
            df_copy = pandas.concat([df_copy, r])
            duplicates = int(row[duplicate_column_name])
            for j in range(0, duplicates):
                df_copy = pandas.concat([df_copy, r])
                dup_count += 1
        self.vcf_df = df_copy
        print(f'There were {dup_count} duplicates. The new VCF (with duplicates) is {len(self.vcf_df)} long.')
        return

    def df_to_csv(self, csv_name: str):
        """
        Loads the dataframe into an exportable CSV file

        Parameters
        ----------
        csv_name: string
            Handle where CSV VCF will be written
        """

        self.vcf_df.to_csv(csv_name)
        return

    def df_to_excel(self, excel_name: str):
        """
        Loads the dataframe into a readable Excel file

        Parameters
        ----------
        excel_name: str
            Handle where the Excel file will be written
        """
        self.vcf_df.to_excel(excel_name, sheet_name='immuVCF')
        return
