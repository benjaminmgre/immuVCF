"""
Functions for formatting lists to be printed on the VCF file.
"""


def print_list(lst):
    """
    Takes a list, and returns a string with every entry in the list.
    """
    lst = [str(i) for i in lst]
    print_string = ', '.join(lst)
    return print_string


def list_print_list(master_list):
    """
    Takes a list of lists, and returns a list of strings to be printed on a dataframe column.
    """
    string_list = []
    for lst in master_list:
        if len(lst) == 0:
            print_string = ''
        else:
            lst = [str(i) for i in lst]
            print_string = ', '.join(lst)
            # print(print_string)
        # print(print_string)
        string_list.append(print_string)
    return string_list
