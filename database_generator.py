#################################################################
# Code to load initial target data from csv and generate clean  #
# SQL database with initial transit forecast inside, but no     #
# telescope information                                         #
# Called by run_sim.py                                          #
#                                                               #
# Hamish Caines 07-2019                                         #
#################################################################


def load_data_from_csv(filename):
    """
    Read csv file and load data, storing data as a list of lists, and storing column headers as a list of strings
    :param filename: location of file to be read: str
    :return columns: list of column headers
    :return rows: list of data entries
    """
    rows = []
    try:
        print(filename)
        file = open(filename, 'r')
        # read column headers from first line
        line = file.readline().split('\n')[0].split(',')
        columns = line
        # read data from the remainder of the file and add to rows
        while True:
            line = file.readline().split('\n')[0].split(',')
            if line == ['']:
                break
            rows.append(line)
        return columns, rows

    except FileNotFoundError:
        raise FileNotFoundError


def obtain_types(values):
    """
    Obtain data types for columns from an array of data
    :param values: Row of data from csv file
    :return: Array of types
    """
    types = []
    # loop through array of values, except for final column
    for value in values[:-1]:
        if value != '':
            try:
                float(value)  # check for numeric value
                types.append('DECIMAL(16,8)')
            except ValueError:
                types.append('VARCHAR(20)')
        else:
            types.append('DECIMAL(16,8)')
    # final column for targets is boolean <------ Need to fix for reading telescopes, currently giving aperture boolean
    types.append('BOOLEAN')
    return types


def table_string_builder(columns, types, table_name):
    """
    Build SQL statement to be executed to create a table
    :param columns: List of column headers
    :param types: List of column types
    :param table_name: Name of table: str
    :return: SQL statement to be executed: str
    """
    string = 'CREATE TABLE IF NOT EXISTS ' + table_name + '('
    # add each name then type
    for i in range(0, len(columns)):
        string += columns[i] + ' ' + types[i] + ', '
    string = string[:-2]
    string += ', UNIQUE (Name))'
    return string


def data_string_builder(columns, data, table_name):
    """
    Build SQL statement to be execute to insert a new data row
    :param columns: List of column headers
    :param data: List of data arrays
    :param table_name: Name of table to be inserted into
    :return: SQL statement to be executed: str
    """
    string = 'INSERT INTO ' + table_name + ' ('
    # add columns
    for column in columns:
        string += column + ', '
    string = string[:-2] + ') VALUES ('
    # add values, inserting escape characters for alphanumeric strings
    for value in data:
        if value == '':
            string += 'NULL, '
        elif value[0].isdigit():
            if not value.isalnum():
                string += '\"' + value + '\", '
            else:
                string += '\'' + value + '\', '
        else:
            string += '\'' + value + '\', '
    string = string[:-2] + ')'
    return string


def generate_sql_table_from_csv(filename, table_name, cursor):
    """
    Generate SQL table populated with data read from csv file
    :param filename: Location of csv file to be read: str
    :param table_name: Name of table: str
    :param cursor: Cursor connected to database
    :return: list of data arrays
    """
    import sqlite3
    columns, data = load_data_from_csv(filename)
    types = obtain_types(data[0])
    table_string = table_string_builder(columns, types, table_name)
    cursor.execute(table_string)

    for single in data:
        data_string = data_string_builder(columns, single, table_name)
        try:
            cursor.execute(data_string)
        except sqlite3.IntegrityError:
            pass

    return data


def main():
    import sqlite3
    database_name = 'clean.db'
    start_file = 'realfake_data2.csv'

    db = sqlite3.connect(database_name)
    cursor = db.cursor()
    generate_sql_table_from_csv(start_file, 'TARGET_DATA', cursor)
    db.commit()

    cursor.execute('SELECT Name FROM TARGET_DATA')
    names = cursor.fetchall()
    for name in names:
        cursor.execute('CREATE TABLE IF NOT EXISTS \'' + name[
            0] + '\'(ObID REAL, Epoch REAL, ObsCenter DECIMAL(16,8), ObsCenterErr DECIMAL(16,8), TrueCenter DECIMAL(16,'
                 '8), TrueCenterErr DECIMAL(16,8), ObsDepth DECIMAL(16,8), ObsDepthErr DECIMAL(16,8), ObsDuration DECIM'
                 'AL(16,8), ObsDurationErr DECIMAL(16,8), Source VARCHAR(25), UNIQUE(ObsCenter))')
    db.commit()

    import actions
    from datetime import datetime
    database = actions.Database(database_name, None, None, None)
    database.run_queries()



    start_date = datetime.today()
    #database.transit_forecast(start_date, start_date+timedelta(days=28))

    print('Database initialised')


if __name__ == '__main__':
    main()
