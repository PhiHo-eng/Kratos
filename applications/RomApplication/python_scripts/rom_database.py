import KratosMultiphysics
import sqlite3
import hashlib
import json
from pathlib import Path
import numpy as np
try:
    import pandas as pd
    from xlsxwriter import Workbook
    missing_pandas = False
    missing_xlsxwriter = False
except ImportError as e:
    if 'pandas' in str(e):
        missing_pandas = True
        missing_xlsxwriter = True
    elif 'xlsxwriter' in str(e):
        missing_xlsxwriter = True
        missing_pandas = False


class RomDatabase(object):

    def __init__(self, general_rom_manager_parameters, mu_names, database_name = "rom_database.db"):
        self.general_rom_manager_parameters = general_rom_manager_parameters
        self.database_root_directory = Path(self.general_rom_manager_parameters["ROM"]["rom_basis_output_folder"].GetString() + "/rom_database"  )
        self.database_name = self.database_root_directory / database_name
        self.npys_directory = self.database_root_directory / "npy_files"
        self.xlsx_directory = self.database_root_directory / "xlsx_files"
        self.create_directories([self.database_root_directory, self.npys_directory, self.xlsx_directory])
        self.mu_names = mu_names   # A single list with the name of each parameter to be used (for reference) in the database, e.g. ["alpha", "permeability", "mach"].
                                   # self.mu_names are expected to match the position and the number of inputs in the lists mu_train, mu_test, mu_run
        self.set_up_or_get_data_base()


    def create_directories(self, directories):
        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)


    def set_up_or_get_data_base(self):
        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()
        table_definitions = {
            "FOM": '''CREATE TABLE IF NOT EXISTS FOM
                        (id INTEGER PRIMARY KEY, parameters TEXT, file_name TEXT)''',
            "ROM": '''CREATE TABLE IF NOT EXISTS ROM
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL,  type_of_projection TEXT,  file_name TEXT)''',
            "HROM": '''CREATE TABLE IF NOT EXISTS HROM
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, type_of_projection TEXT, file_name TEXT)''',
            "HHROM": '''CREATE TABLE IF NOT EXISTS HHROM
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, type_of_projection TEXT, file_name TEXT)''',
            "RightBasis": '''CREATE TABLE IF NOT EXISTS RightBasis
                        (id INTEGER PRIMARY KEY, tol_sol REAL, file_name TEXT)''',
            "LeftBasis": '''CREATE TABLE IF NOT EXISTS LeftBasis
                        (id INTEGER PRIMARY KEY, tol_sol REAL, basis_strategy TEXT, include_phi INTEGER, tol_pg REAL, solving_technique TEXT, monotonicity_preserving INTEGER , file_name TEXT)''',
            "PetrovGalerkinSnapshots": '''CREATE TABLE IF NOT EXISTS PetrovGalerkinSnapshots
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, basis_strategy TEXT, include_phi INTEGER, tol_pg REAL, solving_technique TEXT, monotonicity_preserving INTEGER , file_name TEXT)''',
            "ResidualsProjected": '''CREATE TABLE IF NOT EXISTS ResidualsProjected
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, type_of_projection TEXT, file_name TEXT)''',
            "HROM_Elements": '''CREATE TABLE IF NOT EXISTS HROM_Elements
                        (id INTEGER PRIMARY KEY, tol_sol REAL, tol_res REAL, type_of_projection TEXT, file_name TEXT)''',
            "HROM_Weights": '''CREATE TABLE IF NOT EXISTS HROM_Weights
                        (id INTEGER PRIMARY KEY, tol_sol REAL, tol_res REAL, type_of_projection TEXT, file_name TEXT)'''
        }
        self.table_names = table_definitions.keys()
        for table_name, table_sql in table_definitions.items():
            try:
                cursor.execute(table_sql)
                conn.commit()
                print(f"Table {table_name} created successfully.")
            except sqlite3.OperationalError as e:
                print(f"Error creating table {table_name}: {e}")
        conn.close()


    def hash_parameters(self, *args):
        params_str = '_'.join(str(arg) if isinstance(arg, str) else f"{arg:.10f}" if isinstance(arg, float) else str(arg) for arg in args)
        return hashlib.sha256(params_str.encode()).hexdigest()

    def get_hashed_mu_for_table(self, table_name, mu):
        if self.identify_list_type(mu)=="mu":
            serialized_mu = self.serialize_mu(self.make_mu_dictionary(mu))
        elif self.identify_list_type(mu)=="complete_mu":
            serialized_mu = self.serialize_entire_mu_train(mu)
        else:
            err_msg = f'Error: {self.identify_list_type(mu)}'
            raise Exception(err_msg)
        tol_sol, tol_res, projection_type, pg_data1_str, pg_data2_bool, pg_data3_double, pg_data4_str, pg_data5_bool = self.get_curret_params()
        if table_name == 'FOM':
            hash_mu = self.hash_parameters(serialized_mu)
        elif table_name == 'ROM':
            hash_mu = self.hash_parameters(serialized_mu, tol_sol, projection_type)
        elif table_name == 'HROM':
            hash_mu = self.hash_parameters(serialized_mu, tol_sol, tol_res, projection_type,'HROM')
        elif table_name == 'ResidualsProjected':
            hash_mu = self.hash_parameters(serialized_mu, tol_sol, tol_res, projection_type,'ResidualsProjected')
        elif table_name == 'PetrovGalerkinSnapshots':
            hash_mu = self.hash_parameters(serialized_mu, tol_sol, pg_data1_str,pg_data2_bool,pg_data3_double,pg_data4_str,pg_data5_bool,'PetrovGalerkinSnapshots')
        elif table_name == 'RightBasis':
            hash_mu = self.hash_parameters(serialized_mu, tol_sol)
        elif table_name == 'LeftBasis':
            hash_mu = self.hash_parameters(serialized_mu, tol_sol, pg_data1_str,pg_data2_bool,pg_data3_double,pg_data4_str,pg_data5_bool,'LeftBasis')
        elif table_name == "HROM_Elements":
            hash_mu= self.hash_parameters(serialized_mu, tol_sol,tol_res,projection_type, 'z')
        elif table_name == "HROM_Weights":
            hash_mu= self.hash_parameters(serialized_mu, tol_sol,tol_res,projection_type, 'w')

        return hash_mu, serialized_mu


    def get_curret_params(self):
        tol_sol = self.general_rom_manager_parameters["ROM"]["svd_truncation_tolerance"].GetDouble()
        tol_res =  self.general_rom_manager_parameters["HROM"]["element_selection_svd_truncation_tolerance"].GetDouble()
        projection_type = self.general_rom_manager_parameters["projection_strategy"].GetString()
        pg_data1_str = self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"]["basis_strategy"].GetString()
        pg_data2_bool = self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"]["include_phi"].GetBool()
        pg_data3_double = self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"]["svd_truncation_tolerance"].GetDouble()
        pg_data4_str = self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"]["solving_technique"].GetString()
        pg_data5_bool = self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"]["monotonicity_preserving"].GetBool()
        return tol_sol, tol_res, projection_type, pg_data1_str, pg_data2_bool, pg_data3_double, pg_data4_str, pg_data5_bool



    def check_if_in_database(self, table_name, mu):
        file_name, _ = self.get_hashed_mu_for_table(table_name, mu)
        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()
        cursor.execute(f'SELECT COUNT(*) FROM {table_name} WHERE file_name = ?', (file_name,))
        count = cursor.fetchone()[0]
        conn.close()
        return count > 0, file_name



    def add_to_database(self, table_name, mu, numpy_array ):
        file_name, serialized_mu = self.get_hashed_mu_for_table(table_name, mu)
        tol_sol, tol_res, projection_type, pg_data1_str, pg_data2_bool, pg_data3_double, pg_data4_str, pg_data5_bool = self.get_curret_params()

        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()

        if table_name == 'FOM':
            cursor.execute('INSERT INTO FOM (parameters, file_name) VALUES (?, ?)',
                        (serialized_mu, file_name))
        elif table_name == 'ROM':
            cursor.execute('INSERT INTO ROM (parameters, tol_sol , type_of_projection, file_name) VALUES (?, ?, ?, ?)',
                        (serialized_mu, tol_sol, projection_type, file_name))
        elif table_name == 'HROM':
            cursor.execute('INSERT INTO HROM (parameters, tol_sol , tol_res , type_of_projection, file_name) VALUES (?, ?, ?, ?, ?)',
                        (serialized_mu, tol_sol, tol_res, projection_type, file_name))
        elif table_name == 'ResidualsProjected':
            cursor.execute('INSERT INTO ResidualsProjected (parameters, type_of_projection, tol_sol , tol_res , file_name) VALUES (?, ?, ?, ?, ?)',
                        (serialized_mu, projection_type, tol_sol, tol_res, file_name))
        elif table_name == 'PetrovGalerkinSnapshots':
            cursor.execute('INSERT INTO PetrovGalerkinSnapshots (parameters, tol_sol , basis_strategy, include_phi, tol_pg, solving_technique, monotonicity_preserving, file_name) VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
                        (serialized_mu, tol_sol, pg_data1_str,pg_data2_bool, pg_data3_double, pg_data4_str,  pg_data5_bool, file_name))
        elif table_name == 'RightBasis':
            cursor.execute('INSERT INTO RightBasis (tol_sol, file_name) VALUES (?, ?)',(tol_sol, file_name))
        elif table_name == 'LeftBasis':
            cursor.execute('INSERT INTO LeftBasis (tol_sol, basis_strategy, include_phi, tol_pg, solving_technique, monotonicity_preserving, file_name) VALUES (?, ?, ?, ?, ?, ?, ?)',
                        (tol_sol, pg_data1_str,pg_data2_bool, pg_data3_double, pg_data4_str,  pg_data5_bool, file_name))
        elif table_name == 'HROM_Elements':
            cursor.execute('INSERT INTO HROM_Elements (tol_sol , tol_res , type_of_projection, file_name) VALUES (?, ?, ?, ?)',
                        (tol_sol, tol_res, projection_type, file_name))
        elif table_name == 'HROM_Weights':
            cursor.execute('INSERT INTO HROM_Weights  (tol_sol , tol_res , type_of_projection, file_name) VALUES (?, ?, ?, ?)',
                        (tol_sol, tol_res, projection_type, file_name))
        else:
            err_msg = f'Error: table_name: {table_name} not available. Available options are: {self.table_names}'
            raise Exception(err_msg)

        conn.commit()
        conn.close()

        self.save_as_npy(numpy_array, file_name)



    def identify_list_type(self, lst):
        if not lst:
            return "Empty list"
        if isinstance(lst[0], list):
            if any(isinstance(subitem, list) for item in lst for subitem in item):
                return "Invalid structure - nested lists found"
            return "complete_mu"
        else:
            if any(isinstance(item, list) for item in lst):
                return "Invalid structure - mixed types"
            return "mu"



    def serialize_mu(self, parameters):
        return json.dumps(parameters)


    def save_as_npy(self, result, hash_mu):
        file_path = self.npys_directory / f"{hash_mu}.npy"
        np.save(file_path, result)
        print(f"numpy saved to {file_path}")


    def make_mu_dictionary(self, mu):
        if self.mu_names is None:
            self.mu_names = []
            for i in range(len(mu)):
                self.mu_names.append(f'generic_name_{i}')
        return dict(zip(self.mu_names , mu))


    def get_snapshots_matrix_from_database(self, mu_list, table_name='FOM'):
        unique_tuples = set(tuple(item) for item in mu_list)
        mu_list_unique = [list(item) for item in unique_tuples] #unique members in mu_lust
        SnapshotsMatrix = []
        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()
        for mu in mu_list_unique:
            hash_mu, _ = self.get_hashed_mu_for_table(table_name, mu)
            cursor.execute(f"SELECT file_name FROM {table_name} WHERE file_name = ?", (hash_mu,))
            result = cursor.fetchone()
            if result:
                SnapshotsMatrix.append(self.get_single_numpy_from_database(result[0]))
            else:
                print(f"No entry found for hash {hash_mu}")
        conn.close()
        return np.block(SnapshotsMatrix) if SnapshotsMatrix else None



    def serialize_entire_mu_train(self, mu_train):
        """
        Serializes the entire mu_train list of lists into a string.
        """
        serialized_mu_train = json.dumps(mu_train, sort_keys=True)
        return serialized_mu_train


    def get_single_numpy_from_database(self, hashed_name):
        return np.load(self.npys_directory / (hashed_name + '.npy'))


    def generate_excel(self, full_tables = False, number_of_terms=5):
        if missing_pandas == True:
            KratosMultiphysics.Logger.PrintWarning('\x1b[1;31m[MISSING LIBRARY] \x1b[0m'," pandas library not installed. No excel file generated")
        if missing_xlsxwriter == True:
            KratosMultiphysics.Logger.PrintWarning('\x1b[1;31m[MISSING LIBRARY] \x1b[0m'," xlsxwriter library not installed. No excel file generated")
        if missing_pandas==False and missing_xlsxwriter==False:
            conn = sqlite3.connect(self.database_name)
            query = "SELECT name FROM sqlite_master WHERE type='table';"
            tables = pd.read_sql_query(query, conn)
            if full_tables==False:
                output_file_path = self.xlsx_directory / 'Summary.xlsx'
            elif full_tables==True:
                output_file_path = self.xlsx_directory / 'DataBaseCompleteDump.xlsx'
            with pd.ExcelWriter(output_file_path, engine='xlsxwriter') as writer:
                for table_name in tables['name']:
                    if full_tables==True:
                        df = pd.read_sql_query(f"SELECT * FROM {table_name}", conn)
                    else:
                        df = pd.read_sql_query(f"SELECT * FROM {table_name}  LIMIT {number_of_terms}", conn)
                    df.to_excel(writer, sheet_name=table_name, index=False)
                    print(f"Exported {table_name} to Excel sheet.")
            conn.close()
            print(f"Data exported to Excel file {output_file_path} successfully.")

    def generate_database_summary(self):
        number_of_terms = 5 #TODO user-defined
        self.generate_excel(full_tables=False, number_of_terms=number_of_terms)

    def dump_database_as_excel(self):
        self.generate_excel(full_tables=True)