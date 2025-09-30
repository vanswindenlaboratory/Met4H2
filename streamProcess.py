import warnings
import re

import pandas as pd

class streamProcess:
    """This class reads information about the stream, such as the
    flow rate, temperature and pressure, and composition.

    At the point where the measurements are read, units are converted to
    SI-units by calling a function from m4h2, where it covers the
    conversion.
    """

    def extract_between_brackets(self,text):
        match = re.search(r'\[(.*?)\]', text)
        if match:
            return match.group(1)
        return None

    def read_composition(self,m4h2,xls):
        """This function reads the gas composition from the xls object.
        Allowed input units for the composition are cmol/mol and mol/mol.
        The units of the composition will be changed to mol/mol and stored
        in the stream object."""

        warnings.simplefilter('ignore')
        for i in xls.sheet_names:
            if re.search(r'Composition',i):
                if not re.search(r'Unc',i) and not re.search(r'unc',i):
                    composition = pd.read_excel(xls,i)
        warnings.resetwarnings()

        

        allowed_units = ['cmol_mol-1','\'cmol_mol-1']
        wanted_units = ['mol_mol-1','\'mol_mol-1']

        # Convert the units in composition frame to mol/mol (this is a
        # case not covered by the SI-conversion method).
        # Are strings too 'magical'? Put on discussion for input file.
        component_names = composition.drop(
            columns=['DateTime','Unit']).columns.tolist()
        comp_list = list(m4h2.component_dictionary.keys())
        for component in component_names:
            if 'nC' in component: component = component.replace('nC','n-C')
            if 'iC' in component: component = component.replace('iC','i-C')
            if 'neoC5' in component: component = component.replace('neoC5','neo-C5')
            if 'C6' in component: component = component.replace('C6','n-C6')
            if not component.strip() in comp_list:
                print(component+' is not a recognized component.',
                    'Please add only recognized components.',
                    'The following components are recognized: \n',
                    comp_list)

        composition_shaddow = composition.copy()
        composition_shaddow['Unit'] = 'mol_mol-1'
        composition_shaddow.update(
            composition_shaddow[component_names]/100)
        for i in range(len(composition)):
            if composition["Unit"].iloc[i] not in wanted_units:
                if composition["Unit"].iloc[i] in allowed_units:
                    composition.iloc[i] = composition_shaddow.iloc[i]
                else: raise Exception(
                        'Please change the composition'+
                        ' units to either cmol_mol-1 or mol_mol-1')

        # Also make sure that the composition column names are nice.
        composition = composition.rename(
            columns = {'Unit':'composition_unit'})
        new_column_names = []
        for x in component_names:
            new_column_names.append(
                x.strip().replace(
                'nC','n-C').replace(
                'iC','i-C').replace(
                'neoC5','neo-C5').replace(
                'C6','n-C6'))
        self.composition = composition.rename(
            columns =dict(zip(component_names,new_column_names)))

        # Align time stamps.
        aligned_df = pd.merge(self.process_data,self.composition,on='DateTime',how='outer')
        aligned_df = aligned_df.reset_index()
        aligned_df = aligned_df.sort_values(by=['DateTime'])
        aligned_df[['composition_unit']+new_column_names] = aligned_df[['composition_unit']+new_column_names].apply(lambda col: col.ffill())
        aligned_df = aligned_df.dropna(subset=['composition_unit']+new_column_names)
        aligned_df = aligned_df.dropna(subset=self.process_data.columns)
        aligned_df = aligned_df.reset_index(drop=True)
        self.composition = aligned_df.drop(columns=self.process_data.columns[1:])
        self.composition = self.composition.drop(columns=['index'])

    def set_composition(self,xls=None,read_from_file=True):
        """An unfinished method. Is supposed to provide alternative ways
        to get the composition."""
        if read_from_file == True:
            if isinstance(xls,pd.io.excel._base.ExcelFile):
                self.read_composition(xls)
            else:
                raise Exception(
                    'Provide a pandas.io.excel._base.ExcelFile object'
                    +'if you want to read compostions from file.')

    def read_process_data(self,ucalc,m4h2,xls):
        """This function needs access to the ucalc and m4h2 object in order
        to do SI conversion. It also needs access to the xls object in order
        to read the process information (quantity or rate, and pressure and
        temperature).The process information is saved to the object in a
        data frame.
        """
        warnings.simplefilter('ignore')
        for i in xls.sheet_names:
            if re.search(r'ProcessData',i):
                process_data = pd.read_excel(xls,i)
        warnings.resetwarnings()

        ## Do uptime h multiplication if it is a rate that is measured.
        ## I.e., output will not be rate, but total quantity since last
        ## measurement.
        #uptime_multiply = ['Volume flow rate [m3/h]',
        #                'Standard volume flow rate [Sm3/h]',
        #                'Normal volume flow rate [Nm3/h]',
        #                'Mass flow rate [kg/h]']
        #for col_name in process_data.columns:
        #    if col_name in uptime_multiply:
        #        process_data[col_name] = process_data[
        #            col_name]*process_data['Uptime [h]']
        #        process_data = process_data.rename(
        #            columns={col_name:
        #                    re.sub('/[^]]+','',col_name).replace(
        #                        'flow rate ','')})
        #        print('Input is given as a flow rate. ',
        #            'Output will be amount of gas.')

        # prepare and then convert to SI.
        nominal_df = process_data.drop(
            columns=['DateTime']).filter(regex ='^((?!ptime).)*$' )

        self.input_unit = self.extract_between_brackets(nominal_df.columns[0])
        nominal_df = nominal_df.rename(columns={nominal_df.columns[0]:'['+self.input_unit+']'})
        if 'kg' not in self.input_unit:
            if "h" in self.input_unit:
                new_unit = 'Volume flow rate '+nominal_df.columns[0].replace(self.input_unit,'m3/h')
            else:
                new_unit = 'Volume '+nominal_df.columns[0].replace(self.input_unit,'m3')
            nominal_df.rename(columns={nominal_df.columns[0]:new_unit}, inplace=True)
        nominals = nominal_df.columns.tolist()
        process_data.rename(columns={process_data.columns[1]:nominals[0]}, inplace=True)
        
        labels = ['quantity','temperature','pressure']
        process_data = m4h2.SI_conversion(
            ucalc,process_data,nominals,labels,update_colname=True)

        # Make the nominal column values more generic for easier indexing
        # of data frame later on.
        quantities = [
            'Volume [m³]',
            'Standard volume [m³]',
            'Normal volume [m³]',
            '[kg]',
            'Volume flow rate [m³/s]',
            'Standard volume flow rate [m³/s]',
            'Normal volume flow rate [m³/s]',
            '[kg/s]']
        generic_quantities = [
            'quantity_m³',
            'quantity_Sm³',
            'quantity_Nm³',
            'quantity_kg',
            'quantity_(m³)/(s)',
            'quantity_(Sm³)/(s)',
            'quantity_(Nm³)/(s)',
            'quantity_(kg)/(s)']
        for i in range(len(process_data.columns)):
            col_name = process_data.columns[i]
            if col_name in quantities:
                idx = quantities.index(col_name)
                process_data.rename(
                    columns={col_name:generic_quantities[idx]},inplace=True)
        process_data = process_data.rename(
            columns={'P [Pa]':'pressure_Pa','T [K]':'temperature_K'})

        # Legge enhet på self. og i regningen lenger inne i programmet fikse omregning hvis Nm eller Sm
        self.process_data = process_data
