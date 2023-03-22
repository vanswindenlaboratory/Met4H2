import warnings

import pandas as pd

class streamProcess:
    """This class reads information about the stream, such as the 
    flow rate, temperature and pressure, and composition. 
    
    At the point where the measurements are read, units are converted to 
    SI-units by calling a function from m4h2, where it covers the 
    conversion.
    """
    
    def read_composition(self,xls):
        """This function reads the gas composition from the xls object.
        Allowed input units for the composition are cmol/mol and mol/mol.
        The units of the composition will be changed to mol/mol and stored
        in the stream object."""

        warnings.simplefilter('ignore')
        composition = pd.read_excel(xls,'Composition')
        warnings.resetwarnings()

        allowed_units = ['cmol_mol-1','\'cmol_mol-1']
        wanted_units = ['mol_mol-1','\'mol_mol-1']
        
        # Convert the units in composition frame to mol/mol (this is a 
        # case not covered by the SI-conversion method).
        # Are strings too 'magical'? Put on discussion for input file.
        component_names = composition.drop(
            columns=['DateTime','Unit']).columns.tolist()
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
        to do SI conversion. It also needs acces to the xls object in order 
        to read the process information (quantity or rate, and pressure and 
        temperature).The process information is saved to the object in a 
        data frame.
        """
        warnings.simplefilter('ignore')
        process_data = pd.read_excel(xls,'ProcessData')
        warnings.resetwarnings()

       
        # Are strings to 'magical'? Discuss in settings file discussion.
        process_data = process_data.drop(columns=['Case'])
        # Do uptime h multiplication if it is not a rate that is measured.
        uptime_multiply = ['Volume_m3','StandardVolume_Sm3','Mass_kg']
        for col_name in process_data.columns:
            if col_name in uptime_multiply:
                process_data[col_name] = process_data[
                    col_name]*process_data['Uptime_h']
            

        # prepare and then convert to SI.
        nominal_df = process_data.drop(
            columns=['DateTime']).filter(regex ='^((?!ptime).)*$' )
        nominals = nominal_df.columns.tolist()
        labels = ['quantity','temperature','pressure']
        process_data = m4h2.SI_conversion(
            ucalc,process_data,nominals,labels,update_colname=True)
        
        # Make the nominal column values more generic for easier indexing 
        # of data frame later on.
        quantities = [
            'Volume_m³',
            'StandardVolume_m³',
            'Mass_kg',
            'VolumeFlowRate_m³/s',
            'StandardVolumeFlowRate_m³/s',
            'MassFlowRate_kg/s']
        generic_quantities = [
            'quantity_m³',
            'quantity_Sm³',
            'quantity_kg',
            'quantity_(m³)/(s)',
            'quantity_(Sm³)/(s)',
            'quantity_(kg)/(s)']       
        for i in range(len(process_data.columns)):
            col_name = process_data.columns[i]
            if col_name in quantities:   
                process_data.rename(
                    columns={col_name:generic_quantities[i]},inplace=True)      
        process_data = process_data.rename(
            columns={'P_Pa':'pressure_Pa','T_K':'temperature_K'})
        
        self.process_data = process_data
