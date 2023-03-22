import warnings

import pandas as pd
import numpy as np

class metering:
    """This class reads information that is relevant to the flow meter, 
    temperature and pressure transmitters and stores these to its object. 
    For clarity: it does not read information about the stream, such as 
    the flow rate or temperature and pressure. It only reads information 
    specific to the meters, i.e., uncertainties in the measurements, 
    type of flow meter, meter configuration etc.
    
    At the point of reading uncertainties, these are converted to 
    absolute uncertainties. Then a function from met4h2 is used to do 
    unit conversion to SI units."""
    
    def read_general_info(self,xls):
        """Reading general info about the measurement set up."""
        warnings.simplefilter('ignore')
        self.info = pd.read_excel(xls,'Info')
        warnings.resetwarnings() 

    def read_configuration(self,xls):
        """Reading info about the configuration, e.g., meter in series or 
        parallel, etc."""
        warnings.simplefilter('ignore')
        self.configuration = pd.read_excel(xls,'Configuration')
        warnings.resetwarnings()

    def read_flow_meter_uncertainties(
            self,ucalc,stream,m4h2,xls,meter_type,stream_nominal_values):
        """This method reads the measurement uncertainty of the flow meter.
        The method also reads the associated uncertainties of the 
        temperature and pressure transmitters that are assumed connected 
        to the meter.

        Any relative uncertainties are converted to absolute uncertainties. 
        Absolute uncertainties are converted to SI units using a function 
        from the met4h2 object.
        
        The meter information is saved to the metering object in 
        self.transmitter_uncertainties.
        
        Inputs: ucalc:  uncertaintyCalculator object
                stream:                 stream object
                m4h2:                   m4h2 object
                xls:                    excel object
                meter_type:             f.ex 'USM'
                stream_nominal_values:  SI-converted measurements from 
                                        stream
                                        """
            

        # It is discussable whether we need meter_type as input, as this
        # could  be read from self.configuration. But this should be part 
        # of a more comprehensive discussion about the settings file.


        # Reading uncertainties from file. 
        # Are the strings too 'magic'? 
        # Add that to discussion on input file.
        warnings.simplefilter('ignore')
        unc_meter_info = pd.read_excel(xls,meter_type+'_unc')
        unc_temperature = pd.read_excel(xls,'T_unc')
        unc_pressure = pd.read_excel(xls,'P_unc')
        warnings.resetwarnings()

        # Get square sum of meter uncertainty contributions that are 
        # not NaN.
        # Is string too 'magic'? Add that to discussion on input file.
        unc_meter = unc_meter_info[unc_meter_info['Uncertainty'].notnull()]
        unc_info_meter = [
            np.sqrt(sum(np.square(unc_meter['Uncertainty'].to_numpy())))]
        
        # Check that the meter input uncertainty units are consistent.
        # Are strings too 'magic'? Add that to discussion on input file.
        col = ['Unit','k','Distribution']
        for_error_message = ['units','coverage factors','distributions']
        counter = 0
        for c in col:
            if m4h2.all_equal(unc_meter[c].tolist()) == True:
                unc_info_meter.append(unc_meter[c].iloc[0])
            else: raise Exception(
                    'Please make sure that all '+meter_type+' input '
                    + for_error_message[counter]
                    + ' have identical units.')
            counter = counter + 1


        # Add flow meter, P and T, and composition uncertainties 
        # to a common data frame. 
        
        # Note that pressure and temperature uncertainty is read from the 
        # Misc-cell, if this is zero, then that explains the zero std 
        # distribution after monte carlo draws.
        # Are strings too 'magic'? Add that to discussion on input file.
        unc_info_P = [
            unc_pressure['Uncertainty']
            [unc_pressure["Input variable"] == "Misc"].values[0],
            unc_pressure['Unit']
            [unc_pressure["Input variable"] == "Misc"].values[0],
            unc_pressure['k']
            [unc_pressure["Input variable"] == "Misc"].values[0],
            unc_pressure['Distribution']
            [unc_pressure["Input variable"] == "Misc"].values[0]]
        unc_info_T = [
            unc_temperature['Uncertainty']
            [unc_temperature["Input variable"] == "Misc"].values[0],
            unc_temperature['Unit']
            [unc_temperature["Input variable"] == "Misc"].values[0],
            unc_temperature['k']
            [unc_temperature["Input variable"] == "Misc"].values[0],
            unc_temperature['Distribution']
            [unc_temperature["Input variable"] == "Misc"].values[0]]
        unc_info = unc_info_meter + unc_info_P + unc_info_T
        colnames = []
        for m in ['quantity_','pressure_','temperature_']:
            for n in ['','_unit','_coverage','_distribution']:
                colnames.append(m+'unc'+n)
        data = pd.DataFrame([unc_info],columns=colnames)
           
        # Now prepare SI conversion of uncertainties.
        uncertainties = ['quantity_unc','pressure_unc','temperature_unc']
        
        labels = ['quantity','pressure','temperature']
        unc_units = [x+'_unit' for x in uncertainties]

        # Convert any relative uncertainties to absolute uncertainties. 
        stream_nominal_values = stream.process_data.drop(
            columns = ['DateTime','Uptime_h'])
        data = pd.concat([stream.process_data['DateTime'],data],axis=1)
        for i in range(len(uncertainties)):
            if data[[unc_units[i]]].iloc[0,0] == '%':
                nominal_uncertainty = stream_nominal_values.filter(
                    regex='^{}'.format(
                    uncertainties[i].split('_')[0]))*data[uncertainties[i]][0]/100
                data[labels[i]+'_unc'] = nominal_uncertainty
                data[unc_units[i]] = [
                    nominal_uncertainty.columns[0].split('_')[-1]]*len(data)
        data = data.ffill()
        
        # Do SI conversion of absolute uncertainties.     
        self.transmitter_uncertainties = m4h2.SI_conversion(
            ucalc,data,uncertainties,labels)
    
    def calc_composition_uncertainties(self,method='GC'): #### Very tentative.
        """To be defined. method could be 'GC' or 'lab_sample'.
        May want to calculate it given ISO or meter specifications.
        For now it is read from file in read_composition_uncertainties.""" 

    def read_composition_uncertainties(self,xls,stream):
        """This method reads the measurement uncertainty of the gas 
        composition from the xls object.

        Any relative uncertainties are converted to absolute 
        uncertainties. Absolute uncertainties are converted to SI units 
        using a function from the met4h2 object."""   

        warnings.simplefilter('ignore')
        unc_composition = pd.read_excel(xls,'Composition_unc')
        warnings.resetwarnings()

        # Set colnames based on stream.composition, which is a well
        # defined data frame created in metering.read_composition. 
        # Are strings too 'magical'? That is something for the general
        # discussion on input.
        comp_col_names = [
            i+'_unc' for i in stream.composition.drop(
            columns=['DateTime','composition_unit']).columns.tolist()]
        colnames = [
            'composition_unc_unit',
            'composition_unc_coverage',
            'composition_unc_distribution'] + comp_col_names
        
        # Create composition uncertainty data frame and fill columns.
        composition_uncertainties = pd.DataFrame(columns=colnames)  
        
        # Are strings too 'magical'? That is something for the general
        # discussion on the settings file. This sheet certainly has 
        # potential for improvement.
        # Copy the uncertainties from unc_composition to 
        # composition_uncertainties.
        unc_composition = unc_composition.rename(
            columns = dict(zip(
            unc_composition.columns.tolist(),['Comp']+comp_col_names)))
        composition_uncertainties = composition_uncertainties.merge(
            unc_composition[
            unc_composition['Comp']=='Uncertainty'].drop(
            columns =['Comp']),how='right')
        # Get Unit
        unit = unc_composition['H2_unc'][
            unc_composition['Comp']=='Unit'].values[0]
        # If unit is not recognized throw error
        allowed_units = ['cmol_mol-1','\'cmol_mol-1']
        wanted_units = ['mol_mol-1','\'mol_mol-1']
        if unit not in allowed_units+wanted_units:
            raise Exception('Please change the composition'+
                        ' units to cmol_mol-1, mol_mol-1'+
                        'or relative units.')
        # If unit is cmol_mol-1, convert to mol_mol-1
        if unit in allowed_units: 
            composition_uncertainties = composition_uncertainties/100
        # If unit is relative convert to absolute unit
        if unit == '%': 
            """convert to absolute unit"""
            rel_uncertainties = composition_uncertainties.dropna(
                axis=1)
            composition = stream.composition.drop(
                columns=['DateTime','composition_unit'])
            if rel_uncertainties.shape == composition.shape:
                abs_uncertainties = pd.DataFrame(
                    (rel_uncertainties.values/100)*composition.values, 
                    columns=rel_uncertainties.columns, 
                    index=rel_uncertainties.index)
            if len(
                rel_uncertainties) == 1 and rel_uncertainties.shape[
                1] == composition.shape[1]:
                abs_uncertainties = composition.multiply(
                    np.array(rel_uncertainties/100), axis='columns')
            else: raise Exception('Sizes of input composition and',
                                  +'composition uncertainty matrices'
                                  +'do not match.')
            to_rename = dict(zip(
                abs_uncertainties.columns.to_list(),[
                i + '_unc' for i in abs_uncertainties.columns.tolist()]))
            abs_uncertainties = abs_uncertainties.rename(
                columns = to_rename)
            composition_uncertainties = abs_uncertainties.merge(
                composition_uncertainties,how='left')
        # Get distribution
        distribution = unc_composition['H2_unc'][
            unc_composition['Comp']=='Distribution'].values[0]
        composition_uncertainties[
            'composition_unc_distribution'] = [distribution]
        # Get coverage factor
        coverage = unc_composition['H2_unc'][
            unc_composition['Comp']=='k'].values[0]
        composition_uncertainties[
            'composition_unc_coverage'] = [coverage]
        # Add DateTime and unit
        composition_uncertainties[
            'composition_unc_unit'] = ['mol_mol-1']
        self.composition_uncertainites = pd.concat(
            [stream.composition[['DateTime']],
             composition_uncertainties],axis=1)
    
