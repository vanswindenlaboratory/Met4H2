import warnings
import re

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
        for i in xls.sheet_names:
            if re.search(r'Info',i):
                self.info = pd.read_excel(xls,i)
        warnings.resetwarnings()

    def read_configuration(self,xls):
        """Reading info about the configuration, e.g., meter in series or
        parallel, etc."""
        warnings.simplefilter('ignore')
        for i in xls.sheet_names:
            if re.search(r'Configuration',i):
                self.configuration = pd.read_excel(xls,i)
        idx = self.configuration.loc[
            self.configuration['Configuration'] == 'Quantity meter'].index[0]
        self.meter_type = self.configuration['Setting'][idx]
        warnings.resetwarnings()

        idx = self.configuration.loc[
            self.configuration['Configuration'] == 'Number of Monte Carlo draws'].index[0]
        n_monte_carlo_draws = self.configuration['Setting'][idx]
        if n_monte_carlo_draws > 0:
            self.n_monte_carlo_draws = n_monte_carlo_draws
            self.computation_method = 'montecarlo'
            self.nominal = False
        else:
            self.computation_method = 'nominal'
            self.nominal = True

        idxs = ['EoS',
                'EoS tool',
                'Project name',
                'Project info',
                'Meter configuration',
                'No of meters',
                'Correlation coefficient metering',
                'Process analyzer',
                'Performance requirements',
                'Pipe diameter',
                'Base conditions',
                'Composition uncertainty',
                'Normalize all components',
                'By difference component analysis',
                'Gross calorific value temperature',
                'Gross calorific value H2 only']
        variables_to_fill = list()
        for i,k in enumerate(idxs):
            id = self.configuration.loc[
                self.configuration['Configuration'] == k].index[0]
            variables_to_fill.append(self.configuration['Setting'][id])
        (self.EoS,
        self.framework,
        self.project_name,
        self.project_info,
        self.meter_configuration,
        self.no_of_meters,
        self.corr_coeff_met,
        self.process_analyzer,
        self.performance_requirements,
        self.pipe_diameter,
        self.base_conditions,
        self.composition_unc,
        self.norm_component,
        self.by_difference,
        self.gcv_temperature,
        self.gcv_H2) =  tuple(variables_to_fill)
        if self.norm_component == 'no':
            self.normalize_all = False
        else: self.normalize_all = True

        if self.base_conditions == 'Normal':
            self.output_volume_units = 'Nm3'
        elif self.base_conditions == 'Standard':
            self.output_volume_units = 'Sm3'


        allowed_combinations = [('GERG-2008','refprop'),
                                ('GERG-2008','pyforfluids'),
                                ('GERG-2008','TREND'),
                                ('AGA-8','TREND'),
                                ('AGA-8','AGA8Detail')]
        combination = (self.EoS,self.framework)
        if not combination in allowed_combinations:
            raise Exception('The combination of framework and EoS is not supported.',
                            'Please use one of the following combinations: \n',
                            allowed_combinations)

    def get_transmitter_unc_old(self,df):
        """"""

        # Get meter uncertainty
        idx = df['Input variable'][
            df['Input variable'] == 'Uncertainty level'].index[0]
        if df['Input value'].loc[idx] == 'Overall':
            distribution = df['Distribution'].loc[idx+1]
            unit = df['Unit'].loc[idx+1]
            confidence = df['Confidence interval'].loc[idx+1]
            if distribution == 'Rectangular':
                if df['Confidence interval'].loc[idx+1]*100 != 100:
                    print(str(df['Confidence interval'].loc[idx+1]*100),
                        'confidence interval not supported for a',
                        df['Distribution'].loc[idx+1],
                        'distribution. Using 100 % confidence',
                        'instead.')
                df.loc[4,idx+1] = 1
            if distribution == 'Triangular':
                #dvd = np.sqrt(6)
                if df['Confidence interval'].loc[idx+1]*100 != 100:
                    print(str(df['Confidence interval'].loc[idx+1]*100),
                        'confidence interval not supported for a',
                        df['Distribution'].loc[idx+1],
                        'distribution. Using 100 % confidence',
                        'instead.')
                df['Confidence interval'].loc[idx+1] = 1
            lst = [df['Uncertainty'].loc[idx+1]]#[2*df['Uncertainty'][idx+1]/dvd]
            lst.append(unit)
            lst.append(confidence)
            lst.append(distribution)
        else:
            raise Exception('Detailed uncertainty level for '
                            +'flow meter uncertainty is not yet supported.')

        return lst

    def old_way_to_convert_data(self,stream,unc_meter_info,unc_temperature,unc_pressure):
        unc_info_meter = self.get_transmitter_unc_old(unc_meter_info)
        unc_info_P = unc_pressure.iloc[0,:].to_list()[1:-1]
        unc_info_T = unc_temperature.iloc[0,:].to_list()[1:-1]

        unc_info = unc_info_meter + unc_info_P + unc_info_T
        colnames = []
        for m in ['quantity_','pressure_','temperature_']:
            for n in ['','_unit','_confidence','_distribution']:
                colnames.append(m+'unc'+n)
        data = pd.DataFrame([unc_info],columns=colnames)

        # Now prepare SI conversion of uncertainties.
        uncertainties = ['quantity_unc','pressure_unc','temperature_unc']

        labels = ['quantity','pressure','temperature']
        unc_units = [x+'_unit' for x in uncertainties]

        # Convert any relative uncertainties to absolute uncertainties.
        stream_nominal_values = stream.process_data.drop(
            columns = ['DateTime','Uptime [h]'])
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
        return (data, uncertainties, labels)


    def find_closest_larger_index(self, df, column_name, n):
        """Finds the index of the row in the specified column that is closest to n 
            among the numbers larger than n.

        Args:
            df: The pandas DataFrame.
            column_name: The name of the column to search in.
            n: The target number.

        Returns:
            The index of the closest larger number, or None if no such number is found.
        """

        # Filter values larger than n
        larger_values = df[df[column_name] > n][column_name]

        # Check if there are any larger values
        if larger_values.empty:
            return None  # Return None if no larger values are found

        # Find the index of the closest larger value
        closest_larger_index = larger_values.idxmin()  # Efficiently using idxmin for a Series

        return closest_larger_index


    def read_transmitter_uncertainties(
            self,ucalc,stream,m4h2,xls,Z,Z0,meter_type=None):
        """
        The method reads uncertainties of the
        quantity, temperature, and pressure transmitters that are assumed connected
        to the meter.

        Any relative uncertainties are converted to absolute uncertainties.
        Absolute uncertainties are converted to SI units using a function
        from the met4h2 class.

        The meter information is saved to the metering object in
        self.PT_uncertainties.

        Inputs: ucalc:                  uncertaintyCalculator object
                stream:                 stream object
                m4h2:                   m4h2 object
                xls:                    excel object
                meter_type:             f.ex 'USM'
                stream_nominal_values:  SI-converted measurements from
                                        stream
                                        """


        if meter_type == None:
            try:
                meter_type = self.configuration['Setting'][
                    self.configuration['Configuration']=='Type of meter'][0]
            except: raise Exception('Please specify meter type when calling.',
                                    'read_settings_file()')


        # Reading uncertainties from file.
        warnings.simplefilter('ignore')
        for i in xls.sheet_names:
            if re.search(meter_type+'_unc',i):
                unc_meter_info = pd.read_excel(xls,i)
            if re.search('T_unc',i):
                unc_temperature = pd.read_excel(xls,i)
            if re.search('P_unc',i):
                unc_pressure = pd.read_excel(xls,i)
        warnings.resetwarnings()

        stream_nominal_values = stream.process_data.copy()

        ### Old set-up. Still possible to call function but
        # have written a new_one.
        #data,uncertainties,labels = self.old_way_to_convert_data(
        #    stream,unc_meter_info,unc_temperature,unc_pressure)

        uncertainties = ['quantity_unc','pressure_unc','temperature_unc']
        labels = ['quantity','pressure','temperature']
        unc_units = [x+'_unit' for x in uncertainties]

        ## New set-up.
        cols = ['Input variable',
                'range',
                '',
                '_unit',
                '_confidence',
                '_distribution',
                'comment']
        colnames = list()
        for i, var in enumerate(['quantity_unc','pressure_unc','temperature_unc']):
            colns = list()
            for col in cols:
                colns.append(var + col)
            colnames.append(colns)
        unc_meter_info.columns = colnames[0]
        unc_pressure.columns = colnames[1]
        unc_temperature.columns = colnames[2]

        for i, df in enumerate([unc_meter_info,unc_pressure,unc_temperature]):
            select_confidence = [x for x in df.columns if 'confidence' in x]
            select_distribution = [x for x in df.columns if 'distribution' in x]
            select_range = [x for x in df.columns if 'range' in x]
            if i == 0:
                if 's' in stream_nominal_values.columns[1]:
                    df[select_range] = df[select_range]/3600
                    if 'kg' not in stream_nominal_values.columns[1]:
                        if self.base_conditions == 'Normal' or self.base_conditions == 'Standard':
                            norm_flow = stream_nominal_values.iloc[:,1].copy()
                            stream_nominal_values.iloc[:,1] = (stream_nominal_values.iloc[:,1].values*m4h2.P0*stream_nominal_values['temperature_K']*Z.values[:,0])/(stream_nominal_values['pressure_Pa']*m4h2.T0*Z0.values[:,0])
            for _, row in df.iterrows():
                if row[select_distribution][0] in ['Rectangular', 'Triangular'] and row[select_distribution][0] != 1:
                    print(f"Error: {row[select_confidence][0]*100} % confidence interval not supported for a {row[select_confidence][0]} distribution. Using 100 % confidence instead.")
                    row[select_confidence][0] = 1
            select_label = [x for x in stream_nominal_values.columns if labels[i] in x]
            stream_nominal_values[colnames[i][2:6]] = None
            for j, row in stream_nominal_values.iterrows():
                n = row[select_label][0]
                idx = self.find_closest_larger_index(df, select_range, n)
                stream_nominal_values.loc[j,colnames[i][2:6]] = df[colnames[i][2:6]].iloc[idx[0]]

        stream_nominal_values.iloc[:,1] = norm_flow.copy()
        data = stream_nominal_values

        # Convert any relative uncertainties to absolute uncertainties.
        #data = pd.concat([stream.process_data['DateTime'],data],axis=1) # old
        for i in range(len(uncertainties)):
            for t in range(len(data)):
                if data[[unc_units[i]]].iloc[t,0] == '%':
                    col = stream_nominal_values.filter(
                        regex='^{}(?!.*_unc)'.format(
                        uncertainties[i].split('_')[0]))
                    nominal_uncertainty = col.iloc[t,0]*data[uncertainties[i]][t]/100
                    data[labels[i]+'_unc'][t] = nominal_uncertainty
                    data[unc_units[i]][t] = col.columns[0].split('_')[-1]
        data = data.ffill()
        data = data.drop(columns = stream.process_data.columns[1:])


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
        uncertainties. Absolute uncertainties are converted to SI units."""

        warnings.simplefilter('ignore')
        for i in xls.sheet_names:
            if re.search('Composition_unc',i):
                unc_composition = pd.read_excel(xls,i)
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
            'composition_unc_confidence',
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
        unit = unc_composition.iloc[:,1][
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
        distribution = unc_composition[unc_composition.columns[1]][
            unc_composition['Comp']=='Distribution'].values[0]
        composition_uncertainties[
            'composition_unc_distribution'] = [distribution]
        # Get confidence interval
        confidence = unc_composition[unc_composition.columns[1]][
            unc_composition['Comp']=='Confidence interval'].values[0]
        composition_uncertainties[
            'composition_unc_confidence'] = [confidence]
        # Add DateTime and unit
        composition_uncertainties[
            'composition_unc_unit'] = ['mol_mol-1']
        self.composition_uncertainties = pd.concat(
            [stream.composition[['DateTime']],
            composition_uncertainties],axis=1)
