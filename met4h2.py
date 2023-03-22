import time
import os 
import importlib.util 
import sys 

import numpy as np
import pandas as pd
import seaborn as sns

import mcerp as mc 

from trend import trend
from streamProcess import streamProcess
from metering import metering

from modules.uncertaintyCalculator import uncertaintyCalculator



from modules.molecule import *
from modules.AGA8Detail import AGA8Detail




class met4H2:
    """This class has methods that can collect info from streamProcess and
     metering objects. From there it can do monte carlo draws, and can 
     compute gas properties by using either the TREND or AGA8Detail 
     frameworks. The gas calculated gas properties are then stored in 
     self, and finally energy or energy flow of the stream can be 
     calculated by using the uncertaintyCalculator object."""

    def __init__(self):
        self.all_measurement_input = pd.DataFrame()
        self.monte_carlo_draws = {}
        self.gcv_temperature = 25    # Temperature at which to get gross
                                        # calorific value.
        self.T0 = np.array(15+273.15) 
        self.P0 = np.array(101325)   # AGA8Detail can only have numpy as input.
        self.model_Z_noise = 0.001   # Model noise could be added here, or 
        self.model_Hmol_noise = 0    # defined in settings file. It should be
        self.model_Z0_noise = 0      # given as standard uncertainty.
        self.model_MM_noise = 0
        self.model_DM_noise = 0
        self.model_DM0_noise = 0
        self.component_dictionary = dict(zip(['N2','CO2','C1','C2','C3',
                                              'n-C4','i-C4','n-C5','i-C5',
                                              'neo-C5','n-C6','2-Methylpentane',
                                              '3-Methylpentane','2,2-Dimethylbutane',
                                              '2,3-Dimethylbutane','n-C7','n-C8',
                                              'n-C9','n-C10','n-C11','n-C12','n-C13',
                                              'n-C14','n-C15','Ethylene','Propene',
                                              '1-Butene','cis-2-Butene',
                                              'trans-2-Butene','2-Methylpropene',
                                              '1-Pentene','Propadiene',
                                              '1,2-Butadiene','1,3-Butadiene',
                                              'Acetylene','Cyclopentane',
                                              'Methylcyclopentane',
                                              'Ethylcyclopentane','Cyclohexane',
                                              'Methylcyclohexane','Ethylcyclohexane',
                                              'Benzene','Toluene','Ethylbenzene',
                                              'o-Xylene','Methanol','Methanetiol',
                                              'H2','H2O','H2S','NH3','HCN','CO',
                                              'COS','CS2'],
                                             [['Nitrogen','nitrogen'],
                                              ['CO2','co2'],
                                              ['Methane','methane'],
                                              ['Ethane','ethane'],
                                              ['Propane','propane'],
                                              ['n-Butane','butane'],
                                              ['Isobutane','isobutan'],
                                              ['n-Pentane','pentane'],
                                              ['Isopentane','ipentane'],
                                              ['Neopentane','neopentn'],
                                              ['n-Hexane','hexane'],
                                              ['2-Methylpentane'],
                                              ['3-Methylpentane'],
                                              ['2,2-Dimethylbutane'],
                                              ['2,3-Dimethylbutane'],
                                              ['n-Heptane'],
                                              ['n-Octane'],
                                              ['n-Nonane'],
                                              ['n-Decane'],
                                              ['n-Undecane'],
                                              ['n-Dodecane'],
                                              ['n-Tridecane'],
                                              ['n-Tetradecane'],
                                              ['n-Pentadecane'],
                                              ['Ethene'],
                                              ['Propene'],
                                              ['1-Butene'],
                                              ['cis-2-Butene'],
                                              ['trans-2-Butene'],
                                              ['2-Methylpropene'],
                                              ['1-Pentene'],
                                              ['Propadiene'],
                                              ['1,2-Butadiene'],
                                              ['1,3-Butadiene'],
                                              ['Ethyne'],
                                              ['Cyclopentane'],
                                              ['Methylcyclopentane'],
                                              ['Ethylcyclopentane'],
                                              ['Cyclohexane'],
                                              ['Methylcyclohexane'],
                                              ['Ethylcyclohexane'],
                                              ['Benzene'],
                                              ['Toluene'],
                                              ['Ethylbenzene'],
                                              ['o-Xylene'],
                                              ['Methanol'],
                                              ['Methanetiol'],
                                              ['Hydrogen','hydrogen'],
                                              ['Water'],
                                              ['Hydrogen sulfide'],
                                              ['Ammonia'],
                                              ['Hydrogen cyanide'],
                                              ['Carbon monoxide'],
                                              ['Carbonyl sulfide'],
                                              ['Carbon disulfide']]))

    def all_equal(self,iterator):
        """Checking that all items in iterator are identical.
        Returns true if that is the case."""
        iterator = iter(iterator)
        try:
            first = next(iterator)
        except StopIteration:
            return True
        return all(first == x for x in iterator)
      
    def SI_conversion(
            self,ucalc,data,colnames,labels,update_colname=False):
        """The SI converter in uncertaintyCalculator only accepts inputs 
        in a specific dictionary format. These dictionaries are created for m4h2
        here, and passed to the uncertaintyCalculator.
        
        Output: data frames with SI-converted measurements or uncertainties.
        Input:  ucalc:      uncertaintyCalculator object
                data:       data frame containing values to be converted
                colnames:   column names in the data frame containing the values
                labels:     labels of values in colnames
                """

        initialist = []
        unit_dict = {
            'C':'°C','F':'°F','Mpa':'MPa','barg':'bar',
            'm3h-1':'(m^3)/(h)','Sm3h-1':'(m^3)/(h)',
            'm3h-1':'(m^3)/(h)','m3':'m^3','Sm3':'m^3',
            'm³/s':'(m^3)/(s)','Sm³/s':'(m^3)/(s)',
            'kg/s':'(kg)/(s)','Sm^3/s':'(m^3)/(s)','Sm³':'m^3',
            }
        
        for i in range(len(labels)): initialist.append([])
        
        unit_dictionary = dict(zip(labels,initialist))
        units = []
        for t in range(len(data)):
            standard = []
            for i in range(len(colnames)):
                if update_colname == True:
                    unit = colnames[i].split('_')[-1]
                else: unit = data[colnames[i]+'_unit'][t]
                units.append(unit)
                # When looking for Sm3 this is just to flag it as standard m3
                # and not do more about it, since it should be treated the same
                # way as the SI converter in ucalc treats volumes.
                if 'Sm3' in unit or 'Sm³' in unit: standard.append(True)
                else: standard.append(False)
                if unit in unit_dict.keys():
                    unit = unit_dict[unit]
                tuple_list = list(zip(data[colnames[i]].tolist(),len(data)*[unit]))
                unit_dictionary[labels[i]] = [list(x) for x in tuple_list]
        ucalc.input_dictionary['key'] = unit_dictionary
        ucalc.SI_conversion(input_keys=['key'])
        # Put the converted units into the data frame and rename columns accordingly.
        for no in range(len(colnames)):
            data[colnames[no]] = [
                x[0] for x in ucalc.input_dictionary['key'][labels[no]]]
            if update_colname == True:
                data = data.rename(
                    columns={colnames[no]:colnames[no].replace(
                    units[no],ucalc.input_dictionary['key'][labels[no]][0][1])})
                if standard[no] == True:
                    data = data.rename(
                        columns={colnames[no]:colnames[no].replace('m³','Sm³')})
            else:
                data[colnames[no]+'_unit'] = len(data)*[
                    ucalc.input_dictionary['key'][labels[no]][0][1]]
                if standard[no] == True:
                    data[colnames[no]+'_unit'] = len(data)*[
                        ucalc.input_dictionary['key'][labels[no]][0][1].replace(
                        'm³','Sm³')]
        return data

    def read_settings_file(
            self,ucalc, stream, meter, filename, 
            filepath = None, meter_type = 'USM'):
        """Reads settings file and stores the information to data frames 
        stored in the m4h2 object.All units are converted to SI units.
        
        Inputs: ucalc:  uncertaintyCalculator object
                stream: streamProcess object
                meter:  metering object
                filename:   f.ex 'settings_file.xlsx'
                filepath:   replace with os"""

        if filepath == None:
            xls = pd.ExcelFile(filename)
        else:
            xls = pd.ExcelFile(os.path.join(filepath,filename))

        # Getting information about the stream, i.e., composition, 
        # quantity or rate and pressure and temperature at line 
        # conditions.
        stream.read_composition(xls)
        stream.read_process_data(ucalc,self,xls)
        
        # meter.dataframes should be set after stream. 
        meter.read_general_info(xls)
        meter.read_flow_meter_uncertainties(
            ucalc,stream,self,xls,meter_type,stream.process_data)
        meter.read_composition_uncertainties(xls,stream)
        meter.read_configuration(xls)
        
    def collect_measurement_input(self,stream,meter):
        """Information about stream and meter are saved in their
        respective objects. Given the objects, this function
        collects this info in a common data frame and saves to 
        met4H2.all_measurement_input."""

        # Merge the all_measurement_input and composition frame 
        # on DateTime. Fill na with previous entry where needed.
        self.all_measurement_input = stream.process_data.merge(
            meter.transmitter_uncertainties,on='DateTime')
        comp = stream.composition.merge(
            meter.composition_uncertainites,on='DateTime')
        self.all_measurement_input = self.all_measurement_input.merge(
            comp,left_on='DateTime',right_on='DateTime',how='outer')
        self.all_measurement_input = self.all_measurement_input.sort_values(
            by='DateTime',ignore_index=True).fillna(
            method="ffill").fillna(method="bfill")

       
        # Check that all coverage factors are 2 and distributions are 
        # normal. If not raise exception for now, but could add a bit 
        # in the future.
        for no in range(len(
            self.all_measurement_input.filter(
            regex='distribution').columns)):
            a = self.all_measurement_input[
                self.all_measurement_input.filter(
                regex='distribution').columns[no]].to_numpy()
            a = np.array([x.lower() for x in a])

            b = self.all_measurement_input[
                self.all_measurement_input.filter(
                regex='coverage').columns[no]].to_numpy()
            if not (a[0] == a).all() and a[0].lower() == 'normal':
                raise Exception('Please make sure that all distributions'
                                +' are normal.')
            if not(b[0] == b).all() and (b[0].lower() == 2):
                raise Exception('Please make sure that all coverage'
                                +' factors are k=2.') 
         
    def draw_lhs_monte_carlo_samples(
            self, n_draws=10000,frame=None, 
            temp_lim=None, pres_lim=None, quant_lim=None): 
        """Subsets, then operates on met4H2.all_measurement_input to 
        create monte carlo draws for pressure, temperature, quantity and 
        composition by default. Compostion is normalized to 1
        
        The monte carlo draws are saved in a dictionary which can be
        indexed on parameter name: self.monte_carlo_draws.

        If the optional argument frame is set (it must be a data frame)
        then the function will not do default draws for pressure,
        temperature, quantity and composition, but rather create monte
        carlo samples for the inputted frame. One might want to do this
        to create extra noise for instance.
        
        The function stores the results in the object in data frames 
        contained in a dictionary, where columns are MC draws and 
        time points are in rows. 
        
        temp_lim, pres_lim and quant_lim are optional parameters, 
        that if set should be lists where the first entry is the lower 
        acceptable limit of the corresponding parameter, and the second 
        entry is the upper acceptable limit."""


        # Pressure, temperature and quantity do not need normalizing, 
        # that is why the draws are done in a loop,
        # on the first iteration, there is no normalizing, on the second 
        # there is.

        # This is the default way to run this function. It will create
        # monte carlo draws for pressure, temperature, quantity and 
        # composition.
        if frame is None:
            return_out = False
            # Subsetting the big input frame based on column names.
            # Note that these names are set earlier in code, when
            # the input was read, so they are well-defined.

            # Because the composition has to be normalized, make
            # two separate frames, one for pressure, temperature and
            # quantity, and one for composition.    
            unit = self.all_measurement_input['quantity_unc_unit'].iloc[0]
            pressure_temp_quantity_frame = self.all_measurement_input[
                ['quantity_'+unit,
                'quantity_unc',
                'quantity_unc_coverage',
                'quantity_unc_distribution',
                'pressure_Pa',
                'pressure_unc',
                'pressure_unc_coverage',
                'pressure_unc_distribution',
                'temperature_K',
                'temperature_unc',
                'temperature_unc_coverage',
                'temperature_unc_distribution']]
            pressure_temp_quantity_frame=pressure_temp_quantity_frame.rename(
                columns={'quantity_'+unit:'quantity',
                         'pressure_Pa':'pressure',
                         'temperature_K':'temperature'})
            # Now update uncertainty to standard uncertainty. 
            # For now do normal distribution, could easily be added to 
            # in the future.
            for par in ['quantity_unc','pressure_unc','temperature_unc']:
                mask = pressure_temp_quantity_frame[
                    par+'_distribution'].str.startswith('Normal')
                pressure_temp_quantity_frame.loc[
                    mask,par] = pressure_temp_quantity_frame[
                    par]/pressure_temp_quantity_frame[
                    par+'_coverage']
                pressure_temp_quantity_frame.loc[
                    mask,par+'_coverage']  = 1
                
            
            composition_frame = self.all_measurement_input.filter(
                regex='^((?!quantity).)*$').filter(
                regex='^((?!pressure).)*$').filter(
                regex='^((?!temperature).)*$').filter(
                regex='^((?!ptime).)*$').filter(
                regex='^((?!composition).)*$').filter(
                regex='^((?!DateTime).)*$')
            # Now update uncertainty to standard uncertainty. 
            # For now do normal distribution, could easily be added to 
            # in the future.
            comp_unc = composition_frame.filter(regex='_unc').columns
            for par in comp_unc: 
                composition_frame[par][
                    (self.all_measurement_input[
                    'composition_unc_distribution'
                    ] == 'Normal')] = composition_frame[
                    par]/self.all_measurement_input[
                    'composition_unc_coverage']
            frame_list = [pressure_temp_quantity_frame,
                          composition_frame]
        
        # If not running the function with default monte carlo
        # draws for pressure, temperature, quantity and 
        # composition, but rather providing a self defined frame:
        else:
            # Return results if not running function in default mode.
            return_out = True
            frame_list = [frame]

        # frame_list will have length 1 or 2, depending in wether
        # method is run in default mode or not (see above comments).   
        for i in range(len(frame_list)):
            # Get the frame. It contains parameters and their 
            # uncertainties.
            frame = frame_list[i]
            out_frames = list()
            # set number of monte carlo draws to be made.
            mc.npts = n_draws
            # The out_frame will have draws in columns, so number them.
            col_names = list(range(1, n_draws + 1))
            col_names = [str(x) for x in col_names]
            parameter_list = [x for x in frame.columns if not 'unc' in x]
            # A seperate out_frame with monte_carlo draws in columns
            # and potential time points in rows will be made for each
            # parameter in frame.
            for parameter in parameter_list:
                # Extract the parameter and its uncertainty from frame.
                val = pd.concat(
                    [frame[parameter],frame[parameter+'_unc']],axis=1)
                # Create data frame to use for storing the monte carlo 
                # draws for each parameter (eg., pressure, temp, H2,...).
                mc_df = pd.DataFrame(columns=col_names)
                # Create a new draw for each row
                for _, row in val.iterrows():
                    # make draws
                    generic_draws = mc.N(0,1) 
                    # extract draws
                    dist = generic_draws._mcpts
                    # redistribute the draws within 
                    # the uncertainty limits.
                    new_dist = pd.DataFrame(
                        row[parameter] 
                        + (dist - dist.mean())
                        *(row[parameter+'_unc']/dist.std())).T
                    new_dist.columns = col_names
                    # Add the monte carlo draw in row to the mc_df.
                    mc_df = pd.concat([mc_df,new_dist],axis=0)
                # reindex and add the data frame for the particular
                # parameter to the list of out_frames.
                mc_df.index=list(range(len(val)))
                out_frames.append(mc_df)

            # If run in default mode, the second of two possible
            # iterations will cause us to make draws for the
            # composition. This we want to normalize.
            if i == 1:
                sum_across_frames = sum(out_frames)
                for no in range(len(out_frames)):
                    # Set limit between 0 and 1.
                    out_frames[no][out_frames[no] < 0] = 0
                    out_frames[no][out_frames[no] > 1] = 1
                    self.monte_carlo_draws[
                        parameter_list[no]] = (
                        out_frames[no]/sum_across_frames)
                    self.monte_carlo_draws[
                        parameter_list[no]] = out_frames[
                        no].add_prefix(parameter_list[no]+'_') 
            # If the draws are not for composition, which always
            # is the case if i=0 (i can maximally be 1), truncate
            # the draw distribution if limits are set in optional
            # arguments.
            else:
                limits = [quant_lim,pres_lim,temp_lim]
                for no in range(len(out_frames)):
                    # Set lim between those potentially given 
                    # with temp_lim, pres_lim and quant_lim
                    if isinstance(limits[no],list):
                        out_frames[no][
                            out_frames[no]<limits[no][0]]=limits[no][0]
                        out_frames[no][
                            out_frames[no]>limits[no][1]]=limits[no][1]
                    # The out_frames for each parameter are added to 
                    # a dictionary containing the montecarlo draws, 
                    # which can later be indexe on parameter name.
                    self.monte_carlo_draws[
                        parameter_list[no]] = out_frames[no].add_prefix(
                        parameter_list[no]+'_')             
        
        # Only return if not running funtion in default mode.
        # In that case out_frames will have length=1.
        if return_out == True:
            return out_frames[0]
    
    def get_pressure(self,nominal,standard,nominal_values=None):
        """This function retrieves pressure at line conditions 
        from a data frame, and returns them along with pressure
        at reference conditions as defined in __init__.
        
        Inputs: nominal=False or nominal=True:
                indicates whether one is about to calculate nominal 
                gas properties or if one is about to calculate gas 
                properties for a set of monte carlo draws.

                nominal_values:  
                should be a data frame, containing
                nominal parameter values if nominal=True
                
                standard: 
                indicates which algorithm to use for gas property
                calculation. The standards require different 
                pressure units."""
        
        if standard == 'AGA8Detail': divide_by = 1000 # go to kPa
        if standard == 'TREND': divide_by = 1 # keep Pa
        if nominal == True:
            pressure = nominal_values.filter(
                regex='^pressure')/divide_by 
            pressure = pressure.rename(
            columns={pressure.columns[0]:'pressure_1'})
        else:
            pressure = self.monte_carlo_draws[
                'pressure'].iloc[nominal_values.index,:]/divide_by 
        P0 = self.P0/divide_by
        return np.array(P0), pressure

    def get_temperature(self,nominal, nominal_values=None):
        """This function retrieves temperature at line conditions 
        from a data frame, and returns them along with temperature
        at reference conditions as defined in __init__.
        
        Inputs: nominal=False or nominal=True:
                indicates whether one is about to calculate nominal 
                gas properties or if one is about to calculate gas 
                properties for a set of monte carlo draws.

                nominal_values:  
                should be a data frame, containing
                nominal parameter values if nominal=True"""
        
        # If calculating nominal gas properties, get the temperature
        # values straight from the data frame with nominal values.
        if nominal == True: 
            temperature = nominal_values.filter(regex='^temperature')
            temperature = temperature.rename(
                columns={temperature.columns[0]:'temperature_1'})
        else: temperature = self.monte_carlo_draws[
                'temperature'].iloc[nominal_values.index,:]
        return np.array(self.T0), temperature

    def get_gross_calorific_value(
            self,nominal,component_names,n_draws=10000):
        """Returns H0_frame, which is a data frame containint gross
        calorific values (units: J/mol) at reference temperature. 
        Reference temperature self.gcv_temperature set in __init__.
        
        Inputs: nominal:
                Boolean, True or False, indicates whether about to
                calculate nominal gas properties or calculating 
                monte carlo gas properties.
                
                component_names: list of gas component names
                
                n_draws: number of monte carlo draws."""
        
        # To not have to read extra files the gross calorific value 
        # at different temperatures is hardcoded.
        gcv_df = pd.DataFrame({"Component":[
            'Methane','Ethane','Propane','n-Butane','2-Methylpropane',
            'n-Pentane','2-Methylbutane','2,2-Dimethylpropane','n-Hexane',
            '2-Methylpentane','3-Methylpentane','2,2-Dimethylbutane',
            '2,3-Dimethylbutane','n-Heptane','n-Octane','n-Nonane',
            'n-Decane','n-Undecane','n-Dodecane','n-Tridecane',
            'n-Tetradecane','n-Pentadecane','Ethene','Propene','1-Butene',
            'cis-2-Butene','trans-2-Butene','2-Methylpropene','1-Pentene',
            'Propadiene','1,2-Butadiene','1,3-Butadiene','Ethyne',
            'Cyclopentane','Methylcyclopentane','Ethylcyclopentane',
            'Cyclohexane','Methylcyclohexane','Ethylcyclohexane',
            'Benzene','Toluene','Ethylbenzene','o-Xylene','Methanol',
            'Methanetiol','Hydrogen','Water','Hydrogen sulfide','Ammonia',
            'Hydrogen cyanide','Carbon monoxide',
            'Carbonyl sulfide','Carbon disulfide'],

                   "Name":[
            'C1','C2','C3','n-C4','i-C4','n-C5','i-C5','neo-C5','n-C6',
            '2-Methylpentane','3-Methylpentane','2,2-Dimethylbutane',
            '2,3-Dimethylbutane','n-C7','n-C8','n-C9','n-C10','n-C11',
            'n-C12','n-C13','n-C14','n-C15','Ethylene','Propene',
            '1-Butene','cis-2-Butene','trans-2-Butene','2-Methylpropene',
            '1-Pentene','Propadiene','1,2-Butadiene','1,3-Butadiene',
            'Acetylene','Cyclopentane','Methylcyclopentane',
            'Ethylcyclopentane','Cyclohexane','Methylcyclohexane',
            'Ethylcyclohexane','Benzene','Toluene','Ethylbenzene',
            'o-Xylene','Methanol','Methanetiol','H2','H2O',
            'H2S','NH3','HCN','CO','COS','CS2'],

                   "molwt(g/mol)":[
            16.04246,30.06904,44.09562,58.12220,58.12220,72.14878,
            72.14878,72.14878,86.17536,86.17536,86.17536,86.17536,
            86.17536,100.20194,114.22852,128.25510,142.28168,156.30826,
            170.33484,184.36142,198.38800,212.41458,28.05316,42.07974,
            56.10632,56.10632,56.10632,56.10632,70.13290,40.06386,
            54.09044,54.09044,26.03728,70.13290,84.15948,98.18606,
            84.15948,98.18606,112.21264,78.11184,92.13842,106.16500,
            106.16500,32.04186,48.10746,2.01588,18.01528,34.08088,
            17.03052,27.02534,28.01010,60.07510,76.14070],

                   "Hcmol_0":[
            892.92,1564.35,2224.03,2883.35,2874.21,3542.91,3536.01,
            3521.75,4203.24,4195.64,4198.27,4185.86,4193.68,4862.88,
            5522.41,6182.92,6842.69,7502.22,8162.43,8821.88,9481.71,
            10141.65,1413.55,2061.57,2721.57,2714.88,2711.09,2704.88,
            3381.32,1945.26,2597.15,2544.14,1301.86,3326.14,3977.05,
            4637.2,3960.68,4609.33,5272.76,3305.12,3952.77,4613.16,
            4602.18,766.6,1241.64,286.64,45.064,562.93,384.57,671.92,
            282.8,548.01,1104.05],

                   "Hcmol_15":[
            891.51,1562.14,2221.1,2879.76,2870.58,3538.6,3531.68,3517.44,
            4198.24,4190.62,4193.22,4180.83,4188.61,4857.18,5516.01,
            6175.82,6834.9,7493.73,8153.24,8811.99,9471.12,10130.23,
            1412.12,2059.43,2718.71,2711.94,2708.26,2702.06,3377.76,
            1943.97,2595.12,2542.11,1301.37,3222.19,3972.46,4631.93,
            3956.02,4604.08,5266.9,3302.9,3949.83,4609.54,4598.64,
            765.09,1240.28,286.15,44.431,562.38,383.51,671.67,282.91,
            548.14,1104.32],

                   "Hcmol_15_55":[
            891.46,1562.06,2220.99,2879.63,2870.45,3538.45,3531.52,
            3517.28,4198.06,4190.44,4193.04,4180.65,4188.43,4856.98,
            5515.78,6175.56,6834.62,7493.42,8152.91,8811.63,9470.73,
            10129.82,1412.07,2059.35,2718.6,2711.83,2708.16,2701.96,
            3377.63,1943.92,2595.05,2542.03,1301.35,3222.05,3972.29,
            4631.74,3955.85,4603.89,5266.69,3302.81,3949.72,4609.4,
            4598.52,765.03,1240.23,286.13,44.408,562.36,383.47,671.66,
            282.91,548.15,1104.33],

                   "Hcmol_20":[
            891.05,1561.42,2220.13,2878.58,2869.39,3537.19,3530.25,
            3516.02,4196.6,4188.97,4191.56,4179.17,4186.94,4855.31,5513.9,
            6173.48,6832.33,7490.93,8150.21,8808.73,9467.63,10126.52,
            1411.65,2058.73,2717.76,2710.97,2707.33,2701.13,3376.59,
            1943.54,2594.46,2541.44,1301.21,3320.89,3970.95,4630.2,
            3954.49,4602.36,5264.97,3302.16,3948.86,4608.34,4597.48,
            764.59,1239.84,285.99,44.222,562.19,383.16,671.58,282.95,
            548.19,1104.4],

                   "Hcmol_25":[
            890.58,1560.69,2219.17,2877.4,2868.2,3535.77,3528.83,3514.61,
            4194.95,4187.32,4189.9,4177.52,4185.28,4853.43,5511.8,6171.15,
            6829.77,7488.14,8147.19,8805.48,9464.15,10122.82,1411.18,
            2058.02,2716.82,2710,2706.4,2700.2,3375.42,1943.11,2593.79,
            2540.77,1301.05,3319.59,3969.44,4628.47,3952.96,4600.64,
            5263.05,3301.43,3947.89,4607.15,4596.31,764.09,1239.39,
            285.83,44.013,562.01,382.81,671.5,282.98,548.23,1104.49],

                   "uHcmol":[
            0.19,0.51,0.51,0.72,0.72,0.23,0.23,0.25,0.32,0.53,0.53,0.48,
            0.46,0.67,0.76,0.81,0.87,1.54,1.13,1.21,1.32,1.44,0.21,0.34,
            0.39,0.5,0.47,0.42,0.73,0.6,0.4,0.41,0.32,0.36,0.56,0.71,
            0.32,0.71,0.95,0.27,0.51,0.66,0.76,0.13,0.32,0.02,0.004,0.23,
            0.18,1.26,0.06,0.24,0.43]}) # This bit are standard uncs.

        df = gcv_df.set_index("Name").reindex(component_names).fillna(0) 
    
        if nominal == False:
            H0_df = df[
                ["Hcmol_"+str(self.gcv_temperature),"uHcmol"]].rename(
                columns={
                "uHcmol":"Hcmol_"+str(self.gcv_temperature)+"_unc"}) 
            # This is a frame that's a bit different from the others 
            # frames that can be supplied by the same method. 
            # It has components, not time along the index.
            H0_frame = self.draw_lhs_monte_carlo_samples(
                n_draws=n_draws,frame=H0_df) 
            H0_frame = H0_frame.rename(
                index=dict(zip(list(range(
                len(H0_frame))),component_names))).T 
        else: 
            H0_frame = df[
                ["Hcmol_"+str(self.gcv_temperature)]].T.set_index(
                pd.Index([1]))
            H0_frame.columns.name = None
        
        # Multipy by 1000 to get the gross calorific value in J/mol.
        return H0_frame*1000
    
    def call_AGA8Detail(
            self, T0, P0, temperature, pressure, 
            component_names, component_values):
        """This function returns a list containing the compressibility
        factor at line conditions (Z), the compressibility factor (Z0)
        at reference condtitions, total molecular weight (total_mw), 
        molar density (DM), molar density at reference conditions(DM0),
        and gross calorific value (Hmol), in that order. 
        
        Inputs: T0:             temperature at reference conditions.
                P0:             pressure at reference conditions.
                temperature:    temperature at line conditions.
                pressure:       pressure at line conditions.
                component_names:list of component names.
                component values: numpy array of component values."""

        AGA8_detail_order = [
            'NA','C1','N2','CO2','C2','C3','i-C4','n-C4','i-C5','n-C5',
            'n-C6','n-C7','n-C8','n-C9','n-C10','H2','O2','CO','H2O',
            'H2S','He','Ar'] 

        draws = component_values.shape[1]
        if draws > 1: nominal = False
        else: nominal = True

        
        # Initialize data frames prior to calculation. 
        if nominal == False: frame_size = draws+1
        else: frame_size = 2
        frame_list = []
        out_pars = ['Z_','Z0_','MM_','DM_','DM0_','Hmol_']
        for a in range(len(out_pars)):
            cols = []
            for i in list(range(1, frame_size)):
                cols.append(out_pars[a]+str(i))
            frame_list.append(pd.DataFrame(
                index=np.arange(component_values.shape[0]), 
                columns=cols))
            
        # Unpack empty frames.
        (Z_frame, 
         Z0_frame, 
         MM_frame, 
         DM_frame, 
         DM0_frame, 
         Hmol_frame) = frame_list

        # Get gross calorific values.
        H0_frame = self.get_gross_calorific_value(
            nominal,component_names,n_draws=draws)

        
        for t in range(component_values.shape[0]):
            for draw in range(draws):
        
                # create a data frame (x_array) with the gas composition.
                stacked_list = [
                    [i] for i in list(component_values[t,draw,:])]
                x_array = pd.DataFrame(
                    dict(zip(component_names,stacked_list)))
                    
                # Drop columns that are not contained in AGA8.
                x_array = x_array[x_array.columns.intersection(
                    AGA8_detail_order)]
                # Create columns that contain zeros for components 
                # in AGA which are not in H0_frame. 
                for component in AGA8_detail_order:
                    if not component in x_array.columns:
                        x_array[component] = 0

                # Re-order columns in H0_frame so they match the order 
                # of components in AGA8Detail.
                x_array = x_array[AGA8_detail_order] 
                # Then do AGA8Detail calculation
                Z_frame.iloc[t,draw] = AGA8Detail(
                    P=np.array(pressure.iloc[t,draw]), 
                    T = np.array(temperature.iloc[t,draw]), 
                    x=x_array.iloc[0,:].to_numpy()).run().Z
                Z0_frame.iloc[t,draw] = AGA8Detail(
                    P = P0, 
                    T = T0, 
                    x = x_array.iloc[0,:].to_numpy()).run().Z 
                MM_frame.iloc[t,draw] = AGA8Detail(
                    P=P0, 
                    T = T0, 
                    x = x_array.iloc[0,:].to_numpy()).run().MM
                DM_frame.iloc[t,draw] = AGA8Detail(
                    P=np.array(pressure.iloc[t,draw]), 
                    T = np.array(temperature.iloc[t,draw]), 
                    x=x_array.iloc[0,:].to_numpy()).run().D
                DM0_frame.iloc[t,draw] = AGA8Detail(
                    P=P0, 
                    T = T0, 
                    x = x_array.iloc[0,:].to_numpy()).run().D 

                # Get the columns that are in both x_array and H0_frame.
                x_array = x_array[x_array.columns.intersection(
                    H0_frame.columns)] 
                for component in H0_frame.columns.to_list():
                    if not component in x_array.columns:
                        x_array[component] = 0
                x_array = x_array[H0_frame.columns.to_list()] 

                # Get Hmol by multiplying the concentration of each 
                # component (x_array) by gross calorific value (H0) 
                # and sum across components.
                Hmol_frame.iloc[t,draw] = (
                    H0_frame.iloc[draw,:]*x_array).sum(axis = 1)[0]
        
        # Divide by 1000 to get correct units.
        to_return = (Z_frame, Z0_frame, MM_frame/1000, DM_frame*1000, 
                     DM0_frame*1000, Hmol_frame)    
        return to_return
    
    def call_TREND(
            self,T0, P0, temperature, pressure, 
            component_names, component_values,
            eqn_of_state=None):
        """This function returns a list containing the compressibility
        factor at line conditions (Z), the compressibility factor (Z0)
        at reference condtitions, total molecular weight (total_mw), 
        molar density (DM), molar density at reference conditions(DM0)
        in that order, and gross calorific value (Hmol) in that order. 
        
        Inputs: T0:             temperature at reference conditions.
                P0:             pressure at reference conditions.
                temperature:    temperature at line conditions.
                pressure:       pressure at line conditions.
                component_names:list of component names.
                component values: numpy array of component values."""

        # Get equation of state and match to recquired input digit in
        # the TREND framework.
        if eqn_of_state == None or eqn_of_state == 'GERG2008': 
            eqn_of_state = 1
            mix_ind = 1
        elif eqn_of_state == 'AGA8': 
            eqn_of_state = 7
            mix_ind = 7
        else: raise Exception(
                'Please make sure to use a recognized equation of state.')
        
        
             
        # Determine wether gas property computation is performed
        # for nominal values, or whether the input is a data frame
        # with monte carlo draws.
        n_draws = component_values.shape[1]
        if n_draws > 1: nominal = False
        else: nominal = True

        # Initialize data frames prior to calculation.
        if nominal == False: frame_size = n_draws+1
        else: frame_size = 2
        frame_list = []
        out_pars = ['Z_','Z0_','MM_','DM_','DM0_','Hmol_']
        for a in range(len(out_pars)):
            frame_list.append(pd.DataFrame(
                index=np.arange(component_values.shape[0]), 
                columns=[out_pars[a]+str(i) for i in list(
                range(1,frame_size))]))
        

        # Get gross calorific values.
        H0_frame = self.get_gross_calorific_value(
            nominal,component_names,n_draws=n_draws)
        
        # REMOVE LATER - FOR TESTING PURPOSES
        #draws = 1

        # Some of the algorithms that will be used later only read 
        # the component names in specific formats. Therefore, creating
        # lists of component names.
        coolProp_component_names = [
            self.component_dictionary[i][0] for i in component_names]
        TREND_component_names = [
            self.component_dictionary[i][1] for i in component_names]
        
        # For each draw use the tr object to calculate gas properties.
        for draw in range(n_draws): 
            component_values_draw = component_values[:,draw,:]

           # REMOVE LATER - FOR TESTING PURPOSES.
            #component_values_draw = np.vstack(([0.6,0.4], [0.6,0.4], [0.6,0.4]))
            #TREND_component_names = ['Methane','Ethane']

            # Create the fluid input dictionary.
            fluid_input={
                'input' : 'TP',   # temperature and pressure
                'calctype' : 'D', # default density molar
                'fluids' : TREND_component_names, 
                'moles' : component_values_draw, # molfraction, sum=1 
                'eos_ind' : len(TREND_component_names)*[eqn_of_state], 
                'mix_ind' : mix_ind,
                'path': self.trend_path,
                'unit': 'molar', # unit
                'dll_path' : self.trend_dll_path}

            # set object
            tr = trend(fluid_input)
            # set temperatures and pressures
            tr.cp_fluids = coolProp_component_names
            tr.P0 = P0 #Pa
            tr.T0 = T0 #K
            press = pressure[pressure.columns[draw]].to_numpy()
            temp = temperature[temperature.columns[draw]].to_numpy()

            # REMOVE LATER - FOR TESTING PURPOSES
            #press=np.array([101325,10000000,10000000]) #Pa     
            #temp=np.array([273.15+25,273.15+60,273.15+60]) #K
            #tr.cp_fluids = ['Methane','Ethane']
            
            # Calculate properties.
            tr.calc_properties(temp,press)

            # REMOVE LATER - FOR TESTING PURPOSES.
            #H0_frame = pd.DataFrame(H0_frame.iloc[0,[3,4]]).T
            #component_values_draw = np.vstack(([0.6,0.4], [0.6,0.4], [0.6,0.4]))

            # Get Hmol by multiplying the concentration of each 
            # component (component_values) by gross calorific value (H0) 
            # and sum across components.
            for t in range(len(component_values_draw)):
                frame_list[-1].iloc[t,draw] = (
                    component_values_draw[t,:]*H0_frame.values).sum(
                    axis=1)[0]

            # Put non-TREND-specific results into return frames.
            trend_results = [tr.Z,tr.Z0,tr.total_mw,tr.DM,tr.DM0] 
            for n in range(len(frame_list[:-1])):
                frame = frame_list[:-1][n]
                frame[frame.columns[draw]] = trend_results[n]

        return frame_list
    
    def compute_gas_properties(
            self,method='TREND',eqn_of_state=None,nominal=False,):
        """The main output of this function are compressibility
        factors and calorific values, however, other gas properties,
        such as densities and molar densitites are also computed. 

        method = 'TREND', or method = 'AGA8Detail', this refers
        to the algorithm that the function will call to perform the 
        calculations of gas properties.

        If method = 'TREND' one can choose to set eqn_of_state to
        'AGA8' or 'GERG2008'. If method is not set, it will be 
        'GERG2008' if method = 'TREND' and 'AGA8' if method ='AGA8Detail'.
        
        If nominal=False, gas properties are calculated using 
        monte carlo draws and output is saved to self.main_m4h2_output.

        If nominal=True, monte carlo is not used and the output is saved
        to self.main_nominal_m4h2_output. 

        Note that if nominal=True computed in the latter case, 
        uncertainties will not be computed and model uncertainties 
        will not be added to the nominal values."""

        if method not in ['AGA8Detail','TREND']: 
            raise Exception(
                'Please make sure to apply a supported method for'
                 +'computing gass properties.')
        if nominal == False: 
            # Extracting the number of monte_carlo_draws that have
            # been performed, 'temperature' is one of the paramters that
            # must have been drawn prior to gas property calculation.
            n_draws = self.monte_carlo_draws['temperature'].shape[1]
        else: 
            n_draws = 1
        

        # get a data frame with nominal values. 
        # Column names were defined when saving the input to objects, 
        # so these are well known.
        nominal_values = self.all_measurement_input.filter(
                regex='^((?!ptime).)*$').filter(
                regex='^((?!composition).)*$').filter(
                regex='^((?!DateTime).)*$').filter(
                regex='^((?!unc).)*$')
        
        
        # get gas component names.
        component_names = []
        for comp in nominal_values.columns.to_list():
            if not any(
                non_comp in comp for non_comp in [
                'quantity','pressure','temperature']):
                component_names.append(comp)
        

        # Find all duplicate rows in nominal values and create a 
        # dictionary with keys that keep the first index, and values 
        # that are lists of duplicate rows. This will improve speed
        # because we reduce number of calculations.
        kept_nominal_values = nominal_values.drop_duplicates(
            keep='first')
        keys = nominal_values.drop_duplicates(
            keep='first').index.to_list()
        values = nominal_values[
            nominal_values.duplicated()].index.to_list()
        duplicates = {}
        for keyno in range(len(keys)):
            duplicates[keys[keyno]] = []
            for idx in values: 
                if nominal_values.iloc[idx,:].equals(
                    nominal_values.iloc[keys[keyno],:]):
                    duplicates[keys[keyno]].append(idx)

        
        # Concatenate data frames to 3D numpy arrays for easy 
        # collection of desired arrays of the composition.
        composition_frames = []
        for comp in component_names:
            if nominal == False: 
                composition_frames.append(self.monte_carlo_draws[comp])
                stacked = np.zeros([len(nominal_values),n_draws])
            else: 
                composition_frames.append(
                    pd.DataFrame(nominal_values[comp]))
                stacked = np.zeros([len(nominal_values),1])
        for comp_no in range(len(component_names)):
            stacked = np.dstack(
                (stacked,composition_frames[comp_no].to_numpy()))
        stacked = stacked[:,:,1:]
        stacked = stacked[list(duplicates.keys()),:,:]

        # Get temperature and pressure.
        T0, temperature = self.get_temperature(
            nominal,nominal_values=kept_nominal_values)
        P0, pressure = self.get_pressure(
            nominal,method,nominal_values=kept_nominal_values)
        # Get components of composition.
        if nominal == True: 
            reshaped = kept_nominal_values[
                component_names].to_numpy().reshape(
                (kept_nominal_values[component_names].to_numpy().shape[0], 
                 kept_nominal_values[component_names].to_numpy().shape[1], 1))
            component_values = np.moveaxis(reshaped, [2,1], [1,2])
        else: component_values = stacked

        # Calculations typically doesn't work with neopentane in the mix. 
        # According to ISO 20765-2, neo-pentane can be substituted
        # by n-pentane. The uncertainty will be dominated by n-pentane.
        if 'neo-C5' in component_names:
            component_values[
                :,:,component_names.index('n-C5')] = component_values[
                :,:,component_names.index('neo-C5')]+component_values[
                :,:,component_names.index('n-C5')]
            component_values = np.delete(
                component_values, component_names.index('neo-C5'),axis=2)
            del component_names[component_names.index('neo-C5')]

        # Normalize, as the distribution could have been truncated after
        # monte carlo draws.
        for t in range(component_values.shape[0]):
            for draw in range(component_values.shape[1]):
                component_values[t,draw,:] = component_values[
                    t,draw,:]/component_values[t,draw,:].sum()

        # To check time efficiency of the computation below.
        start = time.time()

        # Perform calculation of gas properties according to
        # chosen standard.
        if method == 'TREND': 
            (Z_frame, 
            Z0_frame, 
            MM_frame, 
            DM_frame, 
            DM0_frame, 
            Hmol_frame) = self.call_TREND(
                T0, P0, temperature, pressure,
                component_names, component_values,
                eqn_of_state=eqn_of_state)
        if method == 'AGA8Detail':
            (Z_frame, 
            Z0_frame, 
            MM_frame, 
            DM_frame, 
            DM0_frame, 
            Hmol_frame) = self.call_AGA8Detail(
                T0, P0, temperature, pressure, 
                component_names, component_values)

    
        # outputs common to AGA8Detail and TREND are dealt with here.

        # Add model noise to parameters 
        # (noise level determined in __init__(self))
        noisy_list = [
            self.model_Z_noise, 
            self.model_Z0_noise, 
            self.model_MM_noise, 
            self.model_DM_noise, 
            self.model_DM0_noise,
            self.model_Hmol_noise]
        frame_list = [
            Z_frame, 
            Z0_frame, 
            MM_frame, 
            DM_frame, 
            DM0_frame, 
            Hmol_frame]
        for i in range(len(noisy_list)):
            if nominal == False: 
                if noisy_list[i] > 0:
                    noise_unc = []
                    for n in frame_list[i].mean(axis=1).to_list():
                        noise_unc.append(n*noisy_list[i])
                    noise_frame = pd.DataFrame(
                        {'noise':len(frame_list[i])*[0],
                         'noise_unc':noise_unc})
                    noise = self.draw_lhs_monte_carlo_samples(
                        n_draws=n_draws,frame=noise_frame)
                else: 
                    noise = pd.DataFrame(
                        0, 
                        index=range(len(frame_list[i])), 
                        columns=list(range(1,n_draws+1)))
                    noise = noise.astype('object')
                # Keep only numbers in column number to be able to 
                # use the dataframes like we want to below.
                counter = 0
                colnums = []
                for x in frame_list[i].columns.to_list():
                    colnums.append(x.split('_')[-1])
                    counter = counter + 1
                frame_list[i].columns = colnums
                # Adding noise (it is either made in 
                # the above if, or else zero.)
                frame_list[i] = pd.DataFrame(
                    frame_list[i].values+noise.values) 
            else: frame_list[i].columns = ['0']
        # Unpacking frame list.
        _, _, MM_df, DM_df, DM0_df, Hmol_df = frame_list


        # Calculate densities and calorific values.
        D_df=pd.DataFrame(DM_df.values*MM_df.values)   # D = DM*MM      
                                                       #density kg/m3 
        D0_df=pd.DataFrame(DM0_df.values*MM_df.values) # D0 = DM0*MM    
                                                       #density kg/m3
        Hv_df=pd.DataFrame(Hmol_df.values/(MM_df.values/D0_df.values))  
                                            # Hv = (Hmol/(MM/D0)) 
                                            # calorific value J/m3
        Hm_df=pd.DataFrame(Hmol_df.values/MM_df.values)
                                            # Hm = Hmol/MM
                                            # calorific value J/m3

        # Now fillna according to duplicated input and add prefixes 
        # to be able to understand what's in the data frames later. 
        calculated = frame_list+[D_df,D0_df,Hv_df,Hm_df]
        frame_names = [
            'Z_','Z0_','MM_','DM_','DM0_','Hg_','D_','D0_','Hv_','Hm_']
        for n in range(len(frame_names)):
            df = pd.DataFrame(
                columns=calculated[n].columns, 
                index=list(
                range(
                max(
                [max(list(duplicates.keys())),
                 max(max(duplicates.values()))])+1))).add_prefix(
                frame_names[n])
            for keyno in range(len(duplicates)):
                df.iloc[list(
                    duplicates.keys())[keyno],:] = calculated[n].iloc[
                    keyno,:]
                for idx in list(duplicates.values())[keyno]:
                    df.iloc[idx,:] = calculated[n].iloc[keyno,:]
            calculated[n] = df

        # To check time efficiency of calculation above.
        elapsed_time = (time.time() - start) 

        # Put output in dictionary and save to object.
        if nominal == False: 
            self.main_m4h2_output = dict(zip(
                ['Z','Z0','MM','DM','DM0','Hg','D','D0','Hv','Hm'],
                calculated))
            return self.main_m4h2_output
        else: 
            self.main_nominal_m4h2_output = dict(zip(
                ['nominal_Z',
                'nominal_Z0',
                'nominal_MM',
                'nominal_DM',
                'nominal_DM0',
                'nominal_Hg',
                'nominal_D',
                'nominal_D0',
                'nominal_Hv',
                'nominal_Hm'],calculated))
            return self.main_nominal_m4h2_output

    def initialize_uncertainty_calculator(
            self,ucalc,gas_properties,
            computation_method='lhs_montecarlo',n_mc_draws=1000):
        
        """uncertaintyCalculator has a number of methods that help 
        storing the information it needs to perform the uncertainty
        calculation to the uncertaintyCalculator object. 
        
        This initialize_uncertainty_calculator passes information from
        the met4H2 object to the uncertainty calculator for energy and
        energy flow calculations."""

        # Get what sort of flow we have. This will affect
        # the choice of functional relationship. 
        unitstring = self.all_measurement_input.filter(
            regex='^((?!unc).)*$').filter(
            regex='^quantity').columns.str.split('_')[0][-1]
        # Settings file is quite restricted in terms of what inputs
        # are allowed.
        possible_units = {
            'm³/s':'volume_flow',
            'Sm³/s':'standard_vol_flow',
            'kg/s':'mass_flow','m³':'volume',
            'Sm³':'standard_vol','kg':'mass'}
        flowtype = possible_units[unitstring]

                                            
        if flowtype == 'volume_flow':                                                                                            
            output_label_and_unit = ['qE|J/s']  
            # Define input variables and units for each input variable.
            labels_and_units = [
                'P0|Pa',
                'T0|K',
                'P|Pa',
                'T|K',
                'quantity|m³/s',
                'Z|-',
                'Z0|-',
                'Hv|J/m³'] 
            # Define functional relationship.
            functional_relationship = '(P*T0*Z0*Hv*quantity)/(P0*T*Z)'         # Define functional relationship                                                        # Define output variable and corresponding unit
        if flowtype == 'standard_vol':
            output_label_and_unit = ['E|J'] 
            labels_and_units = ['P0|Pa',
                                'T0|K',
                                'P|Pa',
                                'T|K',
                                'quantity|Sm³',
                                'Z|-',
                                'Z0|-',
                                'Hv|J/m³']     
            # Need to convert standard volume to volume at line 
            # conditions in the functional relationship.
            # Consider if it is better to run calculator twice in this
            # case, since it affects the sensititivy coefficients.
            fr='(P*T0*Z0*Hv*((P0*T*Z*quantity)/(P*T0*Z0)))/(P0*T*Z)'
            functional_relationship=fr        
        if flowtype != 'volume_flow' and flowtype != 'standard_vol':
            raise Exception('Code isn\'t finished, please use volume'
                            +' flow as quantity input.')       

                                                     
        
        # Define values for each input variable and time point
        P0 = [float(self.P0)]*len(self.all_measurement_input)
        T0 = [float(self.T0)]*len(self.all_measurement_input)
        P = self.all_measurement_input.filter(
            regex='^((?!unc).)*$').filter(
            regex='^pressure').iloc[:,0].to_list()
        T = self.all_measurement_input.filter(
            regex='^((?!unc).)*$').filter(
            regex='^temperature').iloc[:,0].to_list()
        quantity = self.all_measurement_input.filter(
            regex='^((?!unc).)*$').filter(
            regex='^quantity').iloc[:,0].to_list()
        Z = gas_properties['Z'].mean(axis=1).to_list()
        Z0 = gas_properties['Z0'].mean(axis=1).to_list()
        Hv = gas_properties['Hv'].mean(axis=1).to_list()
        values = [P0,T0,P,T,quantity,Z,Z0,Hv]      

        # Define timestamps for each time point
        time_stamps = self.all_measurement_input['DateTime']

        # Set input and output quantities and functional relationship in 
        # uncertaintyCalculator. It is important to define the 
        # input quantities first.
        ucalc.set_input_quantities(
            labels_and_units,values,time_stamps=time_stamps)
        ucalc.set_output_labels_and_units(output_label_and_unit) 
        ucalc.set_functional_relationship(functional_relationship) 

        # Define uncertainty distributions
        P0_u = [0]*len(self.all_measurement_input)
        T0_u = [0]*len(self.all_measurement_input)
        P_u = self.all_measurement_input['pressure_unc'].to_list()
        T_u = self.all_measurement_input['temperature_unc'].to_list()
        quantity_u = self.all_measurement_input['quantity_unc'].to_list()
        Z_u = gas_properties['Z'].std(axis=1).to_list()
        Z0_u = gas_properties['Z0'].std(axis=1).to_list()
        Hv_u = gas_properties['Hv'].std(axis=1).to_list()
        input_expanded_uncertainties = [
            P0_u,T0_u,P_u,T_u,quantity_u,Z_u,Z0_u,Hv_u]  

        distributions =  ['normal']*len(labels_and_units)

        # Set uncertainty distributions in uncertainty calculator
        ucalc.set_input_expanded_uncertainties(
            input_expanded_uncertainties,
            absolute=True,
            distributions=distributions)

        # Define correlation matrix in uncertainty calculator
        forbidden_correlations=[
            ['P','T'],['P0','T0'],['P0','T'],['P0','P'],
            ['P0','Hv'],['P0','Z0'],['P0','Z'],['T0','T'],
            ['T0','P'],['T0','Hv'],['T0','Z0'],['T0','Z'],
            ['quantity','T0'],['quantity','T'],['quantity','P0'],
            ['quantity','Hv'],['quantity','Z0'],['quantity','Z'],
            ['pressure','temperature'],['P0','temperature'],
            ['P0','pressure'],['T0','temperature'],
            ['T0','pressure'],['quantity','temperature']]
        P0_frame = pd.DataFrame(
            float(self.P0),
            index = np.arange(len(self.all_measurement_input)), 
            columns=[
            'P0_'+str(i) for i in range(len(gas_properties['Z'].columns))
            ])
        T0_frame = pd.DataFrame(
            float(self.T0),
            index = np.arange(len(self.all_measurement_input)), 
            columns=[
            'T0_'+str(i) for i in range(len(gas_properties['Z'].columns))
            ])
        frame_list = [
            P0_frame,
            T0_frame,
            self.monte_carlo_draws['pressure'],
            self.monte_carlo_draws['temperature'],
            self.monte_carlo_draws['quantity'],
            gas_properties['Z'],
            gas_properties['Z0'],
            gas_properties['Hv']] 
        (ucalc.input_dictionary['list_correlation_matrix'], 
         _) = self.make_input_correlation_matrices(
            frame_list ,forbidden_correlations)

        # Set computational method in uncertainty calculator
        ucalc.set_computational_method(
            computation_method=computation_method,
            monte_carlo_draws=n_mc_draws)

        # Perform SI conversion and a few other partial calculations 
        # and add these to the object.
        # It is not necessary to do SI_conversion, as we have already 
        # done that when we set the stream and meter objects.
        ucalc.initialize_calculator(SI_conversion=False)  
        
    def make_input_correlation_matrices(
            self, frame_list, forbidden_correlations=None):
        """There will be correlations between some of the gas properties,
        since some of them depend on others.
        
        This method finds correlaitons empirically, then forces zero 
        correlation between properties that should not be correlated.
        Outputs upper triangular matrix for each time point contained 
        in the data frames in frame_list"""

        if forbidden_correlations == None: forbidden_correlations = list()
        out_matrices = list()

        # Calculate correlation at each time point t.
        for t in range(len(frame_list[0])):
            out_labels = list()
            out_matrix = list()

            # Loop through parameters.
            for elm_no in range(len(frame_list)):
                # Create a frame list with which to correlate
                # elements in the origninal list.
                new_frame_list = frame_list[elm_no+1:]
                for new_elm_no in range(len(new_frame_list)):
                    # Add the parameters between which the
                    # correlation is calculated to out_labels.
                    lst = list([
                        frame_list[elm_no].columns[0][
                        :frame_list[
                        elm_no].columns[0].rindex('_')],
                        new_frame_list[new_elm_no].columns[0][
                        :new_frame_list[
                        new_elm_no].columns[0].rindex('_')]])
                    out_labels.append(lst)
                    # do correlation only if allowed.
                    if lst in forbidden_correlations or list(
                        reversed(lst)) in forbidden_correlations:
                        out_matrix.append(0)
                    else:
                        # Converting time row to a numpy array
                        frame_numpy = frame_list[
                            elm_no].iloc[t,:].to_numpy()
                        new_frame_numpy = new_frame_list[
                            new_elm_no].iloc[t,:].to_numpy()
                        corr = pd.Series(
                            frame_numpy).astype('float64').corr(
                            pd.Series(new_frame_numpy).astype('float64'))
                        if np.isnan(corr): 
                            corr = 0
                        out_matrix.append(corr)
            out_matrices.append(out_matrix)
        return out_matrices, out_labels

   
    # Just for QA purposes. Remove later.
    def plotter(self,dict1, dict2):
        parameters = list(dict1.keys()) + ['dict_no']
        frame = pd.DataFrame(columns=parameters)
        for parameter in parameters[:-1]:
            frame[parameter] = dict1[parameter].iloc[0,:].to_list() + dict2[parameter].iloc[0,:].to_list()
        frame['dict_no'] = len(dict1[parameter].columns)*[1] + len(dict1[parameter].columns)*[2]
        # Standardizing 
        for parameter in parameters[:-1]:
            frame[parameter] = (frame[parameter]-frame[parameter].mean())/frame[parameter].std()
        # plotting
        dd=pd.melt(frame,id_vars=['dict_no'],value_vars=parameters[:-1],var_name='gas properties')
        sns.boxplot(x='dict_no',y='value',data=dd,hue='gas properties')

        


if __name__ == '__main__':
    
    
    meter = metering()
    stream = streamProcess()
    m4h2 = met4H2()
    ucalc = uncertaintyCalculator()
    
    m4h2.trend_path = os.path.join(
        os.getcwd(),'Tools','TREND','TREND 5.0')
    m4h2.trend_dll_path = os.path.join(
        m4h2.trend_path,'TREND_x64.dll')
    filepath = os.getcwd()
    filename = 'settings_example_20230202.xlsx'
    n_draws = 100

    # Reading the settings file is one of potentially several 
    # ways to create node, stream and meter objects.
    m4h2.read_settings_file(
        ucalc, stream, meter, filename, 
        filepath = filepath, meter_type = 'USM')

    # When one is happy with node, info from stream and 
    # meter are combined in a big data frame.
    m4h2.collect_measurement_input(stream,meter)

    # Based on the m4hs.all_measurement_input frame, 
    # one can get lhs montecarlo draws.
    m4h2.draw_lhs_monte_carlo_samples(n_draws=n_draws) 


    # Decide whether to do nominal calculations with no uncertainties 
    # or not. If using TREND, default eqn of state is 'GERG2008', 
    # but this can be modified to 'AGA8'
    # These are all the ways one could call compute_gas_properties.
    AGA8_TREND_properties = m4h2.compute_gas_properties(
        method='TREND',eqn_of_state='AGA8') 
    #GERG2008_TREND_properties = m4h2.compute_gas_properties(
    #    method='TREND',eqn_of_state='GERG2008')
    AGA8_Detail_properties = m4h2.compute_gas_properties(
        method='AGA8Detail') 

    #nominal_AGA8_TREND_properties = m4h2.compute_gas_properties(
    # method='TREND',eqn_of_state='AGA8',nominal=True) 
    #nominal_GERG2008_TREND_properties = m4h2.compute_gas_properties(
    # method='TREND',eqn_of_state='GERG2008',nominal=True)
    #nominal_AGA8_Detail_properties = m4h2.compute_gas_properties(
    # method='AGA8Detail',nominal=True) 
    
    #m4h2.plotter(AGA8_TREND_properties, AGA8_Detail_properties)

    m4h2.initialize_uncertainty_calculator(
        ucalc,AGA8_TREND_properties,
        computation_method='analytical',
        n_mc_draws=1000) 
    
    ucalc.perform_uncertainty_analysis()
    ucalc.output_to_excel(os.getcwd(),'tester.xlsx')
 