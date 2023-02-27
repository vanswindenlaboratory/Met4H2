import warnings
import time

import numpy as np
import pandas as pd
import matplotlib as plt
import seaborn as sns

import mcerp as mc # scipy i stedet. salib
import CoolProp.CoolProp as cp

from modules.AGA8Detail import AGA8Detail
from modules.molecule import *
from modules.uncertaintyCalculator import uncertaintyCalculator

from modules.fluid import fluid


# trend: Todo
# - Legg inn testing, warnings og feilmeldinger


class trend:
    def __init__(self,fluid_input):
        self.fluid_input=fluid_input
        self.D=[]
        self.D0=[]
        self.Z=[]
        self.Z0=[]
        self.dielectric_constant=[]
        self.P0 = np.array(101325)
        self.T0 = np.array(273.15+15)
        self.cp_fluids = []
        
    def make_fluid(self,fluid_input):
        fld=fluid(fluid_input['input'],fluid_input['calctype'],fluid_input['fluids'],fluid_input['moles'],fluid_input['eos_ind'],fluid_input['mix_ind'],fluid_input['path'],fluid_input['unit'],fluid_input['dll_path'])
        return fld

    def calc_EOS(self,var,temperature,pressure):
        
        vals=[]
        self.fluid_input['calctype']=var
        self.fluid_copy = self.fluid_input.copy()

        for counter in range(len(self.composition)): 
            self.fluid_copy['moles']=self.composition[counter]
            self.fld=self.make_fluid(self.fluid_copy)
            if len(temperature) == len(self.composition):
                val,_ = self.fld.TREND_EOS(temperature[counter],pressure[counter])
            else: val,_ = self.fld.TREND_EOS(temperature[0],pressure[0])
            vals.append(val)
        return vals
  
    def cal_molecular_weight(self):
        fluids=self.cp_fluids
        
        mw = []
        MW=[]
        
        for _ in range(len(self.composition)):
            mwt=[]
            vals=[]
            for fluid in range(len(fluids)):
                val=cp.PropsSI('molar_mass',fluids[fluid]) #endret til å bruke CoolProp 
                mwt.append(val) 
                val=mwt[fluid]*self.fld.moles[fluid] 
                vals.append(val) 
            MW.append(sum(vals)) 
            mw.append(mwt)
        
        self.MW=MW # unit: kg/mol bare ett tidspunkt, sum elm
        self.mw=mw # unit: kg/mol individuelle elm.
        return MW, mw 

    def calc_properties(self,temperature,pressure):
        
        self.pressure=pressure
        self.temperature=temperature
        self.composition=self.fluid_input['moles']
        
        pressure = pressure/1000000 # go to MPa
        P0=[self.P0/1000000] # go to MPa
        T0=[float(self.T0)]
        

        self.DM=self.calc_EOS('D',temperature,pressure)# unit: mol/m3
        self.DM0=self.calc_EOS('D',T0,P0)         # unit: mol/m3
        self.Z=self.calc_EOS('Z',temperature,pressure) # unit: -
        self.Z0=self.calc_EOS('Z',T0,P0)               # unit: -
        self.dielectric_constant=self.calc_EOS('DE',temperature,pressure) # unit: -
        self.cal_molecular_weight()
        self.D=[self.MW[i]*self.DM[i] for i in range(len(self.D))]         # unit: kg/m3
        self.D0=[self.MW[i]*self.DM0[i] for i in range(len(self.D0))]   # unit: kg/m3
    
class node:
    def all_equal(self,iterator):
        iterator = iter(iterator)
        try:
            first = next(iterator)
        except StopIteration:
            return True
        return all(first == x for x in iterator)
    
    def read_general_info(self,xls):
        warnings.simplefilter('ignore')
        self.info = pd.read_excel(xls,'Info')
        warnings.resetwarnings() 

    def SI_conversion(self,ucalc,data,nominals,unit_dict,labels,update_colname=False):
        initialist = []
        for i in range(len(labels)): initialist.append([])
        
        unit_dictionary = dict(zip(labels,initialist))
        units = []
        for t in range(len(data)):
            standard = []
            for i in range(len(nominals)):
                if update_colname == True:
                    unit = nominals[i].split('_')[-1]
                else: unit = data[nominals[i]+'_unit'][t]
                units.append(unit)
                if 'Sm3' in unit or 'Sm³' in unit: standard.append(True)
                else: standard.append(False)
                if unit in unit_dict.keys():
                    unit = unit_dict[unit]
                tuple_list = list(zip(data[nominals[i]].tolist(),len(data)*[unit]))
                unit_dictionary[labels[i]] = [list(x) for x in tuple_list]
        ucalc.input_dictionary['key'] = unit_dictionary
        ucalc.SI_conversion(input_keys=['key'])
        # Put the converted units into the data frame and rename columns accordingly.
        for no in range(len(nominals)):
            data[nominals[no]] = [x[0] for x in ucalc.input_dictionary['key'][labels[no]]]
            if update_colname == True:
                data = data.rename(columns={nominals[no]:nominals[no].replace(units[no],ucalc.input_dictionary['key'][labels[no]][0][1])})
                if standard[no] == True:
                    data = data.rename(columns={nominals[no]:nominals[no].replace('m³','Sm³')})
            else:
                data[nominals[no]+'_unit'] = len(data)*[ucalc.input_dictionary['key'][labels[no]][0][1]]
                if standard[no] == True:
                    data[nominals[no]+'_unit'] = len(data)*[ucalc.input_dictionary['key'][labels[no]][0][1].replace('m³','Sm³')]
        return data
            
class metering:
    
    
    def read_configuration(self,xls):
        warnings.simplefilter('ignore')
        self.configuration = pd.read_excel(xls,'Configuration')
        warnings.resetwarnings()

    def read_meter_data(self,ucalc,stream,mpoint,xls,meter_type,stream_nominal_values):
        """This method assumes that there is always a temperature and pressure measurement connected to the meter."""
        warnings.simplefilter('ignore')
        unc_meter = pd.read_excel(xls,meter_type+'_unc',index_col='Input variable')
        unc_temperature = pd.read_excel(xls,'T_unc',index_col='Input variable')
        unc_pressure = pd.read_excel(xls,'P_unc',index_col='Input variable')
        warnings.resetwarnings()

        # Add meter, P and T uncertainties. 
        printer = ['units','coverage factors','distributions']
        col = ['Unit','k','Distribution']
        unc_info_meter = [sum(np.square(unc_meter['Uncertainty'][2:].to_numpy()))]
        for no in range(len(col)):
            if mpoint.all_equal(unc_meter[col[no]][2:].tolist()) == True:
                unc_info_meter.append(unc_meter[col[no]][3])
            else: raise Exception('Please make sure that all'+meter_type+'input '+printer[no]+' have identical units.')
        # Note that pressure and temperature uncertainty is read from the Misc-cell, if this is zero, then that explains the
        # the zero std distribution after monte carlo draws.
        unc_info_P = [
            unc_pressure.loc['Misc','Uncertainty'],unc_pressure.loc['Misc','Unit'],
            unc_pressure.loc['Misc','k'],unc_pressure.loc['Misc','Distribution']]
        unc_info_T = [
            unc_temperature.loc['Misc','Uncertainty'],unc_temperature.loc['Misc','Unit'],
            unc_temperature.loc['Misc','k'],unc_temperature.loc['Misc','Distribution']]
        unc_info = unc_info_meter + unc_info_P + unc_info_T
        colnames = []
        for m in ['quantity_','pressure_','temperature_']:
            for n in ['','_unit','_coverage','_distribution']:
                colnames.append(m+'unc'+n)
        data = pd.DataFrame([unc_info],columns=colnames)
           
        # Now do SI conversion of uncertainties.
        uncertainties = ['quantity_unc','pressure_unc','temperature_unc']
        unit_dict = {
            'C':'°C','F':'°F','Mpa':'MPa','barg':'bar',
            'm3h-1':'(m^3)/(h)','Sm3h-1':'(m^3)/(h)','m3h-1':'(m^3)/(h)','m3':'m^3','Sm3':'m^3',
            'm³/s':'(m^3)/(s)','Sm³/s':'(m^3)/(s)','kg/s':'(kg)/(s)','Sm^3/s':'(m^3)/(s)','Sm³':'m^3',
            }
        labels = ['quantity','pressure','temperature']
        unc_units = [x+'_unit' for x in uncertainties]

        # Convert any relative uncertainties to absolute uncertainties, then do SI conversion.
        stream_nominal_values = stream.process_data.drop(columns = ['DateTime','Uptime_h'])
        data = pd.concat([stream.process_data['DateTime'],data],axis=1)
        for i in range(len(uncertainties)):
            nominal_uncertainty = stream_nominal_values.filter(regex='^{}'.format(uncertainties[i].split('_')[0]))*data[uncertainties[i]][0]
            data[labels[i]+'_unc'] = nominal_uncertainty
            data[unc_units[i]] = [nominal_uncertainty.columns[0].split('_')[-1]]*len(data)
              
        self.transmitter_uncertainties = mpoint.SI_conversion(ucalc,data,uncertainties,unit_dict,labels)
        
        
class streamProcess:
    
    
    def read_composition_and_composition_uncertainties(self,mpoint,xls):
        warnings.simplefilter('ignore')
        composition = pd.read_excel(xls,'Composition')
        unc_composition = pd.read_excel(xls,'Composition_unc',index_col='Comp')
        warnings.resetwarnings()

        # Add the composition uncertainties to the composition frame.
        unc_composition = unc_composition.add_suffix('_unc')
        unc_Composition_info = [
            unc_composition.loc['Unit','H2_unc'], unc_composition.loc['k','H2_unc'],
            unc_composition.loc['Distribution','H2_unc']]
        unc_composition = unc_composition.loc[['Uncertainty']]
        unc_composition = unc_composition.reset_index().drop(columns='Comp')
        unc_composition = pd.concat([unc_composition,pd.DataFrame(
            [unc_Composition_info],columns=['composition_unc_unit','composition_unc_coverage','composition_unc_distribution'])],axis=1)
        composition = composition.rename(columns={"Unit": "composition_unit"})
        components = composition.columns.to_list()[2:]
        unc_components = [x+'_unc' for x in components]
        composition = pd.concat([composition,unc_composition],axis=1)

        # Convert the columns in composition frame to cmol/mol and all_measurement_input to SI.
        # Composition:
        units = ['composition_unit','composition_unc_unit']
        print_units = ['composition', 'composition uncertainty']
        wanted_units = ['\'cmol_mol-1','cmol_mol-1']
        allowed_units = ['\'mol_mol-1','mol_mol-1']
        comps = [components, unc_components]
        if mpoint.all_equal(composition['composition_unit']) == False: 
            raise Exception('Please make sure that all input composition units are the same.')
        if mpoint.all_equal(unc_composition['composition_unc_unit']) == False: 
            raise Exception('Please make sure that all input composition uncertainty units are the same.')
        for i in range(len(units)):
            if composition[units[i]][0] not in wanted_units:
                if composition[units[i]][0] in allowed_units:
                    composition[comps[i]] = composition[comps[i]]*100
                    composition[units[i]] = wanted_units[0]
                else: raise Exception('Please change the '+print_units[i]+' units to either cmol_mol-1 or mol_mol-1')
        # Also make sure that the composition column names are nice.
        old_column_names = composition.drop(columns=['DateTime']).filter(regex='^((?!composition).)*$').columns.to_list()
        new_column_names = []
        for x in old_column_names:
            new_column_names.append(x.strip().replace('nC','n-C').replace('iC','i-C').replace('neoC5','neo-C5').replace('C6','n-C6'))
        new_column_names = ["".join(x.split()) for x in new_column_names]
        self.composition = composition.rename(columns =dict(zip(old_column_names,new_column_names)))

    def set_composition_and_composition_uncertainties(self,xls=None,read_from_file=True):
        if read_from_file == True:
            if isinstance(xls,pd.io.excel._base.ExcelFile): self.read_composition(xls)
            else: raise Exception('Provide a pandas.io.excel._base.ExcelFile object if you want to read compostions from file.')
    
    def read_process_data(self,ucalc,mpoint,xls):
        warnings.simplefilter('ignore')
        process_data = pd.read_excel(xls,'ProcessData')
        warnings.resetwarnings()

        ### Do some cleaning.
        process_data = process_data.drop(columns=['Case'])
        # Do uptime h multiplication.
        uptime_multiply = ['Volume_m3','StandardVolume_Sm3','Mass_kg']
        for col_name in process_data.columns:
            if col_name in uptime_multiply:
                process_data[col_name] = process_data[col_name]*process_data['Uptime_h']
            

        # convert to SI
        # Create a dictionary that can be passed to ucalc and call ucalc.
        nominal_df = process_data.drop(columns=['DateTime']).filter(regex ='^((?!ptime).)*$' )
        nominals = nominal_df.columns.tolist()
        unit_dict = {
            'C':'°C','F':'°F','Mpa':'MPa',
            'm3h-1':'(m^3)/(h)','Sm3h-1':'(m^3)/(h)','m3h-1':'(m^3)/(h)','m3':'m^3','Sm3':'m^3'}
        labels = ['quantity','temperature','pressure']
        process_data = mpoint.SI_conversion(ucalc,process_data,nominals,unit_dict,labels,update_colname=True)
        
        # Make the nominal column values more generic for easier indexing of data frame later on.
        quantities = ['Volume_m³','StandardVolume_m³','Mass_kg','VolumeFlowRate_m³/s','StandardVolumeFlowRate_m³/s','MassFlowRate_kg/s']
        SI_converter_compatible = ['quantity_m³','quantity_Sm³','quantity_kg','quantity_(m³)/(s)','quantity_(Sm³)/(s)','quantity_(kg)/(s)']       
   
        
        for i in range(len(process_data.columns)):
            col_name = process_data.columns[i]
            if col_name in quantities:   
                process_data.rename(columns={col_name:SI_converter_compatible[i]},inplace=True)
               
        process_data = process_data.rename(columns={'P_Pa':'pressure_Pa','T_K':'temperature_K'})
        #if 'StandardVolume_Sm3' in nominals: process_data = process_data.rename(columns={'quantity_m³':'quantity_Sm^3'})
        #if 'Volume_m3' in nominals: process_data = process_data.rename(columns={'Volume_m3':'quantity_m^3'}) # Need to check that it works for rates as well.
        self.process_data = process_data


class met4H2:
    def __init__(self):
        self.all_measurement_input = pd.DataFrame()
        self.monte_carlo_draws = {}
        self.gcv_temperature = 25
        self.T0 = np.array(15+273.15) 
        self.P0 = np.array(101325)   # AGA8Detail spiser bare numpy.
        self.model_Z_noise = 0.001
        self.model_Hmol_noise = 0
        self.model_Z0_noise = 0
        self.model_MM_noise = 0
        self.model_DM_noise = 0
        self.model_DM0_noise = 0
        self.component_dictionary = dict(zip(['N2','CO2','C1','C2','C3','n-C4','i-C4','n-C5','i-C5','neo-C5','n-C6','2-Methylpentane','3-Methylpentane','2,2-Dimethylbutane','2,3-Dimethylbutane','n-C7','n-C8','n-C9','n-C10','n-C11','n-C12','n-C13','n-C14','n-C15','Ethylene','Propene','1-Butene','cis-2-Butene','trans-2-Butene','2-Methylpropene','1-Pentene','Propadiene','1,2-Butadiene','1,3-Butadiene','Acetylene','Cyclopentane','Methylcyclopentane','Ethylcyclopentane','Cyclohexane','Methylcyclohexane','Ethylcyclohexane','Benzene','Toluene','Ethylbenzene','o-Xylene','Methanol','Methanetiol','H2','H2O','H2S','NH3','HCN','CO','COS','CS2'],
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
    
    def read_settings_file(self,ucalc, stream, meter, mpoint, filename, filepath = None, meter_type = 'USM'):
        """Reads settings file and stores the information to data frames stored in the m4h2 object.
        All units are converted to SI units."""

        if filepath == None:
            xls = pd.ExcelFile(filename)
        else:
            xls = pd.ExcelFile(filepath  + "/"+ filename)

        mpoint.read_general_info(xls)
        stream.read_composition_and_composition_uncertainties(mpoint,xls)
        stream.read_process_data(ucalc,mpoint,xls)
        
        # meter.dataframes must be set after stream. Consider whether it is tidier to 
        # pass the input as argument or if one should access stream in meter. 
        meter.read_meter_data(ucalc,stream,mpoint,xls,meter_type,stream.process_data)
        meter.read_configuration(xls)
        
    def collect_measurement_input(self,stream,meter):
        # Merge the all_measurement_input and composition frame on DateTime. 
        # Fill na with previous entry where needed.
        self.all_measurement_input = stream.process_data.merge(meter.transmitter_uncertainties,left_on='DateTime',right_on='DateTime',how='outer')
        self.all_measurement_input = self.all_measurement_input.merge(stream.composition,left_on='DateTime',right_on='DateTime',how='outer')
        self.all_measurement_input = self.all_measurement_input.sort_values(by='DateTime',ignore_index=True).fillna(method="ffill").fillna(method="bfill")

       
        # Check that all coverage factors are 2 and distributions are normal. If not raise exception for now, but could add a bit in the future.
        for no in range(len(self.all_measurement_input.filter(regex='distribution').columns)):
            a = self.all_measurement_input[self.all_measurement_input.filter(regex='distribution').columns[no]].to_numpy()
            a = np.array([x.lower() for x in a])

            b = self.all_measurement_input[self.all_measurement_input.filter(regex='coverage').columns[no]].to_numpy()
            if not (a[0] == a).all() and a[0].lower() == 'normal':
                raise Exception('Please make sure that all distributions are normal.')
            if not(b[0] == b).all() and (b[0].lower() == 2):
                raise Exception('Please make sure that all coverage factors are k=2.')  
           
    def make_monte_carlo_draws(self, monte_carlo_draws=10000,frame=None, temp_lim=None, pres_lim=None, quant_lim=None): # rename this function.
        """Subsets, then operates on m4h2.all_measurement_input to create monte carlo draws
        for pressure, temperature, quantity and composition. 
        The function outputs stores the results in the object in data frames contained in a dictionary, 
        where columns are MC draws and time points are in rows. Compostion is normalized to 100.
        temp_lim, pres_lim and quant_lim are optional parameters, that if set should be lists where
        the first entry is the lower acceptable limit of the corresponding parameter, and the second 
        entry is the upper acceptable limit."""

        # Decided to not use the code in uncertaintyCalculator, because a) the code in uncertaintyCalculator
        # depends on a lot of things that are stored in the uncertaintyCalculator object, b) P,T,quantity and 
        # composition are independent measurements, and we don't need all the fancy correlation matrix stuff here.

        # Pressure, temperature and quantity do not need normalizing, that is why the draws are done in a loop,
        # on the first iteration, there is no normalizing, on the second there is.

        
        if frame is None:
            return_out = False
            # Subsetting the big input frame based on column names.
            frame1 = pd.concat([
                self.all_measurement_input.filter(regex='quantity').drop(columns=['quantity_unc_unit','quantity_unc_coverage','quantity_unc_distribution']),
                self.all_measurement_input[['pressure_Pa','pressure_unc','temperature_K','temperature_unc']]],axis=1)
            frame1.columns = ['quantity','quantity_unc','pressure','pressure_unc','temperature','temperature_unc']
            frame2 = self.all_measurement_input.filter(
                regex='^((?!quantity).)*$').filter(
                regex='^((?!pressure).)*$').filter(
                regex='^((?!temperature).)*$').filter(
                regex='^((?!ptime).)*$').filter(
                regex='^((?!composition).)*$').filter(
                regex='^((?!DateTime).)*$')
            frame_list = [frame1,frame2]
        else:
            return_out = True
            frame_list = [frame]
            
        for i in range(len(frame_list)):
            frame = frame_list[i]
            out_frames = list()
            mc.npts = monte_carlo_draws
            col_names = list(range(1, monte_carlo_draws + 1))
            col_names = [str(i) for i in col_names]
            parameter_list = [x for x in frame.columns if not 'unc' in x]
            for parameter in parameter_list:
                val = pd.concat([frame[parameter],frame[parameter+'_unc']],axis=1)
                mc_df = pd.DataFrame(columns=col_names)
                for _, row in val.iterrows():
                    generic_draws = mc.N(0,1) # Gjøre uavhengig av senere beregning lagre i objektet. 
                    dist = generic_draws._mcpts
                    new_dist = pd.DataFrame(row[parameter] + (dist - dist.mean())*(row[parameter+'_unc']/dist.std())).T
                    new_dist.columns = col_names
                    mc_df = pd.concat([mc_df,new_dist],axis=0)
                mc_df.index=list(range(len(val)))
                out_frames.append(mc_df)
            # Then add the contents in the out_frames list to a dictionary that belongs to the object.
            # If a composition was drawn, then normalize.
            if i == 1:
                sum_across_frames = sum(out_frames)
                for no in range(len(out_frames)):
                    # Set limit between 0 and 100
                    out_frames[no][out_frames[no] < 0] = 0
                    out_frames[no][out_frames[no] > 100] = 100
                    self.monte_carlo_draws[parameter_list[no]] = (out_frames[no]/sum_across_frames)*100
                    self.monte_carlo_draws[parameter_list[no]] = out_frames[no].add_prefix(parameter_list[no]+'_') 
            else:
                limits = [quant_lim,pres_lim,temp_lim]
                for no in range(len(out_frames)):
                    # Set lim between those potentially given with temp_lim, pres_lim and quant_lim
                    if isinstance(limits[no],list):
                        out_frames[no][out_frames[no] < limits[no][0]] = limits[no][0]
                        out_frames[no][out_frames[no] > limits[no][1]] = limits[no][1]
                    self.monte_carlo_draws[parameter_list[no]] = out_frames[no].add_prefix(parameter_list[no]+'_')             
        
        if return_out == True:
            return out_frames[0]
    
    def get_pressure(self,nominal,standard,nominal_values=None):
        if standard == 'AGA8Detail': divide_by = 1000 # go to kPa
        if standard == 'TREND': divide_by = 1 # keep Pa
        if nominal == True:
            pressure = nominal_values.filter(regex='^pressure')/divide_by 
        else:
            pressure = self.monte_carlo_draws['pressure'].iloc[nominal_values.index,:]/divide_by 
        pressure = pressure.rename(columns={pressure.columns[0]:'pressure_1'})
        P0 = self.P0/divide_by
        return np.array(P0), pressure

    def get_temperature(self,nominal, nominal_values=None):
        if nominal == True: temperature = nominal_values.filter(regex='^temperature')
        else: temperature = self.monte_carlo_draws['temperature'].iloc[nominal_values.index,:]
        temperature = temperature.rename(columns={temperature.columns[0]:'temperature_1'})
        return np.array(self.T0), temperature

    def get_gross_calorific_value(self,nominal,component_names,monte_carlo_draws=10000):

        # To not have to read extra files the gross calorific value is hardcoded.
        gcv_df = pd.DataFrame({"Component":['Methane','Ethane','Propane','n-Butane','2-Methylpropane','n-Pentane','2-Methylbutane','2,2-Dimethylpropane','n-Hexane','2-Methylpentane','3-Methylpentane','2,2-Dimethylbutane','2,3-Dimethylbutane','n-Heptane','n-Octane','n-Nonane','n-Decane','n-Undecane','n-Dodecane','n-Tridecane','n-Tetradecane','n-Pentadecane','Ethene','Propene','1-Butene','cis-2-Butene','trans-2-Butene','2-Methylpropene','1-Pentene','Propadiene','1,2-Butadiene','1,3-Butadiene','Ethyne','Cyclopentane','Methylcyclopentane','Ethylcyclopentane','Cyclohexane','Methylcyclohexane','Ethylcyclohexane','Benzene','Toluene','Ethylbenzene','o-Xylene','Methanol','Methanetiol','Hydrogen','Water','Hydrogen sulfide','Ammonia','Hydrogen cyanide','Carbon monoxide','Carbonyl sulfide','Carbon disulfide'],
                   "Name":['C1','C2','C3','n-C4','i-C4','n-C5','i-C5','neo-C5','n-C6','2-Methylpentane','3-Methylpentane','2,2-Dimethylbutane','2,3-Dimethylbutane','n-C7','n-C8','n-C9','n-C10','n-C11','n-C12','n-C13','n-C14','n-C15','Ethylene','Propene','1-Butene','cis-2-Butene','trans-2-Butene','2-Methylpropene','1-Pentene','Propadiene','1,2-Butadiene','1,3-Butadiene','Acetylene','Cyclopentane','Methylcyclopentane','Ethylcyclopentane','Cyclohexane','Methylcyclohexane','Ethylcyclohexane','Benzene','Toluene','Ethylbenzene','o-Xylene','Methanol','Methanetiol','H2','H2O','H2S','NH3','HCN','CO','COS','CS2'],
                   "molwt(g/mol)":[16.04246,30.06904,44.09562,58.12220,58.12220,72.14878,72.14878,72.14878,86.17536,86.17536,86.17536,86.17536,86.17536,100.20194,114.22852,128.25510,142.28168,156.30826,170.33484,184.36142,198.38800,212.41458,28.05316,42.07974,56.10632,56.10632,56.10632,56.10632,70.13290,40.06386,54.09044,54.09044,26.03728,70.13290,84.15948,98.18606,84.15948,98.18606,112.21264,78.11184,92.13842,106.16500,106.16500,32.04186,48.10746,2.01588,18.01528,34.08088,17.03052,27.02534,28.01010,60.07510,76.14070],
                   "Hcmol_0":[892.92,1564.35,2224.03,2883.35,2874.21,3542.91,3536.01,3521.75,4203.24,4195.64,4198.27,4185.86,4193.68,4862.88,5522.41,6182.92,6842.69,7502.22,8162.43,8821.88,9481.71,10141.65,1413.55,2061.57,2721.57,2714.88,2711.09,2704.88,3381.32,1945.26,2597.15,2544.14,1301.86,3326.14,3977.05,4637.2,3960.68,4609.33,5272.76,3305.12,3952.77,4613.16,4602.18,766.6,1241.64,286.64,45.064,562.93,384.57,671.92,282.8,548.01,1104.05],
                   "Hcmol_15":[891.51,1562.14,2221.1,2879.76,2870.58,3538.6,3531.68,3517.44,4198.24,4190.62,4193.22,4180.83,4188.61,4857.18,5516.01,6175.82,6834.9,7493.73,8153.24,8811.99,9471.12,10130.23,1412.12,2059.43,2718.71,2711.94,2708.26,2702.06,3377.76,1943.97,2595.12,2542.11,1301.37,3222.19,3972.46,4631.93,3956.02,4604.08,5266.9,3302.9,3949.83,4609.54,4598.64,765.09,1240.28,286.15,44.431,562.38,383.51,671.67,282.91,548.14,1104.32],
                   "Hcmol_15_55":[891.46,1562.06,2220.99,2879.63,2870.45,3538.45,3531.52,3517.28,4198.06,4190.44,4193.04,4180.65,4188.43,4856.98,5515.78,6175.56,6834.62,7493.42,8152.91,8811.63,9470.73,10129.82,1412.07,2059.35,2718.6,2711.83,2708.16,2701.96,3377.63,1943.92,2595.05,2542.03,1301.35,3222.05,3972.29,4631.74,3955.85,4603.89,5266.69,3302.81,3949.72,4609.4,4598.52,765.03,1240.23,286.13,44.408,562.36,383.47,671.66,282.91,548.15,1104.33],
                   "Hcmol_20":[891.05,1561.42,2220.13,2878.58,2869.39,3537.19,3530.25,3516.02,4196.6,4188.97,4191.56,4179.17,4186.94,4855.31,5513.9,6173.48,6832.33,7490.93,8150.21,8808.73,9467.63,10126.52,1411.65,2058.73,2717.76,2710.97,2707.33,2701.13,3376.59,1943.54,2594.46,2541.44,1301.21,3320.89,3970.95,4630.2,3954.49,4602.36,5264.97,3302.16,3948.86,4608.34,4597.48,764.59,1239.84,285.99,44.222,562.19,383.16,671.58,282.95,548.19,1104.4],
                   "Hcmol_25":[890.58,1560.69,2219.17,2877.4,2868.2,3535.77,3528.83,3514.61,4194.95,4187.32,4189.9,4177.52,4185.28,4853.43,5511.8,6171.15,6829.77,7488.14,8147.19,8805.48,9464.15,10122.82,1411.18,2058.02,2716.82,2710,2706.4,2700.2,3375.42,1943.11,2593.79,2540.77,1301.05,3319.59,3969.44,4628.47,3952.96,4600.64,5263.05,3301.43,3947.89,4607.15,4596.31,764.09,1239.39,285.83,44.013,562.01,382.81,671.5,282.98,548.23,1104.49],
                   "uHcmol":[0.19,0.51,0.51,0.72,0.72,0.23,0.23,0.25,0.32,0.53,0.53,0.48,0.46,0.67,0.76,0.81,0.87,1.54,1.13,1.21,1.32,1.44,0.21,0.34,0.39,0.5,0.47,0.42,0.73,0.6,0.4,0.41,0.32,0.36,0.56,0.71,0.32,0.71,0.95,0.27,0.51,0.66,0.76,0.13,0.32,0.02,0.004,0.23,0.18,1.26,0.06,0.24,0.43]})

        df = gcv_df.set_index("Name").reindex(component_names).fillna(0) 
    
        if nominal == False:
            H0_df = df[["Hcmol_"+str(self.gcv_temperature),"uHcmol"]].rename(columns={"uHcmol":"Hcmol_"+str(self.gcv_temperature)+"_unc"}) 
            # This is a frame that's a bit different from the others frames that can be supplied by the same method. It has components, not time along the index.
            H0_frame = self.make_monte_carlo_draws(monte_carlo_draws=monte_carlo_draws,frame=H0_df) 
            H0_frame = H0_frame.rename(index=dict(zip(list(range(len(H0_frame))),component_names))).T 
        else: 
            H0_frame = df[["Hcmol_"+str(self.gcv_temperature)]].T.set_index(pd.Index([1]))
            H0_frame.columns.name = None
        return H0_frame
    
    def call_AGA8Detail(self, T0, P0, temperature, pressure, component_names, component_values):


        AGA8_detail_order = ['NA','C1','N2','CO2','C2','C3','i-C4','n-C4','i-C5','n-C5','n-C6','n-C7','n-C8','n-C9','n-C10','H2','O2','CO','H2O','H2S','He','Ar'] 

        draws = component_values.shape[1]
        if draws > 1: nominal = False
        else: nominal = True

        
        # Initialize data frames prior to calculation. Could extend this if we want to get further output from TREND.
        if nominal == False: frame_size = draws+1
        else: frame_size = 2
        frame_list = []
        out_pars = ['Z_','Z0_','MM_','DM_','DM0_','Hmol_']
        for a in range(len(out_pars)):
            frame_list.append(pd.DataFrame(index=np.arange(component_values.shape[0]), columns=[out_pars[a]+str(i) for i in list(range(1, frame_size))]))
        # Unpack empty frames.
        Z_frame, Z0_frame, MM_frame, DM_frame, DM0_frame, Hmol_frame = frame_list

        # Get gross calorific values.
        H0_frame = self.get_gross_calorific_value(nominal,component_names,monte_carlo_draws=draws)

        
        for t in range(component_values.shape[0]):
            for draw in range(draws):
                # Parallel process here?
        
                # create a data frame (x_array) with the gas composition.
                stacked_list = [[i] for i in list(component_values[t,draw,:])]
                x_array = pd.DataFrame(dict(zip(component_names,stacked_list)))
                    
                # Drop columns that are not contained in AGA_8.  QA: What about neo-pentane?
                x_array = x_array[x_array.columns.intersection(AGA8_detail_order)]
                # Create columns that contain zeros for components in AGA which are not in H0_frame. 
                for component in AGA8_detail_order:
                    if not component in x_array.columns:
                        x_array[component] = 0
                # Re-order columns in H0_frame so they match the order of components in AGA8Detail.
                x_array = x_array[AGA8_detail_order] 
                # Then do AGA8Detail calculation
                Z_frame.iloc[t,draw] = AGA8Detail(P=np.array(pressure.iloc[t,draw]), T = np.array(temperature.iloc[t,draw]), x=x_array.iloc[0,:].to_numpy()).run().Z
                Z0_frame.iloc[t,draw] = AGA8Detail(P=P0, T = T0, x=x_array.iloc[0,:].to_numpy()).run().Z 
                MM_frame.iloc[t,draw] = AGA8Detail(P=P0, T = T0, x=x_array.iloc[0,:].to_numpy()).run().MM
                DM_frame.iloc[t,draw] = AGA8Detail(P=np.array(pressure.iloc[t,draw]), T = np.array(temperature.iloc[t,draw]), x=x_array.iloc[0,:].to_numpy()).run().D
                DM0_frame.iloc[t,draw] = AGA8Detail(P=P0, T = T0, x=x_array.iloc[0,:].to_numpy()).run().D 

                x_array = x_array[x_array.columns.intersection(H0_frame.columns)] # Get the columns that are in both x_array and H0_frame.
                for component in H0_frame.columns.to_list():
                    if not component in x_array.columns:
                        x_array[component] = 0
                x_array = x_array[H0_frame.columns.to_list()] 

                # Get Hmol by multiplying the concentration of each component (x_array) by gross calorific value (H0) and sum across components.
                Hmol_frame.iloc[t,draw] = (H0_frame.iloc[draw,:]*x_array).sum(axis = 1)[0]
            
        return Z_frame, Z0_frame, MM_frame/1000, DM_frame*1000, DM0_frame*1000, Hmol_frame
    
    def call_TREND(self,T0, P0, temperature, pressure, component_names, component_values,eqn_of_state=None):
        if eqn_of_state == None or eqn_of_state == 'GERG2008': eqn_of_state = 1
        elif eqn_of_state == 'AGA8': eqn_of_state = 7
        else: raise Exception('Please make sure to use a recognized equation of state.')
        
        # TREND doesn't always work with neopentane in the mix. So for now remove neopentane from the mix, while figuring out something smarter.
        if 'neo-C5' in component_names:
            component_values = np.delete(component_values, component_names.index('neo-C5'),axis=2)
            del component_names[component_names.index('neo-C5')]

        # Normalize
        for t in range(component_values.shape[0]):
            for draw in range(component_values.shape[1]):
                component_values[t,draw,:] = component_values[t,draw,:]/component_values[t,draw,:].sum()
             

        draws = component_values.shape[1]
        if draws > 1: nominal = False
        else: nominal = True

        # Initialize data frames prior to calculation.
        if nominal == False: frame_size = draws+1
        else: frame_size = 2
        frame_list = []
        out_pars = ['Z_','Z0_','MM_','DM_','DM0_','Hmol_']
        for a in range(len(out_pars)):
            frame_list.append(pd.DataFrame(index=np.arange(component_values.shape[0]), columns=[out_pars[a]+str(i) for i in list(range(1, frame_size))]))
        

        # Get gross calorific values.
        H0_frame = self.get_gross_calorific_value(nominal,component_names,monte_carlo_draws=draws)
        
        # REMOVE LATER - FOR TESTING PURPOSES
        #draws = 1

        coolProp_component_names = [self.component_dictionary[i][0] for i in component_names]
        TREND_component_names = [self.component_dictionary[i][1] for i in component_names]
        for draw in range(draws): 
            component_values_draw = component_values[:,draw,:]

           # REMOVE LATER - FOR TESTING PURPOSES.
            #component_values_draw = np.vstack(([0.6,0.4], [0.6,0.4], [0.6,0.4]))
            #TREND_component_names = ['Methane','Ethane']

            # Create the fluid input dictionary.
            fluid_input={
                'input' : 'TP',   # temperature and pressure
                'calctype' : 'D', # default density molar
                'fluids' : TREND_component_names, 
                'moles' : component_values_draw, # molfraksjon sum=1 
                'eos_ind' : len(TREND_component_names)*[1], # eqn state. 1=GERG2008/Helmholtz. 
                'mix_ind' : 1,
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

            # Get Hmol by multiplying the concentration of each component (component_values) by gross calorific value (H0) and sum across components.
            for t in range(len(component_values_draw)):
                frame_list[-1].iloc[t,draw] = (component_values_draw[t,:]*H0_frame.values).sum(axis=1)[0]

            # Put non-TREND-specific results into return frames.
            trend_results = [tr.Z,tr.Z0,tr.MW,tr.DM,tr.DM0] 
            for n in range(len(frame_list[:-1])):
                frame = frame_list[:-1][n]
                frame[frame.columns[draw]] = trend_results[n]
        
        
        # Save TREND specific output directly to object (check with Kjetil which ones). Return the rest to calculate some more stuff.
        # dielktrisk konst, osv. 

        return frame_list
    
    def compute_gas_properties(self,method='TREND',nominal=False,eqn_of_state=None):
        """If nominal=False, Z-s and H-es are calculated using monte carlo draws and output is saved to self.main_m4h2_output. 
        If nominal=True, monte carlo is not used and the output is saved to self.main_nominal_m4h2_output. 
        Note that uncertainties in Z and H are NOT computed in the latter case, not even model uncertainties are added to the nominal values."""

        if not method in ['AGA8Detail','TREND']: raise Exception('Please make sure to apply a supported method for computing gass properties.')
        if nominal == False: monte_carlo_draws = self.monte_carlo_draws['temperature'].shape[1]
        else: monte_carlo_draws = 1
        

        # get a data frame with nominal values.
        nominal_values = self.all_measurement_input.filter(
                regex='^((?!ptime).)*$').filter(
                regex='^((?!composition).)*$').filter(
                regex='^((?!DateTime).)*$').filter(
                regex='^((?!unc).)*$')
        
        
        # get gas component names.
        component_names = [comp for comp in nominal_values.columns.to_list() if not any(non_comp in comp for non_comp in ['quantity','pressure','temperature'])]
            
        

        # Find all duplicate rows in nominal values and create a dictionary with keys that keep the first index, and values that are lists of duplicate rows.
        new_nominal_values = nominal_values.drop_duplicates(keep='first')
        keys = nominal_values.drop_duplicates(keep='first').index.to_list()
        values = nominal_values[nominal_values.duplicated()].index.to_list()
        duplicates = {}
        for keyno in range(len(keys)):
            duplicates[keys[keyno]] = []
            for idx in values: 
                if nominal_values.iloc[idx,:].equals(nominal_values.iloc[keys[keyno],:]):
                    duplicates[keys[keyno]].append(idx)

        
        # Concatenate data frames to 3D numpy arrays for easy collection of desired arrays.
        composition_frames = []
        for comp in component_names:
            if nominal == False: 
                composition_frames.append(self.monte_carlo_draws[comp])
                stacked = np.zeros([len(nominal_values),monte_carlo_draws])
            else: 
                composition_frames.append(pd.DataFrame(nominal_values[comp]))
                stacked = np.zeros([len(nominal_values),1])
        for comp_no in range(len(component_names)):
            stacked = np.dstack((stacked,composition_frames[comp_no].to_numpy()))
        stacked = stacked[:,:,1:]/100 # because composition is cmol/mol in the all_measurement_input
        stacked = stacked[list(duplicates.keys()),:,:]

        # To check time efficiency of the computation below.
        start = time.time()

        # temperature and pressure.
        T0, temperature = self.get_temperature(nominal,nominal_values=new_nominal_values)
        P0, pressure = self.get_pressure(nominal,method,nominal_values=new_nominal_values)
        if nominal == True: 
            reshaped = new_nominal_values[component_names].to_numpy().reshape((new_nominal_values[component_names].to_numpy().shape[0], new_nominal_values[component_names].to_numpy().shape[1], 1))
            component_values = np.moveaxis(reshaped, [2,1], [1,2])
        else: component_values = stacked

        if method == 'TREND': 
            Z_frame, Z0_frame, MM_frame, DM_frame, DM0_frame, Hmol_frame = self.call_TREND(T0, P0, temperature, pressure, component_names, component_values,eqn_of_state=eqn_of_state)
        if method == 'AGA8Detail':
            Z_frame, Z0_frame, MM_frame, DM_frame, DM0_frame, Hmol_frame = self.call_AGA8Detail(T0, P0, temperature, pressure, component_names, component_values)

    
        # outputs common to AGA8Detail and TREND are dealt with here.
        # Add model noise to parameters (noise level determined in __init__(self))
        noisy_list = [self.model_Z_noise, self.model_Z0_noise, self.model_MM_noise, self.model_DM_noise, self.model_DM0_noise,self.model_Hmol_noise]
        frame_list = [Z_frame, Z0_frame, MM_frame, DM_frame, DM0_frame, Hmol_frame]
        for i in range(len(noisy_list)):
            if nominal == False: 
                if noisy_list[i] > 0:
                    noise_unc = []
                    for n in frame_list[i].mean(axis=1).to_list():
                        noise_unc.append(n*noisy_list[i])
                    noise_frame = pd.DataFrame({'noise':len(frame_list[i])*[0],'noise_unc':noise_unc})
                    noise = self.make_monte_carlo_draws(monte_carlo_draws=monte_carlo_draws,frame=noise_frame)
                else: 
                    noise = pd.DataFrame(0, index=range(len(frame_list[i])), columns=list(range(1,monte_carlo_draws+1)))
                    noise = noise.astype('object')
                # Keep only numbers in column number to be able to use the dataframes like we want to below.
                frame_list[i].columns = [x.split('_')[-1] for x in frame_list[i].columns.to_list()]
                frame_list[i] = pd.DataFrame(frame_list[i].values+noise.values) # Adding noise (it is either made in the above if, or else zero.)
            else: frame_list[i].columns = ['0']
        # Unpacking frame list.
        _, _, MM_df, DM_df, DM0_df, Hmol_df = frame_list

        # Calculate densities and calorific values.
        D_df=pd.DataFrame(DM_df.values*MM_df.values)   # D = DM*MM      #density kg/m3 
        D0_df=pd.DataFrame(DM0_df.values*MM_df.values) # D0 = DM0*MM    #density kg/m3
        Hv_df=pd.DataFrame(Hmol_df.values/(MM_df.values/D0_df.values))  # Hv = (Hmol/(MM/D0)) 
        Hm_df=pd.DataFrame(Hmol_df.values/MM_df.values)# Hm = Hmol/MM

        # Now fillna according to duplicated input and add prefixes to be able to understand what's in the data frames later. 
        calculated = frame_list+[D_df,D0_df,Hv_df,Hm_df]
        frame_names = ['Z_','Z0_','MM_','DM_','DM0_','Hg_','D_','D0_','Hv_','Hm_']
        for n in range(len(frame_names)):
            df = pd.DataFrame(columns=calculated[n].columns, index=list(range(max([max(list(duplicates.keys())),max(max(duplicates.values()))])+1))).add_prefix(frame_names[n])
            for keyno in range(len(duplicates)):
                df.iloc[list(duplicates.keys())[keyno],:] = calculated[n].iloc[keyno,:]
                for idx in list(duplicates.values())[keyno]:
                    df.iloc[idx,:] = calculated[n].iloc[keyno,:]
            calculated[n] = df

        elapsed_time = (time.time() - start) 
        # Put output in dictionary and save to object.
        if nominal == False: 
            self.main_m4h2_output = dict(zip(['Z','Z0','MM','DM','DM0','Hg','D','D0','Hv','Hm'],calculated))
            return self.main_m4h2_output
        else: 
            self.main_nominal_m4h2_output = dict(zip(['nominal_Z','nominal_Z0','nominal_MM','nominal_DM','nominal_DM0','nominal_Hg','nominal_D','nominal_D0','nominal_Hv','nominal_Hm'],calculated))
            return self.main_nominal_m4h2_output
           
    def initialize_uncertainty_calculator(self,ucalc,gas_properties,computation_method='lhs_montecarlo',monte_carlo_draws=1000):

        # check what sort of flow 
        unitstring = self.all_measurement_input.filter(regex='^((?!unc).)*$').filter(regex='^quantity').columns.str.split('_')[0][-1]
        possible_units = {'m³/s':'volume_flow','Sm³/s':'standard_vol_flow','kg/s':'mass_flow','m³':'volume','Sm³':'standard_vol','kg':'mass'}
        flowtype = possible_units[unitstring]

        # Define general input for a sample functional relationship with 3 input variables and two time points
        labels_and_units = ['P0|Pa','T0|K','P|Pa','T|K','quantity|(m^3)/(s)','Z|-','Z0|-','Hv|(J)/(m^3)']     # Define input variables and units for each input variable
        functional_relationship = '(P*T0*Z0*Hv*quantity)/(P0*T*Z)'         # Define functional relationship                                         
        if flowtype == 'volume_flow':                                                                                            
            output_label_and_unit = ['qE|(J)/(s)']                                                       # Define output variable and corresponding unit
        if flowtype == 'standard_vol':
            output_label_and_unit = ['E|J'] 
        if flowtype != 'volume_flow' and flowtype != 'standard_vol':
            raise Exception('Code isn\'t finished, please use volume flow as quantity input.')                                                    
        
        # Define values for each input variable and time point
        P0 = [float(self.P0)]*len(self.all_measurement_input)
        T0 = [float(self.T0)]*len(self.all_measurement_input)
        P = self.all_measurement_input.filter(regex='^((?!unc).)*$').filter(regex='^pressure').iloc[:,0].to_list()
        T = self.all_measurement_input.filter(regex='^((?!unc).)*$').filter(regex='^temperature').iloc[:,0].to_list()
        quantity = self.all_measurement_input.filter(regex='^((?!unc).)*$').filter(regex='^quantity').iloc[:,0].to_list()
        Z = gas_properties['Z'].mean(axis=1).to_list()
        Z0 = gas_properties['Z0'].mean(axis=1).to_list()
        Hv = gas_properties['Hv'].mean(axis=1).to_list()
        values = [P0,T0,P,T,quantity,Z,Z0,Hv]      

        # Define timestamps for each time point
        time_stamps = self.all_measurement_input['DateTime']

        # Set input and output quantities and functional relationship in uncertainty calculator.
        # It is important to define the input quantities first.
        ucalc.set_input_quantities(labels_and_units,values,time_stamps=time_stamps)
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
        input_expanded_uncertainties = [P0_u,T0_u,P_u,T_u,quantity_u,Z_u,Z0_u,Hv_u]  

        distributions =  ['normal']*len(labels_and_units)

        # Set uncertainty distributions in uncertainty calculator
        ucalc.set_input_expanded_uncertainties(input_expanded_uncertainties,absolute=True,distributions=distributions)

        # Define correlation matrix in uncertainty calculator
        forbidden_correlations=[
            ['P','T'],['P0','T0'],['P0','T'],['P0','P'],['P0','Hv'],['P0','Z0'],
            ['P0','Z'],['T0','T'],['T0','P'],['T0','Hv'],['T0','Z0'],['T0','Z'],
            ['quantity','T0'],['quantity','T'],['quantity','P0'],['quantity','Hv'],
            ['quantity','Z0'],['quantity','Z']]
        P0_frame = pd.DataFrame(float(self.P0),index = np.arange(len(self.all_measurement_input)), columns=['P0_'+str(i) for i in range(len(gas_properties['Z'].columns))])
        T0_frame = pd.DataFrame(float(self.T0),index = np.arange(len(self.all_measurement_input)), columns=['T0_'+str(i) for i in range(len(gas_properties['Z'].columns))])
        frame_list = [
            P0_frame,T0_frame,self.monte_carlo_draws['pressure'],
            self.monte_carlo_draws['temperature'],self.monte_carlo_draws['quantity'],
            gas_properties['Z'],gas_properties['Z0'],gas_properties['Hv']] 
        ucalc.input_dictionary['list_correlation_matrix'], _ = self.make_input_correlation_matrices(frame_list ,forbidden_correlations)

        # Set computational method in uncertainty calculator
        ucalc.set_computational_method(computation_method=computation_method,monte_carlo_draws=monte_carlo_draws)

        # Perform SI conversion and a few other partial calculations and add these to the object.
        # This function can also be used to add an input dictionary directly to the object, but we already did that in the previous method calls.
        # It is not necessary to do SI_conversion, as we have already done that when we set the stream and meter objects.
        ucalc.initialize_calculator(SI_conversion=True) 
        
    def make_input_correlation_matrices(self, frame_list, forbidden_correlations=None):
        """outputs upper triangular matrix for each time point contained in the data frames in frame_list"""
        if forbidden_correlations == None: forbidden_correlations = list()
        out_matrices = list()
        for t in range(len(frame_list[0])):
            out_labels = list()
            out_matrix = list()
            for elm_no in range(len(frame_list)):
                new_frame_list = frame_list[elm_no+1:]
                for new_elm_no in range(len(new_frame_list)):
                    lst = list([
                        frame_list[elm_no].columns[0][:frame_list[elm_no].columns[0].rindex('_')],
                        new_frame_list[new_elm_no].columns[0][:new_frame_list[new_elm_no].columns[0].rindex('_')]])
                    out_labels.append(lst)
                    # do correlation
                    if lst in forbidden_correlations or list(reversed(lst)) in forbidden_correlations:
                        out_matrix.append(0)
                    else:
                        frame_numpy = frame_list[elm_no].iloc[t,:].to_numpy()
                        new_frame_numpy = new_frame_list[new_elm_no].iloc[t,:].to_numpy()
                        corr = pd.Series(frame_numpy).astype('float64').corr(pd.Series(new_frame_numpy).astype('float64'))
                        if np.isnan(corr): corr = 0
                        out_matrix.append(corr)
            out_matrices.append(out_matrix)
        return out_matrices, out_labels

   
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

        

  ### Kladd  
    def parallel_proc_function_AGA8(self,time,draw,pressure,temperature,P0,T0,stacked,component_names,AGA8_detail_order,frame_list): 
        Z_frame, Z0_frame, MM_frame, DM_frame, DM0_frame = frame_list
        stacked_list = [[i] for i in list(stacked[time,draw,:])]
        x_array = pd.DataFrame(dict(zip(component_names,stacked_list)))
        x_array = x_array[x_array.columns.intersection(AGA8_detail_order)] 
        for component in AGA8_detail_order:
            if not component in x_array.columns:
                x_array[component] = 0
        x_array = x_array[AGA8_detail_order]
        Z_frame.iloc[time,draw] = AGA8Detail(P=np.array(pressure.iloc[time,draw]), T = np.array(temperature.iloc[time,draw]), x=x_array.iloc[0,:].to_numpy()).run().Z
        Z0_frame.iloc[time,draw] = AGA8Detail(P=P0, T = T0, x=x_array.iloc[0,:].to_numpy()).run().Z 
        MM_frame.iloc[time,draw] = AGA8Detail(P=P0, T = T0, x=x_array.iloc[0,:].to_numpy()).run().MM
        DM_frame.iloc[time,draw] = AGA8Detail(P=np.array(pressure.iloc[time,draw]), T = np.array(temperature.iloc[time,draw]), x=x_array.iloc[0,:].to_numpy()).run().D
        DM0_frame.iloc[time,draw] = AGA8Detail(P=P0, T = T0, x=x_array.iloc[0,:].to_numpy()).run().D 

        return Z_frame, Z0_frame, MM_frame, DM_frame, DM0_frame, x_array

 
