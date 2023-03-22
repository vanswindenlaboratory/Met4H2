import numpy as np

import CoolProp.CoolProp as cp

import os 
import importlib.util 
import sys 

for dirpath, dirnames, filenames in os.walk("."):
    for filename in [f for f in filenames if f.endswith("fluid.py")]:
        relative_path = os.path.join(dirpath, filename)

module_name = 'fluid' 
file_path = os.path.join(os.getcwd(),relative_path) 
spec = importlib.util.spec_from_file_location(module_name, file_path) 
module = importlib.util.module_from_spec(spec) 
sys.modules[module_name] = module 
spec.loader.exec_module(module)

from fluid import fluid





class trend:
    """This class contains methods to facilitate calling of the TREND_EOS 
    wrapper and calculation of a few basic gas properties, such as 
    densities, molar densities, compressibility factors and the 
    dielectric constant."""
    
    def __init__(self,fluid_input):
        self.fluid_input=fluid_input  # Dictionary containing info 
                                      # for the EoS.
        self.cp_fluids = []           # Names of relevant fluids in format 
                                      # that is recognized by CoolProp. 
                                      # This is set outside the object, 
                                      # prior to calling functions that
                                      # live in trend.
        self.D=[]                     # D, D0, Z, Z0, and dielectric_const 
                                      # are densities and compressibility
        self.D0=[]                    # factors and diel. const. computed 
        self.Z=[]                     # with the EoS.
        self.Z0=[]
        self.dielectric_constant=[]
        self.P0 = np.array(101325)    # P0 and T0 are pressure (Pa) and 
        self.T0 = np.array(273.15+15) # temperature (K) at reference 
                                      # conditions.
        
        
    def make_fluid(self,fluid_input):
        """To run calculate the equation of state the fluid object must be set first.
        This function creates and returns the fluid object from input stored in a 
        dictionary. The dictionary follows this pattern:
        fluid_input={
                'input' : 'TP',                   # temp and pressure
                'calctype' : 'D',                 # default density molar
                'fluids' : TREND_component_names, # components in mix
                'moles' : component_values_draw,  # molfraction, sum=1 
                'eos_ind' : len(TREND_component_names)*[1], 
                                                  # EoS. 1=GERG2008
                'mix_ind' : 1,                    # defining mixing rule
                'path': trend_path,               # path
                'unit': 'molar',                  # unit
                'dll_path' : trend_dll_path}      # path
        """
        fld=fluid(fluid_input['input'],
                  fluid_input['calctype'],
                  fluid_input['fluids'],
                  fluid_input['moles'],
                  fluid_input['eos_ind'],
                  fluid_input['mix_ind'],
                  fluid_input['path'],
                  fluid_input['unit'],
                  fluid_input['dll_path'])
        return fld


    def calc_EOS(self,calc_property,temperature,pressure):
        """Inputs are a gas property to be calculated (indicated by a 
        string), temperatures and pressures (lists) at line conditions. 
        Output is a list of values of the gas properties at each 
        time point."""

        gas_property=[]

        # Updating the gas property to be calculated in fluid_input,
        # and creating a copy, because the composition could change with 
        # time.
        self.fluid_input['calctype']=calc_property
        self.fluid_copy = self.fluid_input.copy()

        comp_len =len(self.composition)
        for counter in range(comp_len): 
            # Updating composition.
            self.fluid_copy['moles']=self.composition[counter]
            # Creating the fluid object.
            self.fld=self.make_fluid(self.fluid_copy)
            # The length of temperature and composition will match if we 
            # are at line conditions. But at reference conditions 
            # temperature and pressure have length of 1.
            if len(temperature)==comp_len:
                val,_ = self.fld.TREND_EOS(
                    temperature[counter],pressure[counter])
            else: 
                val,_ = self.fld.TREND_EOS(temperature[0],pressure[0])
            gas_property.append(val)
        return gas_property
    
  
    def cal_molecular_weight(self):
        """Reads composition, names of components and number of moles 
        for each component from self. Outputs total molecular weight for 
        the fluid, as well as the molecular weight of one mole of the 
        individual components."""
        
        # Getting names of fluids in format recognized by CoolProp
        fluids=self.cp_fluids
        
        mw = []
        total_mw=[]
        
        for _ in range(len(self.composition)):
            mwt=[]
            vals=[]
            for fluid in range(len(fluids)):
                # Using CoolProp to get the molecular weight.
                val=cp.PropsSI('molar_mass',fluids[fluid]) 
                mwt.append(val) 
                val=mwt[fluid]*self.fld.moles[fluid] 
                vals.append(val) 
            total_mw.append(sum(vals)) 
            mw.append(mwt)
        return total_mw, mw 


    def calc_properties(self,temperature,pressure):
        """Inputs are temperature and pressure (lists) at line conditions. 
        Composition is retrieved from self, and this information is then 
        used to calculate gas properties, which are stored in self."""

        self.pressure=pressure
        self.temperature=temperature
        self.composition=self.fluid_input['moles']
        
        # The input is always in SI units, 
        # but TREND requires pressure in MPa.
        pressure = pressure/1000000 # go to MPa
        P0=[self.P0/1000000] # go to MPa
        T0=[float(self.T0)]
        
        # Calculating gas properties.
        self.DM=self.calc_EOS('D',temperature,pressure)# unit: mol/m3
        self.DM0=self.calc_EOS('D',T0,P0)              # unit: mol/m3
        self.Z=self.calc_EOS('Z',temperature,pressure) # unit: -
        self.Z0=self.calc_EOS('Z',T0,P0)               # unit: -
        self.dielectric_constant=self.calc_EOS('DE',temperature,pressure)
                                                       # unit: -
        
        self.total_mw, self.mw = self.cal_molecular_weight() 
                            # unit total_mw: kg/mol summed over components.
                            # unit mw: kg/mol individual components.
        self.D=[self.MW[i]*self.DM[i] for i in range(len(self.D))]      
                            # unit: kg/m3
        self.D0=[self.MW[i]*self.DM0[i] for i in range(len(self.D0))]   
                            # unit: kg/m3
 