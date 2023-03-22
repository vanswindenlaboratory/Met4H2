from __future__ import print_function
import ctypes as ct
import sys
import platform
import os
import time

input_length = 12
fluid_length = 30
path_length = 255
unit_length = 20

class fluid:
    def __init__(self,input,calctype,fluids,moles,eos_ind,mix_ind,path,unit,dll_path):
        input = str.encode(input)
        self.input = ct.create_string_buffer(input.ljust(input_length),input_length)
        calctype = str.encode(calctype)
        self.calctype = ct.create_string_buffer(calctype.ljust(input_length),input_length)

        #fluid
        fluids_type = (ct.c_char * fluid_length) * fluid_length
        fluids_tmp = fluids_type()

        for i in range(0,fluid_length):
            fluids_tmp[i] = ct.create_string_buffer(b" ".ljust(fluid_length),fluid_length)

        for fid,fluid in enumerate(fluids):
            fluids_tmp[fid] = ct.create_string_buffer(str.encode(fluid.ljust(fluid_length)),fluid_length)
        self.fluid = fluids_tmp
        
        #moles
        moles_tmp = (fluid_length * ct.c_double)()
        for mid,mole in enumerate(moles):
            moles_tmp[mid] =  mole
        self.moles = moles_tmp

        #eos_ind
        eos_ind_tmp = (fluid_length * ct.c_int)()
        for eid,eos_id in enumerate(eos_ind):
            eos_ind_tmp[eid] =  eos_id
        self.eos_ind = eos_ind_tmp
         
        self.mix_ind = ct.c_int(mix_ind)

        path = str.encode(path)
        self.path = ct.create_string_buffer(path.ljust(path_length),path_length)
        unit = str.encode(unit)
        self.unit = ct.create_string_buffer(unit.ljust(unit_length),unit_length)

        if fid==mid and mid==eid:
            pass
        else:
            raise ValueError('NOT SAME LENGTHS OF INPUTS')
        
        self.handle_ptr = ct.POINTER(ct.c_int)()
        self.dll_path = dll_path
        assert(os.path.exists(self.dll_path))
        self.TrendDLL = ct.windll.LoadLibrary(self.dll_path)

    def TREND_EOS(self,pr1,pr2):
        errorflag = ct.c_int(0)
        self.TrendDLL.TREND_EOS_STDCALL.restype = ct.c_double # !beware required line otherwise you get unsensible output
        return self.TrendDLL.TREND_EOS_STDCALL(self.calctype, self.input, ct.byref(ct.c_double(pr1)), ct.byref(ct.c_double(pr2)), self.fluid, self.moles, self.eos_ind, ct.byref(self.mix_ind), self.path, self.unit, ct.byref(errorflag), ct.byref(self.handle_ptr), 12, 12, 30, 255, 20),errorflag
    


