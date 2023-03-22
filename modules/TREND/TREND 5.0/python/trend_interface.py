# TREND 5.0 Python Interface 
# Example Call
# @authors: David Celný , Sven Pohl
# Bochum, 16.05.2019
from fluid import *

# For path seperation use D:\\...\\
path = 'p:\\600\\60020\\104652-Met4H2- EURAMET\\1. Arbeidsfiler\\TREND\\TREND 5.0\\'

# Path to the TREND 5.0 DLL
dll_path = 'p:\\600\\60020\\104652-Met4H2- EURAMET\\1. Arbeidsfiler\\TREND\\TREND 5.0\\TREND_x64.dll'

# Create object of fluid
fld = fluid('TP','D',['water'],[1],[1],1,path,'molar',dll_path)

# Calculate density with temperature and pressure input
dense,err = fld.TREND_EOS(300,1)
print('Density: ',dense)
print('Error: ',err)