<<<<<<< HEAD
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 17:49:23 2025

@author: Nate Whitsett
"""

import pandas as pd
from scipy.stats import vonmises
import numpy as np
import ardor.SPI_Forward_Models.SPI_Simulation as SPI
import ardor.SPI_Forward_Models.Orbit_Model_Library as OML
from matplotlib import pyplot as plt
output = pd.DataFrame()
# kappas = [0,0.25,0.5,1,2,4,8]
# for kappa in kappas:
#     x = np.linspace(-np.pi, np.pi, num=100)
#     outx = np.linspace(0, 1, num=100)
#     y = vonmises.pdf(x, kappa=kappa, loc=0)*2*np.pi
#     output['Data_'+str(kappa) +'x'] = outx
#     output['Data_'+str(kappa) +'y'] = y
#     output['Theta_x'] = x
    
    
ratios = [ 0.5, 1, 10]
es = [0, 0.05, 0.2, 0.5]
for ratio in ratios:
    for e in es:
        outx, y=  SPI.SPI_Cubic(ratio, np.pi, e, 0.05, 5, length=500)
        output['Data_e'+str(e) + '_' + str(ratio) +'x'] = outx
        output['Data_e'+str(e) + '_' + str(ratio) +'y'] = y
        output['Cubic_Theta_x'] = np.linspace(0,2*np.pi,num=len(y))
        plt.plot(outx,y )
plt.show()
output.to_csv('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/Cubic_Data.csv', index=False)
=======
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 17:49:23 2025

@author: Nate Whitsett
"""

import pandas as pd
from scipy.stats import vonmises
import numpy as np
import ardor.SPI_Forward_Models.SPI_Simulation as SPI
import ardor.SPI_Forward_Models.Orbit_Model_Library as OML
output = pd.DataFrame()
# kappas = [0,0.25,0.5,1,2,4,8]
# for kappa in kappas:
#     x = np.linspace(-np.pi, np.pi, num=100)
#     outx = np.linspace(0, 1, num=100)
#     y = vonmises.pdf(x, kappa=kappa, loc=0)*2*np.pi
#     output['Data_'+str(kappa) +'x'] = outx
#     output['Data_'+str(kappa) +'y'] = y
#     output['Theta_x'] = x
    
    
ratios = [1, 10, 100]
es = [0, 0.05, 0.2, 0.5]
star = OML.Star(1, 1, 1, 1, 1, 1, B=1, alfven=0.1)
planet = OML.Planet(1, 5, 0.05, 0.1, 1)
for ratio in ratios:
    for e in es:
        y, outx =  SPI.SPI_Cubic(ratio, 2.14281, e, star, planet, length=500)
        output['Data_e'+str(e) + '_' + str(ratio) +'x'] = outx
        output['Data_e'+str(e) + '_' + str(ratio) +'y'] = y
        output['Cubic_Theta_x'] = np.linspace(0,2*np.pi,num=len(y))
output.to_csv('C:/Users/Nate Whitsett/Desktop/Cubic_Data.csv')
>>>>>>> 843bcfd59f29e0297b65b05360ddcf3ed4f7ad77
