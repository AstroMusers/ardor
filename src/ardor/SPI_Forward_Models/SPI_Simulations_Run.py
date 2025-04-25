<<<<<<< HEAD
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 13:45:06 2024

@author: whitsett.n
"""
import ardor.Flares.Flare as Flare
import ardor.SPI_Forward_Models.SPI_Simulation as SPI
import ardor.SPI_Forward_Models.Orbit_Model_Library as OML
import os
import numpy as np

def Run_Inverse_Cubic_Sim():
    ############# G TYPE SIMULATION #############
    base_dirG = 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/G_Type_LC'
    base_dirM = 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/M_Type_LC'

    G_star = OML.Star(0.896, 1, 1, radius = 0.92165, age = 4.6e9, B = 1, alfven=0.1)
    M_star = OML.Star(0.59, 1, 1, radius = 0.603, age = 4.6e9, B = 100, alfven=0.1)
    lc_dir = 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/G_Type_LC/219101992/tess2019253231442-s0016-0000000219101992-0152-s_lc.fits'
    lc_dir_M= 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/M_Type_LC/219822564/tess2019198215352-s0014-0000000219822564-0150-s_lc.fits'
    G_output = 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/G_Type_Inverse_Cubic/Flares'
    M_output = "C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/M_Type_Inverse_Cubic/Flares"
    lc_count = 25
    lc = Flare.tier0(lc_dir)
    # planet_G = OML.Planet(1, None, 0.055, 0, B=10, orbit_length=10, star = G_star)
    # x, model = SPI.SPI_Cubic(10, 0,0, 0.055, 5, length = 500)
    # ### e = 0.0, G_Type
    # for hosts in range(50):
    #     for iterations in range(lc_count):
    #         lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
    #         flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
    #         Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), 
    #                     output_dir = G_output, obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.csv', planet_epoch= epoch)
    #     print(hosts)
        
    planet_G = OML.Planet(1, None, 0.055, 0.05, B=10, orbit_length=10, star = G_star)
    model, x = SPI.SPI_Cubic(10, 0, 0.05, 0.055, 5, length = 500)
    ### e = 0.05
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), 
                        output_dir = G_output, obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.05_B_10.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_G = OML.Planet(1, None, 0.055, 0.2, B=10, orbit_length=1000, star = G_star)
    model, x = SPI.SPI_Cubic(10, 0, 0.2, 0.055, 5, length = 500)
    ### e = 0.2
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), 
                        output_dir = G_output, obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.2_B_10.csv', planet_epoch= epoch)
        print(hosts)
    
    
    ### e = 0.5
    planet_G = OML.Planet(1, None, 0.055, 0.5, B=10, orbit_length=1000, star = G_star)
    model, x = SPI.SPI_Cubic(10, 0, 0.5, 0.055, 5, length = 500)
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), 
                        output_dir = G_output, obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.5_B_10.csv', planet_epoch= epoch)
        print(hosts)
    
        
    ############# M TYPE SIMULATION #############
    
    planet_M = OML.Planet(1, None, 0.048, 0, B=10, orbit_length=10, star = M_star)
    model, x = SPI.SPI_Cubic(0.1, 0, 0, 0.055, 5, length = 500)
    ### e = 0.0, M_Type
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), 
                        output_dir = M_output, obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0_B_10.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_M = OML.Planet(1, None, 0.048, 0.05, B=10, orbit_length=10, star = M_star)
    model, x = SPI.SPI_Cubic(0.1, 0, 0.05, 0.055, 5, length = 500)
    ### e = 0.05
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), 
                        output_dir = M_output, obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.05_B_10.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_M = OML.Planet(1, None, 0.048, 0.2, B=10, orbit_length=10, star = M_star)
    model, x = SPI.SPI_Cubic(0.1, 0, 0.2, 0.055, 5, length = 500)
    ### e = 0.2
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), 
                        output_dir = M_output, obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.2_B_10.csv', planet_epoch= epoch)
        print(hosts)
    
    
    ### e = 0.5
    planet_M = OML.Planet(1, None, 0.048, 0.5, B=10, orbit_length=10, star = M_star)
    model, x = SPI.SPI_Cubic(0.1, 0, 0.5, 0.055, 5, length = 500)
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), 
                        output_dir = M_output, obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.5_B_10.csv', planet_epoch= epoch)
        print(hosts)
    
    
    ############# G TYPE SIMULATION SATURATED DYNAMO #############
    
        
    planet_G = OML.Planet(1, None, 0.055, 0.05, B=100, orbit_length=10, star = G_star)
    model, x = SPI.SPI_Cubic(100, 0, 0.05, 0.055, 5, length = 500)
    ### e = 0.05
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts),
                        output_dir = G_output, obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.05_B_100.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_G = OML.Planet(1, None, 0.055, 0.2, B=100, orbit_length=10, star = G_star)
    model, x = SPI.SPI_Cubic(100, 0, 0.2, 0.055, 5, length = 500)
    ### e = 0.2
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), 
                        output_dir = G_output, obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.2_B_100.csv', planet_epoch= epoch)
        print(hosts)
    
    
    ### e = 0.5
    planet_G = OML.Planet(1, None, 0.055, 0.5, B=100, orbit_length=1000, star = G_star)
    model, x = SPI.SPI_Cubic(100, 0, 0.5, 0.055, 5, length = 500)
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), 
                        output_dir = G_output, obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.5_B_100.csv', planet_epoch= epoch)
        print(hosts)
    
        
    ############# M TYPE SIMULATION SATURATED DYNAMO #############
        
    planet_M = OML.Planet(1, None, 0.048, 0.05, B=100, orbit_length=1000, star = M_star)
    model, x = SPI.SPI_Cubic(1, 0, 0.05, 0.055, 5, length = 500)
    ### e = 0.05
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), 
                        output_dir = M_output, obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.05_B_100.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_M = OML.Planet(1, None, 0.048, 0.2, B=100, orbit_length=1000, star = M_star)
    model, x = SPI.SPI_Cubic(1, 0, 0.2, 0.055, 5, length = 500)
    ### e = 0.2
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), 
                        output_dir = M_output, obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.2_B_100.csv', planet_epoch= epoch)
        print(hosts)
    
    
    ### e = 0.5
    planet_M = OML.Planet(1, None, 0.048, 0.5, B=100, orbit_length=1000, star = M_star)
    model, x = SPI.SPI_Cubic(1, 0, 0.5, 0.055, 5, length = 500)
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), 
                        output_dir = M_output, obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.5_B_100.csv', planet_epoch= epoch)
        print(hosts)
    
    
    
    
    
    ############# G TYPE SIMULATION WEAK FIELD #############
    
    
    planet_G = OML.Planet(1, None, 0.055, 0.05, B=1, orbit_length=1000, star = G_star)
    model, x = SPI.SPI_Cubic(1, 0, 0.05, 0.055, 5, length = 500)
    ### e = 0.05
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), 
                        output_dir = G_output, obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.05_B_1.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_G = OML.Planet(1, None, 0.055, 0.2, B=1, orbit_length=1000, star = G_star)
    mmodel, x = SPI.SPI_Cubic(1, 0, 0.2, 0.055, 5, length = 500)
    ### e = 0.2
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), 
                        output_dir = G_output, obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.2_B_1.csv', planet_epoch= epoch)
        print(hosts)
    
    
    ### e = 0.5
    planet_G = OML.Planet(1, None, 0.055, 0.5, B=1, orbit_length=1000, star = G_star)
    model, x = SPI.SPI_Cubic(1, 0, 0.5, 0.055, 5, length = 500)
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), 
                        output_dir = G_output, obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.5_B_1.csv', planet_epoch= epoch)
        print(hosts)
    
        
    ############# M TYPE SIMULATION WEAK FIELD #############
        
    planet_M = OML.Planet(1, None, 0.048, 0.05, B=1, orbit_length=1000, star = M_star)
    model, x = SPI.SPI_Cubic(0.01, 0, 0.05, 0.055, 5, length = 500)
    ### e = 0.05
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), 
                        output_dir = M_output, obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.05_B_1.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_M = OML.Planet(1, None, 0.048, 0.2, B=1, orbit_length=1000, star = M_star)
    model, x = SPI.SPI_Cubic(0.01, 0, 0.2, 0.055, 5, length = 500)
    ### e = 0.2
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), 
                        output_dir = M_output, obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.2_B_1.csv', planet_epoch= epoch)
        print(hosts)
    
    
    ### e = 0.5
    planet_M = OML.Planet(1, None, 0.048, 0.5, B=1, orbit_length=1000, star = M_star)
    model, x = SPI.SPI_Cubic(0.01, 0, 0.5, 0.055, 5, length = 500)
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), 
                        output_dir = M_output, obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.5_B_1.csv', planet_epoch= epoch)
        print(hosts)


def Run_Kappa_Sim():
    lc_dir = 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/G_Type_LC/219101992/tess2019253231442-s0016-0000000219101992-0152-s_lc.fits'
    lc_dir_M= 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/M_Type_LC/219822564/tess2019198215352-s0014-0000000219822564-0150-s_lc.fits'

    kappa_list = [0, 0.25, 0.5, 1, 2, 4, 8]
    lc_count = 25
    
    ### M_Type
    for kappa in kappa_list:
        for hosts in range(50):
            loc = 0.5
            for iterations in range(lc_count):
                lc, epoch = SPI.SPI_kappa_flare_injection(lc_dir_M, kappa, loc, 5, sp_type = 'M', flare_type='Flaring')
                flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
                Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = "M_" + str(219101992) + '_' + str(hosts) + '_kappa_' + str(kappa), output_dir = 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/Output/VM/M_Type', obs_time = lc_count*24*60*27, catalog_name = 'Kappa_Sim_M.csv', planet_epoch= epoch)
                
                print(iterations, hosts, kappa)
    ### G_Type
    for kappa in kappa_list:
        for hosts in range(50):
            loc = 0.5
            for iterations in range(lc_count):
                lc, epoch = SPI.SPI_kappa_flare_injection(lc_dir, kappa, loc, 5, sp_type = 'G', flare_type='Flaring')
                flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)

                Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = "G_" +str(219101992) +'_kappa_' + str(hosts) +'_' + str(kappa), output_dir = 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/Output/VM/G_Type', obs_time = lc_count*24*60*27, catalog_name = 'Kappa_Sim_G.csv', planet_epoch= epoch)
                print(iterations, hosts, kappa)


Run_Inverse_Cubic_Sim()
=======
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 13:45:06 2024

@author: Nate Whitsett
"""
import ardor.Flares.Flare as Flare
import ardor.SPI_Forward_Models.SPI_Simulation as SPI
import ardor.SPI_Forward_Models.Orbit_Model_Library as OML
import os
import numpy as np

def Run_Inverse_Cubic_Sim():
    ############# G TYPE SIMULATION #############
    base_dirG = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/G_Type_LC'
    base_dirM = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/M_Type_LC'
    targets = os.listdir(base_dirG)

    G_star = OML.Star(0.896, 1, 1, radius = 0.92165, age = 4.6e9, B = 1, alfven=0.1)
    M_star = OML.Star(0.59, 1, 1, radius = 0.603, age = 4.6e9, B = 100, alfven=0.1)
    lc_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/G_Type_LC/219101992/tess2019253231442-s0016-0000000219101992-0152-s_lc.fits'
    lc_dir_M= 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/M_Type_LC/219822564/tess2019198215352-s0014-0000000219822564-0150-s_lc.fits'

    lc_count = 25
    lc = Flare.tier0(lc_dir)
    planet_G = OML.Planet(1, None, 0.055, 0, B=10, orbit_length=1000, Star = G_star)
    model, x = OML.probability_density(G_star, planet_G, 0.5)
    ### e = 0.0, G_Type
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_G = OML.Planet(1, None, 0.055, 0.05, B=10, orbit_length=1000, Star = G_star)
    model, x = OML.probability_density(G_star, planet_G, 0.5)
    ### e = 0.05
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.05.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_G = OML.Planet(1, None, 0.055, 0.2, B=10, orbit_length=1000, Star = G_star)
    model, x = OML.probability_density(G_star, planet_G, 0.5)
    ### e = 0.2
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.2.csv', planet_epoch= epoch)
        print(hosts)
    
    
    ### e = 0.5
    planet_G = OML.Planet(1, None, 0.055, 0.5, B=10, orbit_length=1000, Star = G_star)
    model, x = OML.probability_density(G_star, planet_G, 0.5)
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.5.csv', planet_epoch= epoch)
        print(hosts)
    
        
    ############# M TYPE SIMULATION #############
    
    planet_M = OML.Planet(1, None, 0.048, 0, B=10, orbit_length=1000, Star = M_star)
    model, x = OML.probability_density(M_star, planet_M, 0.5)
    ### e = 0.0, M_Type
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_M = OML.Planet(1, None, 0.048, 0.05, B=10, orbit_length=1000, Star = M_star)
    model, x = OML.probability_density(M_star, planet_M, 0.5)
    ### e = 0.05
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.05.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_M = OML.Planet(1, None, 0.048, 0.2, B=10, orbit_length=1000, Star = M_star)
    model, x = OML.probability_density(M_star, planet_M, 0.5)
    ### e = 0.2
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.2.csv', planet_epoch= epoch)
        print(hosts)
    
    
    ### e = 0.5
    planet_M = OML.Planet(1, None, 0.048, 0.5, B=10, orbit_length=1000, Star = M_star)
    model, x = OML.probability_density(M_star, planet_M, 0.5)
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.5.csv', planet_epoch= epoch)
        print(hosts)
    
    
    
    
    
    
    
    
    
    
    ############# G TYPE SIMULATION SATURATED DYNAMO #############
    
    
    planet_G = OML.Planet(1, None, 0.055, 0, B=100, orbit_length=1000, Star = G_star)
    model, x = OML.probability_density(G_star, planet_G, 0.5)
    ### e = 0.0, G_Type
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0_B_100.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_G = OML.Planet(1, None, 0.055, 0.05, B=100, orbit_length=1000, Star = G_star)
    model, x = OML.probability_density(G_star, planet_G, 0.5)
    ### e = 0.05
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.05_B_100.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_G = OML.Planet(1, None, 0.055, 0.2, B=100, orbit_length=1000, Star = G_star)
    model, x = OML.probability_density(G_star, planet_G, 0.5)
    ### e = 0.2
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.2_B_100.csv', planet_epoch= epoch)
        print(hosts)
    
    
    ### e = 0.5
    planet_G = OML.Planet(1, None, 0.055, 0.5, B=100, orbit_length=1000, Star = G_star)
    model, x = OML.probability_density(G_star, planet_G, 0.5)
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.5_B_100.csv', planet_epoch= epoch)
        print(hosts)
    
        
    ############# M TYPE SIMULATION SATURATED DYNAMO #############
    
    planet_M = OML.Planet(1, None, 0.048, 0, B=100, orbit_length=1000, Star = M_star)
    model, x = OML.probability_density(M_star, planet_M, 0.5)
    ### e = 0.0, M_Type
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0_B_100.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_M = OML.Planet(1, None, 0.048, 0.05, B=100, orbit_length=1000, Star = M_star)
    model, x = OML.probability_density(M_star, planet_M, 0.5)
    ### e = 0.05
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.05_B_100.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_M = OML.Planet(1, None, 0.048, 0.2, B=100, orbit_length=1000, Star = M_star)
    model, x = OML.probability_density(M_star, planet_M, 0.5)
    ### e = 0.2
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.2_B_100.csv', planet_epoch= epoch)
        print(hosts)
    
    
    ### e = 0.5
    planet_M = OML.Planet(1, None, 0.048, 0.5, B=100, orbit_length=1000, Star = M_star)
    model, x = OML.probability_density(M_star, planet_M, 0.5)
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.5_B_100.csv', planet_epoch= epoch)
        print(hosts)
    
    
    
    
    
    ############# G TYPE SIMULATION WEAK FIELD #############
    
    
    planet_G = OML.Planet(1, None, 0.055, 0, B=1, orbit_length=1000, Star = G_star)
    model, x = OML.probability_density(G_star, planet_G, 0.5)
    ### e = 0.0, G_Type
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0_B_1.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_G = OML.Planet(1, None, 0.055, 0.05, B=1, orbit_length=1000, Star = G_star)
    model, x = OML.probability_density(G_star, planet_G, 0.5)
    ### e = 0.05
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.05_B_1.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_G = OML.Planet(1, None, 0.055, 0.2, B=1, orbit_length=1000, Star = G_star)
    model, x = OML.probability_density(G_star, planet_G, 0.5)
    ### e = 0.2
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.2_B_1.csv', planet_epoch= epoch)
        print(hosts)
    
    
    ### e = 0.5
    planet_G = OML.Planet(1, None, 0.055, 0.5, B=1, orbit_length=1000, Star = G_star)
    model, x = OML.probability_density(G_star, planet_G, 0.5)
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir, G_star, planet_G, fast = False, sp_type = 'G', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219101992) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'G_Type_e_0.5_B_1.csv', planet_epoch= epoch)
        print(hosts)
    
        
    ############# M TYPE SIMULATION WEAK FIELD #############
    
    planet_M = OML.Planet(1, None, 0.048, 0, B=1, orbit_length=1000, Star = M_star)
    model, x = OML.probability_density(M_star, planet_M, 0.5)
    ### e = 0.0, M_Type
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0_B_1.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_M = OML.Planet(1, None, 0.048, 0.05, B=1, orbit_length=1000, Star = M_star)
    model, x = OML.probability_density(M_star, planet_M, 0.5)
    ### e = 0.05
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.05_B_1.csv', planet_epoch= epoch)
        print(hosts)
        
    planet_M = OML.Planet(1, None, 0.048, 0.2, B=1, orbit_length=1000, Star = M_star)
    model, x = OML.probability_density(M_star, planet_M, 0.5)
    ### e = 0.2
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.2_B_1.csv', planet_epoch= epoch)
        print(hosts)
    
    
    ### e = 0.5
    planet_M = OML.Planet(1, None, 0.048, 0.5, B=1, orbit_length=1000, Star = M_star)
    model, x = OML.probability_density(M_star, planet_M, 0.5)
    for hosts in range(50):
        for iterations in range(lc_count):
            lc, epoch = SPI.SPI_cubic_flare_injection(lc_dir_M, M_star, planet_M, fast = False, sp_type = 'M', prior_model = True, phases = x, model = model)
            flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
            Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = str(219822564) + '_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data', obs_time = lc_count*24*60*27, catalog_name = 'M_Type_e_0.5_B_1.csv', planet_epoch= epoch)
        print(hosts)


def Run_Kappa_Sim():
    lc_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/G_Type_LC/219101992/tess2019253231442-s0016-0000000219101992-0152-s_lc.fits'
    lc_dir_M= 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/M_Type_LC/219822564/tess2019198215352-s0014-0000000219822564-0150-s_lc.fits'

    kappa_list = [0, 0.25, 0.5, 1, 2, 4, 8]
    lc_count = 25
    ### M_Type
    for kappa in kappa_list:
        for hosts in range(50):
            loc = 0.5
            for iterations in range(lc_count):
                lc, epoch = SPI.SPI_kappa_flare_injection(lc_dir_M, kappa, loc, 5, sp_type = 'M', flare_type='Flaring')
                flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)
                Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = "M_" + str(219101992) + '_' + str(hosts) + '_kappa_' + str(kappa), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Python Scripts/ardor/src/ardor/SPI_Forward_Models/Sim_Outputs_Kappa/M_Type', obs_time = lc_count*24*60*27, catalog_name = 'Kappa_Sim_M.csv', planet_epoch= epoch)
                print(iterations, hosts, kappa)
    ### G_Type
    for kappa in kappa_list:
        for hosts in range(50):
            loc = 0.5
            for iterations in range(lc_count):
                lc, epoch = SPI.SPI_kappa_flare_injection(lc_dir, kappa, loc, 5, sp_type = 'G', flare_type='Flaring')
                flares, lengths = Flare.tier1(lc.flux, 3, fast = False, injection = True)

                Flare.tier2(lc.time, lc.flux, lc.error, flares, lengths,chi_square_cutoff=20,  csv=False, host_name = "G_" +str(219101992) +'_kappa_' + str(kappa) +'_' + str(hosts), output_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Python Scripts/ardor/src/ardor/SPI_Forward_Models/Sim_Outputs_Kappa/G_Type', obs_time = lc_count*24*60*27, catalog_name = 'Kappa_Sim_G.csv', planet_epoch= epoch)
                print(iterations, hosts, kappa)
Run_Kappa_Sim()
>>>>>>> 843bcfd59f29e0297b65b05360ddcf3ed4f7ad77
