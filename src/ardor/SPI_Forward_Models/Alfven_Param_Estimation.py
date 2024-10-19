# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:26:48 2024

@author: whitsett.n
"""

# planet_list = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/Alfven_Catalog.csv')
# masses = planet_list['st_mass']
# # masses_err = (planet_list['st_masserr1'], planet_list['st_masserr2'])
# radii = planet_list['st_rad']
# # radii_err = (planet_list['st_raderr1'], planet_list['st_raderr2'])
# ages = planet_list['st_age']
# # ages_err = radii_err = (planet_list['st_ageerr1'], planet_list['st_ageerr2'])
# lumin = planet_list['st_lum']
# # lumin_err = (planet_list['st_lumerr1'], planet_list['st_lumerr2'])
# rot = planet_list['st_rotp']
# # rot_err = (planet_list['st_rotperr1'], planet_list['st_rotperr2'])


# alfven_rad = []
# alfven_lower = []
# alfven_upper = []

# B_list = []
# B_l = []
# B_u = []

# Mdot_list = []
# Mdot_l = []
# Mdot_u = []

# energy100_list = []
# energy100_l = []
# energy100_u = []
# energy25_list = []
# energy25_l = []
# energy25_u = []
# energy10_list = []
# energy10_l = []
# energy10_u = []

# time100_list = []
# time100_l = []
# time100_u = []
# time25_list = []
# time25_l = []
# time25_u = []
# time10_list = []
# time10_l = []
# time10_u = []

# lower_phase = []
# upper_phase = []



# a = planet_list['pl_orbsmax']
# # a_err = (planet_list['pl_orbsmaxerr1'], planet_list['pl_orbsmaxerr2'])
# e = planet_list['pl_orbeccen']
# # e_err = (planet_list['pl_orbeccenerr1'], planet_list['pl_orbeccenerr2'])
# per = planet_list['pl_orbper']
# # per_err = (planet_list['pl_orbpererr1'], planet_list['pl_orbpererr2'])
# R_pl = planet_list['pl_radj']
# # R_pl_err = (planet_list['pl_radjerr1'], planet_list['pl_radjerr2'])
# name = planet_list['Host_ID']
# spectral = planet_list['st_spectype']
# count = 0
# for index, hosts in enumerate(name):
#     print(hosts)
    
    
#     alfven_list = []
#     Lower_Phase_List = []
#     Upper_Phase_List = []
#     B_list1 = []
#     Mdot_list1 = []
    
    
#     time_list100_1 = []
#     energy_list100_1 = []
#     time_list25_1 = []
#     energy_list25_1 = []
#     time_list10_1 = []
#     energy_list10_1 = []
    
#     counter = 0 
#     if radii[index] > 3:
#         print('Giant Star')
#         upper_phase.append(np.nan)
#         lower_phase.append(np.nan)
#         continue
#     host = Star(masses[index], 1, lumin[index], radius=radii[index], age=1e9*ages[index], p_rot=rot[index], err = False)   
#     planet100 = Planet(R_pl[index], per[index], a[index], e[index], 100, masses[index], orbit_resolution=0.05)
#     planet25 = Planet(R_pl[index], per[index], a[index], e[index], 25, masses[index], orbit_resolution=0.05)
#     planet10 = Planet(R_pl[index], per[index], a[index], e[index], 10, masses[index], orbit_resolution=0.05)
    
#     # time100_1, energy100_1 = interaction(host, planet100)
#     # time25_1, energy25_1 = interaction(host, planet25)
#     # time10_1, energy10_1 = interaction(host, planet10)


#     alfven1 = host.Alfven
#     B1 = host.B
#     Mdot1 = host.massloss
    
#     if np.isnan(host.Alfven) == True:
#         B_list.append(np.nan)
#         B_l.append(np.nan)
#         B_u.append(np.nan)
#         Mdot_list.append(np.nan)
#         Mdot_l.append(np.nan)
#         Mdot_u.append(np.nan)
#         alfven_rad.append(np.nan)
#         alfven_lower.append(np.nan)
#         alfven_upper.append(np.nan)
#         lower_phase.append(np.nan)
#         upper_phase.append(np.nan)
#         energy100_list.append(np.nan)
#         energy100_l.append(np.nan)
#         energy100_u.append(np.nan)
#         energy25_list.append(np.nan)
#         energy25_l.append(np.nan)
#         energy25_u.append(np.nan)
#         energy10_list.append(np.nan)
#         energy10_l.append(np.nan)
#         energy10_u.append(np.nan)
#         time100_list.append(np.nan)
#         time100_l.append(np.nan)
#         time100_u.append(np.nan)
#         time25_list.append(np.nan)
#         time25_l.append(np.nan)
#         time25_u.append(np.nan)
#         time10_list.append(np.nan)
#         time10_l.append(np.nan)
#         time10_u.append(np.nan)
#         continue
    
#     # for trials in range(5000):
#     #     mass = asymmetric_error(masses[index], masses_err[0][index], masses_err[1][index])
#     #     radius = asymmetric_error(radii[index], radii_err[0][index], radii_err[1][index])
#     #     age = asymmetric_error(ages[index], ages_err[0][index], ages_err[1][index])
#     #     luminosity = asymmetric_error(lumin[index], lumin_err[0][index], lumin_err[1][index], lumin_check = True)
#     #     rotation = asymmetric_error(rot[index], rot_err[0][index], rot_err[1][index])
        
#     #     pl_a = asymmetric_error(a[index], a_err[0][index], a_err[1][index])
#     #     pl_e = asymmetric_error(e[index], e_err[0][index], e_err[1][index])
#     #     pl_rad = asymmetric_error(R_pl[index], R_pl_err[0][index], R_pl_err[1][index])
        
        
#     #     planet100_1 = Planet(pl_rad, per[index], pl_a, pl_e, 100, mass, orbit_resolution=1)
#     #     planet25_1 = Planet(pl_rad, per[index], pl_a, pl_e, 25, mass, orbit_resolution=1)
#     #     planet10_1 = Planet(pl_rad, per[index], pl_a, pl_e, 10, mass, orbit_resolution=1)
#     #     host = Star(mass, 1, luminosity, radius=radius, age=1e9*age, p_rot=rotation, err = True)
        
        
#     #     alfven = host.Alfven
#     #     B = host.B
#     #     Mdot = host.massloss
#     #     time100, energy100 = interaction(host, planet100)
#     #     time25, energy25 = interaction(host, planet25)
#     #     time10, energy10 = interaction(host, planet10)
#     #     if np.isnan(alfven) == True:
#     #         if len(alfven_list) > 5:
#     #             alfven = np.random.normal(loc = alfven1, scale = np.std(alfven_list))
#     #         else:
#     #             alfven = alfven1
#     #     if np.isnan(B) == True:
#     #         if len(B_list1) > 5:
#     #             B = np.random.normal(loc = B1, scale = np.std(B_list))
#     #         else:
#     #             B = B1
#     #     if np.isnan(Mdot) == True:
#     #         if len(Mdot_list1) > 5:
#     #             Mdot = np.random.normal(loc = Mdot1, scale = np.std(Mdot_list1))
#     #         else:
#     #             Mdot = Mdot1
#     #     if np.isnan(time100) == True:
#     #         if len(time_list100_1) > 5:
#     #             time100 = np.random.normal(loc = time100_1, scale = np.std(time_list100_1))
#     #         else:
#     #             time100 = time100_1
#     #     if np.isnan(time25) == True:
#     #         if len(time_list25_1) > 5:
#     #             time25 = np.random.normal(loc = time25_1, scale = np.std(time_list25_1))
#     #         else:
#     #             time25 = time25_1
#     #     if np.isnan(time10) == True:
#     #         if len(time_list10_1) > 5:
#     #             time10 = np.random.normal(loc = time10_1, scale = np.std(time_list25_1))
#     #         else:
#     #             time10 = time10_1
#     #     if np.isnan(energy100) == True:
#     #         if len(energy_list100_1) > 5:
#     #             energy100 = np.random.normal(loc = energy100_1, scale = np.std(energy_list100_1))
#     #         else:
#     #             energy100 = energy100_1
#     #     if np.isnan(energy25) == True:
#     #         if len(energy_list25_1) > 5:
#     #             energy25 = np.random.normal(loc = energy25_1, scale = np.std(energy_list25_1))
#     #         else:
#     #             energy25 = energy25_1
#     #     if np.isnan(energy10) == True:
#     #         if len(energy_list10_1) > 5:
#     #             energy10 = np.random.normal(loc = energy10_1, scale = np.std(energy_list10_1))
#     #         else:
#     #             energy10 = energy10_1
        
#     #     counter += 1
#     #     alfven_list.append(alfven)
#     #     B_list1.append(B)
#     #     Mdot_list1.append(Mdot)
#     #     time_list100_1.append(time100)
#     #     energy_list100_1.append(energy100)
#     #     time_list25_1.append(time25)
#     #     energy_list25_1.append(energy25)
#     #     time_list10_1.append(time10)
#     #     energy_list10_1.append(energy10)
#     # alfven_list = np.array(alfven_list)[~np.isnan(np.array(alfven_list))]
#     alfven = alfven1
#     if np.isnan(alfven1) == False:
#         period = planet100.period
#         orbit = planet100.orbit
#         time = planet100.time
#         diff = orbit - alfven
#         if diff[0] > 0:
#             super_alfvenic = True
#             sub_alfvenic = False
#         elif diff[0] < 0:
#             super_alfvenic = False
#             sub_alfvenic = True
#         check = False
#         count = 0
#         for indices, dist in enumerate(diff):
#             if super_alfvenic == False and sub_alfvenic == True and dist > 0:
#                 lower_phase.append(time[indices])
#                 super_alfvenic = True
#                 sub_alfvenic = False
#                 check = True
#             if super_alfvenic == True and sub_alfvenic == False and dist < 0:
#                 upper_phase.append(time[indices])
#                 sub_alfvenic = True
#                 super_alfvenic = False
#                 check = True
#             if dist < 0:
#                 count += 1
#         if count == len(diff):
#             lower_phase.append(0)
#             upper_phase.append(1)
#         elif check == False:
#             lower_phase.append(np.nan)
#             upper_phase.append(np.nan)
#     print(len(upper_phase), len(lower_phase))
    # alfven_list.sort()
    # B_list1.sort()
    # Mdot_list1.sort()
    # energy_list100_1.sort()
    # energy_list10_1.sort()
    # energy_list25_1.sort()
    # time_list100_1.sort()
    # time_list10_1.sort()
    # time_list25_1.sort()
    
    # alfven_rad.append(alfven1)
    # alfven_lower.append(alfven_list[1700])
    # alfven_upper.append(alfven_list[3400])
    # B_list.append(B1)
    # B_l.append(B_list1[1700])
    # B_u.append(B_list1[3400])
    # Mdot_list.append(Mdot1)
    # Mdot_l.append(Mdot_list1[1700])
    # Mdot_u.append(Mdot_list1[3400])
    # energy100_list.append(energy100_1)
    # energy100_l.append(energy_list100_1[1700])
    # energy100_u.append(energy_list100_1[3400])
    # energy25_list.append(energy25_1)
    # energy25_l.append(energy_list25_1[1700])
    # energy25_u.append(energy_list25_1[3400])
    # energy10_list.append(energy10_1)
    # energy10_l.append(energy_list10_1[1700])
    # energy10_u.append(energy_list10_1[3400])

    # time100_list.append(time100_1)
    # time100_l.append(time_list100_1[1700])
    # time100_u.append(time_list100_1[3400])
    # time25_list.append(time25_1)
    # time25_l.append(time_list25_1[1700])
    # time25_u.append(time_list25_1[3400])
    # time10_list.append(time10_1)
    # time10_l.append(time_list10_1[1700])
    # time10_u.append(time_list10_1[3400])
    # print(len(upper_phase), len(lower_phase))
    # print(len(upper_phase), len(lower_phase), len(alfven_rad), len(B_list), len(B_l), len(B_u), len(Mdot_list), len(Mdot_l), len(Mdot_u))
    # print(len(energy100_list), len(energy10_list), len(energy25_list), len(time100_list), len(time10_list), len(time25_list))
    # ZZ = np.column_stack((hosts, alfven1, alfven_list[1700], alfven_list[3400],B1,B_list1[1700],B_list1[3400],Mdot1,Mdot_list1[1700],Mdot_list1[3400], energy100_1, energy_list100_1[1700], energy_list100_1[3400], energy25_1, energy_list25_1[1700], energy_list25_1[3400], energy10_1, energy_list10_1[1700], energy_list10_1[3400], time100_1, time_list100_1[1700], time_list100_1[3400], time25_1, time_list25_1[1700], time_list25_1[3400], time10_1, time_list10_1[1700], time_list10_1[3400]))
    # with open("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Alfven_Catalog.csv", "a") as f:
    # ZZ = np.column_stack((lower_phase, upper_phase))
        #     np.savetxt(f, ZZ, delimiter=",", fmt='%s')
    #     f.close()
# planet_list['Alfven_Rad'] = alfven_rad
# planet_list['Alfven_Rad_Lower'] = alfven_lower
# planet_list['Alfven_Rad_Upper'] = alfven_upper
# planet_list['Sub_Alfv_lphase'] = lower_phase
# planet_list['Sub_Alfv_uphase'] = upper_phase
# planet_list['B_st'] = B_list
# planet_list['B_st_err1'] = B_l
# planet_list['B_st_err2'] = B_u
# planet_list['Mdot'] = Mdot_list
# planet_list['Mdot_err1'] = Mdot_l
# planet_list['Mdot_err2'] = Mdot_u

# planet_list['Interaction_Energy_100'] = energy100_list
# planet_list['Interaction_Energy_100_l'] = energy100_l
# planet_list['Interaction_Energy_100_u'] = energy100_u

# planet_list['Interaction_Energy_10'] = energy10_list
# planet_list['Interaction_Energy_10_l'] = energy10_l
# planet_list['Interaction_Energy_10_u'] = energy10_u

# planet_list['Interaction_Energy_25'] = energy25_list
# planet_list['Interaction_Energy_25_l'] = energy25_l
# planet_list['Interaction_Energy_25_u'] = energy25_u

# planet_list['Interaction_time_100'] = time100_list
# planet_list['Interaction_time_100_l'] = time100_l
# planet_list['Interaction_time_100_u'] = time100_u

# planet_list['Interaction_time_10'] = time10_list
# planet_list['Interaction_time_10_l'] = time10_l
# planet_list['Interaction_time_10_u'] = time10_u

# planet_list['Interaction_time_25'] = time25_list
# planet_list['Interaction_time_25_l'] = time25_l
# planet_list['Interaction_time_25_u'] = time25_u

# planet_list.to_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/Alfven_Catalog_2.csv', index=False)
# for index, row in planet_list.iterrows():
#     if math.isnan(row['pl_orbsmax']) == True:
#         a = 0.164
#     elif math.isnan(row['pl_orbsmax']) == False:
#         a = row['pl_orbsmax']
#     if math.isnan(row['st_lum']) == True:
#         b = 1
#     elif math.isnan(row['st_lum']) == False:
#         b = row['st_lum']
#     if math.isnan(row['sy_dist']) == True:
#         c = 50
#     elif math.isnan(row['sy_dist']) == False:
#         c = row['sy_dist']

#     star = Star(row['st_mass'], c, b, age=row['st_age']*1e9)
#     planet = Planet(row['pl_radj'], a, row['pl_orbeccen'], 100, arg_periastron=row['pl_orblper'], orbit_resolution=0.0002, inclination=row['pl_orbincl'], stellar_mass = star.mass)
#     if interaction_checker(star, planet) == True:
#         interact, energy, power = interaction(star, planet)
#         curve, index2 = phase_curve(star, planet, interact)
#         probability = probability_density(star, planet, index2)
#         phase = np.arange(0, 1, 1/len(probability))
#         phase = phase[0:31416]
#         time = planet.time
#         parameters = pd.DataFrame(index=range(len(planet.orbit)), columns=[0])
#         parameters.loc[0, 0]= 'semi-major axis: ' + str(planet.a)
#         parameters.loc[1, 0]= 'eccentricity: ' + str(planet.e)
#         parameters.loc[2, 0]= 'True Anomaly: ' + str((planet.true_anomaly/np.pi)*180)
#         parameters.loc[3, 0]= 'Alfven Radius: ' + str(star.Alfven)
#         parameters.loc[4, 0] = 'Periastron_phase: ' + str(planet.periastron-1)
#         parameters.loc[5, 0] = 'Periastron_time: ' + str(planet.time[probability.index(max(probability))])
#         parameters = parameters.values.tolist()
#         export_curve = pd.DataFrame({'phase': phase, 'time': time, 'curve': curve, 'distance': planet.orbit2, 'relative probability': probability, 'parameters': parameters, 'flare_energy': energy, 'flare_power': power})
#         export_curve.to_csv('Desktop/' + row[0]+ '.csv')
