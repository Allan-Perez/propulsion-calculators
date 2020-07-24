#!/usr/bin/env python3
from rocketcea.cea_obj import CEA_Obj
from math import exp, pi, sqrt, tan, radians
from numpy import log
import numpy as np
#import psweep as ps
import itertools as itt
import pandas as pd
from datetime import date
#import RP1_prop as rp_prop
# Rocket CEA Docs: https://rocketcea.readthedocs.io/en/latest/
#Created by Akshay Kulkarni - 02/09/2019
#This script is intended to determine nozzle and combustion chamber dimensions based on specific inputs.
#Please beware that the parameters are based on ideal assumptions and equilibrium calculations based on NASA-CEA outputs.

#Constants
R = 1545.32 #ft-lbf/lbmol-R
g = 32.2 #ft/s^2

input_parameters = ('Chamb_P', 'Optimum_Expansion_Ratio', 'Mass_flow', 'Mass_Ratio', 'Ambient_P')
output_parameters = ('Chamb_T', 'R_Exit', 'Throat_A', 'Throat_Ratio', 'C*', 'Exit_Mach', 'Isp', 'Chamb_Cp', \
                      'Chamb_Molec_Weight', 'Chamb_Density', 'Exit_V', 'Thrust')

def main(params, df_carry):
  pcc = params['pcc']
  pamb = params['pamb']
  mr = params['mr']
  mass_flow = params['mass_flow']
  #CEA Variables
  ispObj = CEA_Obj( oxName='N2O', fuelName='HTPB')
  eps = ispObj.get_eps_at_PcOvPe(pcc, mr, pcc/pamb); #Optimim Expansion Ratio
  isp = ispObj.estimate_Ambient_Isp(pcc, mr, eps, pamb)

  #Chamber Properties
  c_properties = ispObj.get_Chamber_MolWt_gamma(pcc, mr, eps) #tuple of chamber properties
  c_molweight = c_properties[0] #Molecular Weight - in lb/lbmol
  c_gamma = c_properties[1] #Ratio of Specific Heats
  c_temp = ispObj.get_Tcomb(pcc, mr) #Chamber Combustion Temperature (Rankine)
  c_star = ispObj.get_Cstar(pcc, mr)
  c_cp = ispObj.get_Chamber_Cp(pcc, mr, eps)
  c_density = ispObj.get_Chamber_Density(pcc, mr, eps)

  #Throat Properties
  t_prop = ispObj.get_IvacCstrTc_ThtMwGam(pcc, mr, eps) #Gas Properties in the Throat
  t_gamma = t_prop[4]
  t_temp = (2*c_temp)/(t_gamma+1) #Eqn 3-22 from RPE, Using Throat Gamma (Rankine)
  t_pressureratio = ispObj.get_Throat_PcOvPe(pcc, mr) #Chamber to Throat Pressure Ratio
  t_pressure = pcc/t_pressureratio #Pressure at the throat (psi) - verified with RPE pg.57 - Pc/Pt = 0.56

  #Augmented Constants
  R_specific_comb = (R/c_molweight)
  R_specific_comb_in = (R/(12*c_molweight)) #Universal Gas Constant for Chamber

  #Nozzle Parameters
  #Throat
  t_area = (mass_flow/pcc)*sqrt((c_temp*R_specific_comb_in)/c_gamma)*(1+(((c_gamma-1)/2)**((c_gamma+1)/(2*(c_gamma-1))))); #Throat Area (in^2) - refer: https://www.grc.nasa.gov/www/k-12/airplane/rktthsum.html
  #Assuming Chamber Temperature and Chamber Pressure are Stagnation Pressures - velocity in CC assumed to be negligible
  t_radius = sqrt(t_area/pi) #Throat Radius (in)
  t_diameter = 2*2.54*t_radius #Throat Diameter (cm)
  t_area_cm = t_area*6.4516 #Throat Area (cm^2)

  #Exit
  e_pressureratio = ispObj.get_PcOvPe(pcc, mr, eps) #Ratio of Pressure from Combustion Chamber to Exit
  e_pressure_actual = pcc/e_pressureratio #Actual Exit Pressure (psi)
  e_area = eps*t_area #Nozzle Exit Area (in^2)
  e_prop = ispObj.get_IvacCstrTc_exitMwGam(pcc, mr, eps) #Nozzle Exit Properties
  e_molweight = e_prop[3] #Molecular Weight - in lb/lbmol
  e_gamma = e_prop[4] #Nozzle Exit Gamma - this parameter changes significantly from the chamber gamma
  e_mach = ispObj.get_MachNumber(pcc, mr, eps) #Exit Mach Number

  #Augmented Constants
  R_specific_exit = (R/e_molweight) #Universal Gas Constant for Nozzle Exit

  #Thrust Calcs
  e_vel_HH = sqrt(((2*g*e_gamma*R_specific_exit*c_temp)/(e_gamma-1))*(1- ((pamb)/(pcc))**((e_gamma-1)/(e_gamma)))) #H&H Eqn 1.18
  thrust_HH = (mass_flow*(e_vel_HH))/(32.2) + (e_pressure_actual - pamb)*e_area

  """
  print ("Input Parameters")
  print ("Chamber Pressure (psi) = ", pcc)
  print ("Expansion Ratio = ", eps)
  print ("Mass Flow Rate (lb/s) = ", mass_flow)
  print ("Mass Ratio = ", mr)
  print ("Ambient Pressure (psi) = ", pamb)

  print ("Calculated Parameters")
  print ("Chamber Temperature (R) = ", c_temp)
  print ("R for Exit Gases (ft-lbf/lbm R) = ", R_specific_exit)
  print ( "Throat Area (in^2) = ", t_area )
  print ("Throat Radius (in) = ", t_radius)
  print ('C* (ft/s) = ', c_star)
  print ("Exit Mach Number = ", e_mach)
  print ("Isp = ", isp)
  print ("Chamber Cp = ", c_cp)
  print ("Chamber Molecular Weight = ", c_molweight)
  print ("Chamber Density (lbm/ft^3) = ", c_density)

  print ("Outputs")
  print ("Exit Velocity ft(s) = ", e_vel_HH)
  print ("Thrust (lbf)", thrust_HH)
  """
  input_results = [pcc, eps, mass_flow, mr, pamb]
  output_results = [c_temp, R_specific_exit, t_area, t_radius, \
                  c_star, e_mach, isp, c_cp, c_molweight, \
                  c_density, e_vel_HH, thrust_HH]
  new_series = pd.Series(input_results + output_results, index= input_parameters + output_parameters)
  df_carry = df_carry.append(new_series, ignore_index=True)
  return df_carry

def tuples_to_dict(tuple_):
  dict_ = {}
  for el in tuple_:
    dict_[el[0]] = el[1]
  return dict_

def sweep(to_sweep, not_to_sweep):
  # one sweep at a time. 
  # use recurrently
  # itt.product(*to_sweep)
  to_sweep_arr = []
  for el in to_sweep + not_to_sweep:
    name = el[0]
    del el[0]
    if len(el) > 1: 
      to_sweep_arr.append( [(name, el) for el in np.linspace(*el)] )
    else:
      to_sweep_arr.append( [(name, el[0])] )
  cartesian_sweep = list(itt.product(*to_sweep_arr))
  dataframe = pd.DataFrame(None, columns= input_parameters + output_parameters)
  for sim_element in cartesian_sweep:
    sim_element = tuples_to_dict(sim_element)
    print(f"Sim element: {tuples_to_dict(sim_element)}")
    dataframe = main(sim_element, dataframe)
  return dataframe
  

if __name__=="__main__":
  options_nums = range(4)
  options_ref  = (500, 14.7, 6.3, 3.63)
  options_default_min_max_delta = 0.1
  options_default_n = 5
  options_names = ("pcc", "pamb", "mr", "mass_flow")
  #Optimization Variables
  pcc = options_ref[0] #Combustion Chamber Pressure (psia)
  pamb = options_ref[1] #Ambient Atmospheric Pressure (psia)
  mr = options_ref[2] #oxidizer to fuel mass ratio
  mass_flow = options_ref[3] #mass flow rate (lbm/s)

  print(f"""
    0: Combusion Chamber Pressure [{options_names[0]}](psi abs -- psi relative to vacumm) 
    1: Ambient Atmospheric Pressure [{options_names[1]}](psi abs -- psi relative to vacumm) 
    2: Oxidizer to fuel mass ratio [{options_names[2]}](no units) 
    3: Mass flow rate [{options_names[3]}](lbm/s) 
    """)

  sweep_params_options = input("""
    Please write the numbers of the options that you wish 
    to sweep simulate (comma-separated, e.g. 0,1 or 2,3,4): \n""")
  sweep_params_options = sweep_params_options.split(",")
  to_sweep = []
  for option_str in sweep_params_options:
    option_num = int(option_str.strip())
    assert option_num in options_nums, f"{option_num} is an invalid option."
    tmp_ref_d = options_ref[option_num] * options_default_min_max_delta
    tmp_min_ref = options_ref[option_num] - tmp_ref_d
    tmp_max_ref = options_ref[option_num] + tmp_ref_d
    option_min = float(input(f"Minimum value of {options_names[option_num]} to simulate (e.g. {tmp_min_ref}): "))
    option_max = float(input(f"Maximum value of {options_names[option_num]} to simulate (e.g. {tmp_max_ref}): "))
    option_n = int(input(f"Number of {options_names[option_num]} simulations (e.g. {options_default_n}): "))
    to_sweep.append([options_names[option_num], option_min, option_max, option_n])
  to_sweep_names = set([el[0] for el in to_sweep])
  not_to_sweep_names = set(options_names) - to_sweep_names
  not_to_sweep = []
  for option_str in not_to_sweep_names:
    idx_option = options_names.index(option_str)
    tmp = float(input(f"Default value for {option_str} (e.g. {options_ref[idx_option]}): "))
    not_to_sweep.append([option_str, tmp])

  #print("To sweep: ", to_sweep)
  #print("Not to sweep: ", not_to_sweep)
  res = sweep(to_sweep, not_to_sweep)
  print(res)
  res.to_csv(f'{date.today().strftime("%Y%m%d")}_sweep_results.csv')

  #main(pcc, pamb, mr, mass_flow)
