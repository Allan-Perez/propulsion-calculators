import os
import itertools as itt


def linspace(start, end, num):
  l = []
  delta = (end - start)/num
  c = start
  for i in range(num):
    l.append(c)
    c+=delta
  return l
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
      to_sweep_arr.append( [(name, el) for el in linspace(*el)] )
    else:
      to_sweep_arr.append( [(name, el[0])] )
  cartesian_sweep = list(itt.product(*to_sweep_arr))
  # dataframe = pd.DataFrame(None, columns= input_parameters + output_parameters)
  i = 0
  total = len(cartesian_sweep)
  increment = total/40
  for sim_element in cartesian_sweep:
    # print(f"{i}/{t}")
    print("[" + "=" * int(i / increment) +  " " * int((total - i)/ increment) + "]" +  str((i/total)*100) + "%", end='\n' if i+1 == total else '\r')
    # dataframe = main(sim_element, dataframe)
    sim_element = tuples_to_dict(sim_element)
    command = f"./simulation.o {sim_element['initial_tank_pressure']} {sim_element['tank_volume']} {sim_element['orifice_diameter']} {sim_element['orifice_number']} {sim_element['chamber_pressure_bar']} {sim_element['initial_ullage']} {sim_element['orifice_k2_coefficient']} {i}"
    os.system(command)
    i+=1
  # return dataframe
  

if __name__=="__main__":
  options_nums = range(4)
  options_ref  = (57.2265, 0.0084995122110964, 0.0015113, 22, 34.4738, 85, 30.0)
  options_default_min_max_delta = 0.1
  options_default_n = 5
  options_names = ("initial_tank_pressure", "tank_volume", "orifice_diameter", "orifice_number", "chamber_pressure_bar", "initial_ullage", "orifice_k2_coefficient")
  #Optimization Variables
  pcc = options_ref[0] #Combustion Chamber Pressure (psia)
  pamb = options_ref[1] #Ambient Atmospheric Pressure (psia)
  mr = options_ref[2] #oxidizer to fuel mass ratio
  mass_flow = options_ref[3] #mass flow rate (lbm/s)

  print(f"""
    0: initial_tank_pressure [{options_names[0]}](bar)
    1: tank_volume = 0.0084995122110964 [{options_names[1]}](m3)
    2: orifice_diameter = 0.0015113 [{options_names[2]}](m)
    3: orifice_number = 22 [{options_names[0]}]
    4: chamber_pressure_bar = 34.4738 [{options_names[0]}](bar)
    5: initial_ullage = 85 [{options_names[0]}]
    6: orifice_k2_coefficient = [{options_names[0]}]
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
  sweep(to_sweep, not_to_sweep)
  #print(res)
  # res.to_csv(f'{date.today().strftime("%Y%m%d")}_sweep_results.csv')

  #main(pcc, pamb, mr, mass_flow)
