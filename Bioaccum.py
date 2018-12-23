import matplotlib
import numpy as np
matplotlib.use("TkAgg")
import Classes as obj
import prob as pr
from copy import deepcopy
import spatial


# once all the model parameters taken in and are lists, this function loops through them, and
# does all the sampling for each statistical value
def set_all_h_and_s(model_para, all_data):

    count = 0
    for l in all_data:
        for item in l:
            for i in range(len(item)):
                if type(item[i]) == pr.Var:
                    count += 1
                    pr.set_hyper_samp_cube(model_para, item[i])
                if type(item[i]) == list:
                    list_item = item[i]
                    for j in range(len(list_item)):
                        if type(list_item[j]) == pr.Var:
                            count += 1
                            pr.set_hyper_samp_cube(model_para, list_item[j])
                        if type(list_item[j]) == list:
                            nested_list_item = list_item[j]
                            for k in range(len(nested_list_item)):
                                print(type(nested_list_item[k]))
                                if type(nested_list_item[k]) == pr.Var:
                                    count += 1
                                    pr.set_hyper_samp_cube(model_para, nested_list_item[k])
                                    print(nested_list_item[k].lhs)
                                          

    if count >= 1:
        return True
    return False

# Goes into a statistical parameter and sets it as a single value according the the u_count or v_count index of the sample array
def check_inst_non_st(inst, u_count, v_count):
    u_count = int(u_count)
    v_count = int(v_count)
    new_entry = []
    for i in range(len(inst)):
        if type(inst[i]) is pr.Var:
            if inst[i].type == 'U':
                print(inst[i].values, len(inst[i].values))
                new_entry.append(inst[i].values[u_count])
            if inst[i].type == 'V':
                new_entry.append(inst[i].values[v_count])
        elif type(inst[i]) == list:
            new_sub_entry = []
            for j in range(len(inst[i])):
                if type(inst[i][j]) == str or type(inst[i][j]) == float:
                    new_sub_entry.append(inst[i][j])
                if type(inst[i][j]) == pr.Var:
                    if inst[i][j].type == 'U':
                        new_sub_entry.append(inst[i][j].values[u_count])
                    if inst[i][j].type == 'V':
                        new_sub_entry.append(inst[i][j].values[v_count])

            new_entry.append(new_sub_entry)

        else:

            new_entry.append(inst[i])

    return new_entry


# initiates all regions for a single monte carlo sample in a single time period
def init_region(reg_data, temp, u_count, v_count):

    regions = []
    for i in range(len(reg_data)):
        region = reg_data[i]
        region = check_inst_non_st(region, u_count, v_count)
        toadd = obj.Region(region[0], temp[i], region[1], region[2], region[3], region[4], region[5])

        # dealing with cox
        if region[6] != '':
            toadd.set_cox(region[6])
        else:
            toadd.calc_cox()

        # dealing with adoc apoc
        if region[7] != '' and region[8] != '':
            toadd.set_adoc_apoc(region[7], region[8])

        # checking to see if anything is missing
        if toadd.init_check():
            regions.append(toadd)
        else:
            print("Something is wrong with your ", i, " Regional Entry")
    return regions

# initiates chemicals for a specific single region in a single monte carlo sample in a time period
def init_chems(chem_data, r_con_data, region, region_index, u_count, v_count):
    
    r_concentrations = r_con_data[region_index]
    
    chemicals = []
    for i in range(len(chem_data)):
        
        concentrations = r_concentrations[i]
        chemical = chem_data[i]

        chemical = check_inst_non_st(chemical, u_count, v_count)
        concentrations = check_inst_non_st(concentrations, u_count, v_count)
        toadd = obj.Chemical(chemical[0], chemical[1], concentrations[0])

        # dealing with ddoc and dpoc
        if chemical[2] and chemical[3] != '':
            toadd.set_ddoc_dpoc(chemical[2], chemical[3])
        if chemical[4] != '':
            toadd.set_beta1(chemical[4])
        if chemical[5] != '':
            toadd.set_beta2(chemical[5])
        if chemical[6] != '':
            toadd.set_beta3(chemical[6])
        if chemical[7] != '':
            toadd.set_beta4(chemical[7])
        if chemical[8] != '':
            toadd.set_beta5(chemical[8])
        
        # dealing with cwto and cwdo
        if concentrations[1] != '':
            toadd.set_cwto(concentrations[1])
        if concentrations[2] != '':
            toadd.set_cwdo(concentrations[2])
        if concentrations[3] != '':
            toadd.set_cwp(concentrations[3])

        if toadd.Cwp != -1 and toadd.Cwdo != -1 and toadd.Cs == '' and region.Ocs != '':
            toadd.calc_cs(region)
            
        if toadd.Cwp != -1 and toadd.Cwdo == -1 and toadd.Cs == '' and toadd.Cwto != -1:
            toadd.calc_cs(region)
            
        if toadd.Cwdo == -1 and toadd.Cwto != -1 and region.adoc and region.apoc and region.Xdoc and region.Xpoc != '':
            toadd.calc_phi_and_cwdo(region)
            
        if toadd.Cwdo == -1 and (toadd.Cwto == -1 or region.adoc or region.apoc or region.Xdoc or region.Xpoc == ''):
            print('Not enough parameters to solve for Cwdo in ' + chemical[0])
            exit(0)
        if toadd.Cwp == -1 and region.Ocs != '':
            toadd.calc_pore_water(region.Ocs)
        if toadd.Cwp == -1 and (region.Ocs == '' or toadd.Cs == ''):
            print('Not enough parameters to solve for Cwp in' + chemical[0])
            exit(0)
            
        chemicals.append(toadd)


    
    return chemicals


# same as chemicals and regions but for phytos
def init_phyto(phyto_data, chemicals, u_count, v_count, per_step=0):

    phytos = []
    phyto = phyto_data[0]
    phyto = check_inst_non_st(phyto, u_count, v_count)

    toadd = obj.Pplank(phyto[0], phyto[1], len(chemicals), per_step)

    if phyto[2] != '':
        toadd.set_vlp(phyto[2])
    if phyto[3] != '':
        toadd.set_vnp(phyto[3])

    toadd.calc_vwb()
    for i in range(len(chemicals)):

        kow = chemicals[i].Kow
        beta1 = chemicals[i].beta1
        beta2 = chemicals[i].beta2
        beta4 = chemicals[i].beta4
        toadd.calc_k1(kow, i)
        toadd.calc_k2(kow, i, beta1, beta4)

    if toadd.init_check():
        phytos.append(toadd)
    else:
        print("Something is wrong with your ", 1,  "st Phyto Entry")


    return phytos


def init_zoop(zoo_data, region, chemicals, phyto, u_count, v_count, per_step=0):

    zoops = []

    for i in range(len(zoo_data)):
        zoop = zoo_data[i]
        zoop = check_inst_non_st(zoop, u_count, v_count)
        # Potential loop in the future
        toadd = obj.Zooplank(zoop[0], zoop[1], zoop[2], zoop[10], len(chemicals), per_step)

        # inital set
        if zoop[3] != '':
            toadd.set_vnb(zoop[3])
        if zoop[4] != '':
            toadd.set_mp(zoop[4])
        if zoop[5] != '' and zoop[6] != '' and zoop[7] != '':
            toadd.set_el_en_ew(zoop[5], zoop[6], zoop[7])
        if zoop[8] != '':
            toadd.set_gd(zoop[8])
        if zoop[9] != '':
            toadd.set_kg(zoop[9])

        toadd.calc_mo()
        toadd.calc_vwb()
        toadd.calc_diet_per(phyto)
        toadd.calc_gut_per()

        toadd.calc_gv(region.Cox)
        if toadd.gd_set == 0:
            if toadd.flag == 0:
                toadd.calc_gd_no_filter(region.T)
            if toadd.flag == 1:
                toadd.calc_gd_filter(region.Css)
            if toadd.kg_set == 0:
                toadd.calc_kg(region.T)
        toadd.calc_gf()

        for j in range(len(chemicals)):

            kow = chemicals[j].Kow
            ed = chemicals[j].Ed
            ew = chemicals[j].Ew
            beta1 = chemicals[j].beta1
            beta2 = chemicals[j].beta2
            beta3 = chemicals[j].beta3
            beta4 = chemicals[j].beta4
            beta5 = chemicals[j].beta5

            toadd.calc_k1(ew, j)
            toadd.calc_k2(kow, j, beta1, beta5)
            toadd.calc_kgb(kow, j, beta3, beta4 ,beta5)
            toadd.calc_ke(ed, j)
            toadd.calc_kd(ed, j)

        zoops.append(toadd)


    return zoops


def init_single_fish_pre_region(fish_data, region, chemicals, diet, u_count, v_count, per_step=0):

    fish = check_inst_non_st(fish_data, u_count, v_count)
    toadd = obj.Fish(fish[0], fish[1], fish[2], diet, fish[10], len(chemicals), per_step)
    tempadd = obj.Fish(fish[0], fish[1], fish[2], diet, fish[10], len(chemicals), per_step)

    if fish[3] != '':
        toadd.set_vnb(fish[3])
        tempadd.set_vnb(fish[3])

    tempadd.calc_vwb()
    toadd.calc_vwb()

    # inital param
    if fish[6] != '' and fish[7] != '' and fish[8] != '':
        toadd.set_el_en_ew(fish[6], fish[7], fish[8])
    if fish[8] != '':
        toadd.set_gd(fish[8])
    if fish[9] != '':
        toadd.set_kg(fish[9])
    if fish[4] != '':
        toadd.set_mp(fish[4])

    toadd.calc_mo()

    toadd.calc_gv(region.Cox)
    if toadd.gd_set == 0:
        if toadd.flag == 0:
            toadd.calc_gd_no_filter(region.T)
        if toadd.flag == 1:
            toadd.calc_gd_filter(region.Css)
        if toadd.kg_set == 0:
            toadd.calc_kg(region.T)


    return toadd, tempadd


def init_fish_pre_region(fish_data, region, chemicals, phyto, zoops, diet_data, u_count, v_count, per_step=0):

    tempfishs = []
    for zoop in zoops:
        tempfishs.append(zoop)
    tempfishs.append(phyto)

    fishs = []
    for fish in fish_data:

        fish = check_inst_non_st(fish, u_count, v_count)
        # setting up diets
        diet = diet_data[fish[0]]
        toadd = obj.Fish(fish[0], fish[1], fish[2], diet, fish[10], len(chemicals), per_step)
        tempadd = obj.Fish(fish[0], fish[1], fish[2], diet, fish[10], len(chemicals), per_step)

        if fish[3] != '':
            toadd.set_vnb(fish[3])
            tempadd.set_vnb(fish[3])

        tempadd.calc_vwb()
        toadd.calc_vwb()
        tempfishs.append(tempadd)

        # inital param
        if fish[6] != '' and fish[7] != '' and fish[8] != '':
            toadd.set_el_en_ew(fish[6], fish[7], fish[8])
        if fish[8] != '':
            toadd.set_gd(fish[8])
        if fish[9] != '':
            toadd.set_kg(fish[9])
        if fish[4] != '':
            toadd.set_mp(fish[4])

        toadd.calc_mo()

        toadd.calc_gv(region.Cox)
        if toadd.gd_set == 0:
            if toadd.flag == 0:
                toadd.calc_gd_no_filter(region.T)
            if toadd.flag == 1:
                toadd.calc_gd_filter(region.Css)
            if toadd.kg_set == 0:
                toadd.calc_kg(region.T)

        fishs.append(toadd)

    return fishs, tempfishs

# pre and post are necessary to deal with the diet issue


def init_fish_post_region(fishs, tempfishs, region, chemicals):


    for i in range(len(fishs)):

        fishs[i].calc_diet_per(tempfishs, region)
        fishs[i].calc_gut_per()
        fishs[i].calc_gf()

        for j in range(len(chemicals)):
            kow = chemicals[j].Kow
            ed = chemicals[j].Ed
            ew = chemicals[j].Ew
            beta1 = chemicals[j].beta1
            beta2 = chemicals[j].beta2
            beta3 = chemicals[j].beta3
            beta4 = chemicals[j].beta4
            beta5 = chemicals[j].beta5

            fishs[i].calc_k1(ew, j)
            fishs[i].calc_k2(kow, j, beta1, beta5)
            fishs[i].calc_kgb(kow, j, beta3, beta4, beta5)
            fishs[i].calc_ke(ed, j)
            fishs[i].calc_kd(ed, j)


    return fishs


def init_single_fish_post_region(fish, tempfishs, region, chemicals):

    fish.calc_diet_per(tempfishs, region)
    fish.calc_gut_per()
    fish.calc_gf()

    for j in range(len(chemicals)):
        kow = chemicals[j].Kow
        ed = chemicals[j].Ed
        ew = chemicals[j].Ew
        beta1 = chemicals[j].beta1
        beta2 = chemicals[j].beta2
        beta3 = chemicals[j].beta3
        beta4 = chemicals[j].beta4
        beta5 = chemicals[j].beta5
        fish.calc_k1(ew, j)
        fish.calc_k2(kow, j, beta1, beta5)
        fish.calc_kgb(kow, j, beta3, beta4, beta5)
        fish.calc_ke(ed, j)
        fish.calc_kd(ed, j)

    return fish


# Function to reorder fish *has not been tested*
def reorder_fish(fishs):

    new_order = []
    below = []
    # find the root fish
    while True:
        fish = find_base_fish(fishs)
        if type(fish) == obj.Fish:
            new_order.append(fish)
            below.append(fish.name)
            fishs.remove(fish)
        if fish is None:
            break

        below.append('Phytoplankton')
        below.append('Zooplankton')
        below.append('Sediment/Detritus')

    while True:
        nextfish = find_next_fish(fishs, below)
        if type(nextfish) == obj.Fish:
            new_order.append(nextfish)
            below.append(nextfish.name)
            fishs.remove(nextfish)
        if nextfish is None:
            break

    return new_order


def find_base_fish(fishs):

    basefish = None

    for fish in fishs:
        # counter for non 0 non phyto zoop diet
        count = 0
        for entry in fish.diet_frac:

            if (entry[0] != 'Sediment/Detritus' and entry[0] != 'Phytoplankton' and entry[0] != 'Zooplankton')\
                    and entry[1] > 0:
                count += 1

        if count == 0:
            basefish = fish

    return basefish


def find_next_fish(possible_fish, below):

    nextfish = None

    for fish in possible_fish:
        count = 0
        for entry in fish.diet_frac:
            if (entry[0] not in below) and entry[1] > 0:

                count += 1

        if count == 0:
            nextfish = fish

    return nextfish



def solve_steady(region, chemicals, phytos, zoops, inverts, fishs):

    # nested dictonary where can look up first by region then chemical then animal to find concentration
    conc_log = {}
    conc_log[region.name] = {}

    phytolog = conc_log[region.name]

    # assuming one phytoƒƒ
    phyto = phytos[0]
    phytolog[phyto.name] = {}
    for i in range(len(chemicals)):
        cwd = chemicals[i].Cwdo
        conc_in_phyto = phyto.solve_steady_state(cwd, i)
        phytolog[phyto.name][chemicals[i].name] = conc_in_phyto

    zooplog = conc_log[region.name]
    for zoo in zoops:
        Cb = []
        zooplog[zoo.name] = {}
        for i in range(len(chemicals)):
            # since single region all chemical properties index by 0
            phi = chemicals[i].phi
            Cwdo = chemicals[i].Cwdo
            Cwp = chemicals[i].Cwp
            phyto_con = conc_log[region.name]['Phytoplankton'][chemicals[i].name]

            conc_in_zoops = zoo.solve_steady_state(phi, i, Cwp, Cwdo, phyto_con)
            Cb.append(conc_in_zoops)
            zooplog[zoo.name][chemicals[i].name] = conc_in_zoops

        zoo.Cb = Cb
    invertlog = conc_log[region.name]

    for i in range(len(inverts)):
        invertlog[inverts[i].name] = {}
        Cb = []
        for j in range(len(chemicals)):
            phi = chemicals[j].phi
            Cwp = chemicals[j].Cwp
            Cwdo = chemicals[j].Cwdo
            con_in_i = inverts[i].solve_steady_state(phi, j, Cwp, Cwdo, invertlog, chemicals[j])
            Cb.append(con_in_i)
            invertlog[inverts[i].name][chemicals[j].name] = con_in_i
        inverts[i].Cb = Cb


    fishlog = conc_log[region.name]
    for i in range(len(fishs)):
        fishlog[fishs[i].name] = {}
        Cb = []
        for j in range(len(chemicals)):
            phi = chemicals[j].phi
            Cwp = chemicals[j].Cwp
            Cwdo = chemicals[j].Cwdo
            con_in_i = fishs[i].solve_steady_state(phi, j, Cwp, Cwdo, fishlog, chemicals[j])
            Cb.append(con_in_i)
            fishlog[fishs[i].name][chemicals[j].name] = con_in_i
        fishs[i].Cb = Cb

    return conc_log


def init_prior_con_dic(monte_carlo_length,regions, chems, phyto, zoops, inverts, fishs, draws):

    prior_con_dic_nonfish_list = []
    for p in range(monte_carlo_length):
        prior_con_dic_nonfish = {}
        for k in range(len(regions)):
            prior_con_dic_nonfish[regions[k][0]] = {}
            for i in range(len(zoops)):
                prior_con_dic_nonfish[regions[k][0]][zoops[i][0]] = {}
                for j in range(len(chems)):
                    prior_con_dic_nonfish[regions[k][0]][zoops[i][0]][chems[j][0]] = 0
            for s in range(len(inverts)):
                prior_con_dic_nonfish[regions[k][0]][inverts[s][0]] = {}
                for j in range(len(chems)):
                    prior_con_dic_nonfish[regions[k][0]][inverts[s][0]][chems[j][0]] = 0
            prior_con_dic_nonfish[regions[k][0]][phyto[0][0]] = {}
            for j in range(len(chems)):
                prior_con_dic_nonfish[regions[k][0]][phyto[0][0]][chems[j][0]] = 0


        prior_con_dic_nonfish_list.append(prior_con_dic_nonfish)

    prior_con_dic_fish_list = [[[[0 for _ in range(len(chems))] for _ in range(draws)] for _ in range(len(fishs))] for _ in range(monte_carlo_length)]

    return prior_con_dic_nonfish_list, prior_con_dic_fish_list


def solve_zoop_phyto_invert_time_period(region, chemicals, phytos, zoops, inverts, prior_chem_amounts, days):

    new_chem_amounts = deepcopy(prior_chem_amounts)

    phyto = phytos[0]
    for j in range(len(chemicals)):

        # Phytos and zoops
        prior_C = prior_chem_amounts[region.name][phyto.name][chemicals[j].name]
        new_chem_amounts[region.name][phyto.name][chemicals[j].name] = phyto.solve_next_time_step(chemicals[j].Cwdo, j,
                                                                                                  prior_C, days)
        for k in range(len(zoops)):
            prior_C = prior_chem_amounts[region.name][zoops[k].name][chemicals[j].name]
            phi = chemicals[j].phi
            Cwdo = chemicals[j].Cwdo
            Cwp = chemicals[j].Cwp

            phyto_con = (new_chem_amounts[region.name][phyto.name][chemicals[j].name] +
                         prior_chem_amounts[region.name][phyto.name][chemicals[j].name])/2

            new_chem_amounts[region.name][zoops[k].name][chemicals[j].name] =\
                zoops[k].solve_next_time_step(phi, j, Cwp, Cwdo, phyto_con, prior_C)


        fishlog_new = new_chem_amounts[region.name]
        fishlog_old = prior_chem_amounts[region.name]
        for i in range(len(inverts)):
            prior_C = prior_chem_amounts[region.name][inverts[i].name][chemicals[j].name]
            phi = chemicals[j].phi
            Cwdo = chemicals[j].Cwdo
            Cwp = chemicals[j].Cwp
            new_chem_amounts[region.name][inverts[i].name][chemicals[j].name] =\
                inverts[i].solve_next_time_step(phi, j, Cwp, Cwdo, fishlog_old, fishlog_new, chemicals[j], prior_C, 0) + prior_C


    return new_chem_amounts

def not_eating(predator):

    sum = 0
    for i in range(len(predator.diet_frac)):
        sum += predator.diet_frac[i][1]

    if sum == 0:
        return True

    return False

def eats_fish(predator, non_fish_num):
    count = 0

    for i in range(non_fish_num, len(predator.diet_frac)):
        count += predator.diet_frac[i][1] == 0

    if count > 0:
        return True

    return False


# when solving for a fish this function gets the concentrations in the prey of the fish
def get_prey_con(predator, reg_index, len_chems, chem_names, fish_names, non_fish_num, locations, conc_fishs):

    # There are no prior concentrations if fish doesn't eat
    if not_eating(predator) == True:
        return None

    if eats_fish(predator,non_fish_num) == False:
        return None

    # Otherwise find the populations in region previously and now and average there concentrations
    else:
        reg_check = 0
        prior_conc_prey = {}

        #loop over prey
        for i in range(len(predator.diet_frac)):
            prey_chemicals = []

            # if prey is being eaten, and is a fish
            if predator.diet_frac[i][1] != 0 and predator.diet_frac[i][0] in fish_names:

                # make prey section of dictionary
                prior_conc_prey[predator.diet_frac[i][0]] = {}

                fish_index = i - (non_fish_num)

                # get prior prey locations
                prior_prey_locations = locations[fish_index]

                #loop over locations
                for k in range(len(prior_prey_locations)):

                    # if location is this region
                    if prior_prey_locations[k] == reg_index:
                        reg_check = 1
                        old_concs = []
                        #get chemicals from that region
                        for j in range(len_chems):

                            # make entry in dictionary
                            prior_conc_prey[predator.diet_frac[i][0]][chem_names[j]] = 0

                            #get prior chemicals
                            old_concs.append(conc_fishs[fish_index][k][j])
                        prey_chemicals.append(old_concs)

            if len(prey_chemicals) != 0 :
                prey_chemicals = np.asarray(prey_chemicals)
                for j in range(len_chems):
                    average_chemical_con = np.mean(prey_chemicals[:, j])
                    prior_conc_prey[fish_names[fish_index]][chem_names[j]] = average_chemical_con

        if reg_check == 1:
            return prior_conc_prey
        else:
            return None


def solve_single_fish_single_region(region, chemicals, fish, days, prior_prey_cons, new_prey_cons, prior_C, out_check):

    new_cons_by_chemical = []
    for j in range(len(chemicals)):
        prior_c = prior_C[j]
        phi = chemicals[j].phi
        Cwdo = chemicals[j].Cwdo
        Cwp = chemicals[j].Cwp
        new_con = fish.solve_next_time_step(phi, j, Cwp, Cwdo, new_prey_cons, prior_prey_cons, chemicals[j], prior_c, out_check) + prior_c
        new_cons_by_chemical.append(new_con)
    return new_cons_by_chemical

def bio_monte_carlo_loop(model_para, all_data, t, time_per_step, fish_by_region, prior_locations, locations, u_iter, v_iter, p_dic=None):

    # renaming so the brain works
    region_data = all_data[0]
    temperatures = np.asarray(all_data[1])
    temperatures = temperatures[: , t]
    chem_data = all_data[2]
    phyto_data = all_data[3]
    zoop_data = all_data[4]
    invert_data = all_data[5]
    fish_data = all_data[6]
    diet_data = all_data[7]
    r_con_data = all_data[9][t]

    
    f_names = [fish[0] for fish in fish_data]
    # dictionaries will be appended to with the 8
    uv_results = [[],[]]
    steady_state = []
    # enter the monte carlo loops
    inner_count = 0
    u_count = 0
    while u_count < u_iter:

        u_count += 1
        v_count = 0

        while v_count < v_iter:


            if model_para[8] == 'NO':
                # if there are distrubtions then we have many simulations per timestep so we need to get the correct
                # prior concentrations dictionary for non_fish
                prior_conc_nonfish = p_dic[0][inner_count]
                new_concentrations_lower = deepcopy(prior_conc_nonfish)
                # single_setup_bottom_web returns:
                # - list of regions objects
                # - 2d list chemicals objects by all regions
                # - a phyto object
                # - 2d list of zoop objects by all regions
                # - 2d list of invert objects by all regions

                regions, r_chems, r_phytos, r_zoops, r_inverts = single_setup_bottom_web(region_data, temperatures, chem_data, r_con_data,
                                                                                         phyto_data, zoop_data, invert_data,
                                                                                         diet_data, u_count, v_count, time_per_step)

                ############################################################################################################

                # days that have passed used only for phyto cause we solve analytically
                days = (t+1) * time_per_step
                for i in range(len(regions)):
                    new_concentrations_lower = solve_zoop_phyto_invert_time_period(regions[i], r_chems[i], r_phytos[i], r_zoops[i], r_inverts[i],
                                                                             new_concentrations_lower, days)


                p_dic[0][inner_count] = new_concentrations_lower
                #################################### Bottom Feeds Done ######################################################

                r_fish = [[] for _ in range(len(regions))]
                r_tempfish = [[] for _ in range(len(regions))]
                # start making the tempfish so that diet can work
                for i in range(len(r_tempfish)):
                    r_tempfish[i].append(r_phytos[i][0])
                    for k in range(len(r_zoops[i])):
                        r_tempfish[i].append(r_zoops[i][k])
                    for p in range(len(r_inverts[i])):
                        r_tempfish[i].append(r_inverts[i][p])

                # loops over number of fish species... ie do all simulations for a single fish first
                for i in range(len(all_data[6])):
                    #loops over number of regions...do simulation for each region that fish is in
                    for j in range(len(regions)):

                        # The diet of the fish must be adjusted at this point to make sure that the fish doesn't
                        # try to eat anything that isn't in the region
                        # adj_fish_diet is the diet that is used for these iterations

                         adj_fish_diet = spatial.adjust_diet_to_region(fish_data[i][0], j, diet_data,
                                                                       fish_by_region, len(r_inverts[0]) + 2)

                         # initates a single fish in a single region, before we know the foodweb of that region
                         # returns tempfish so that we can by the end of this iteration know the food web, and solve for the
                         # rest of fish that is foodweb dependent

                         fish, tempfish = init_single_fish_pre_region(fish_data[i], regions[j], r_chems[j], adj_fish_diet
                                                                      , u_count, v_count, time_per_step)

                         r_fish[j].append(fish)
                         r_tempfish[j].append(tempfish)


                #now we can go back over each region, and solve things foodweb specific for that region
                for i in range(len(r_fish)):
                    for j in range(len(r_fish[i])):
                        r_fish[i][j] = init_single_fish_post_region(r_fish[i][j], r_tempfish[i], regions[i], r_chems[i])


                prior_conc_fishs = p_dic[1][inner_count]
                new_concentrations_lower = p_dic[0][inner_count]
                new_conc_fishs = deepcopy(prior_conc_fishs)
                for i in range(len(fish_data)):
                     for j in range(len(locations[i])):
                         reg_index = locations[i][j]
                         if reg_index != 1000:
                            if fish_data[i][0] in fish_by_region[reg_index]:
                                prior_conc_fish_prey = get_prey_con(r_fish[reg_index][i], reg_index, len(chem_data), [chem[0] for chem in chem_data],f_names, 3 + len(invert_data), prior_locations, prior_conc_fishs)
                                current_conc_fish_prey = get_prey_con(r_fish[reg_index][i], reg_index, len(chem_data), [chem[0] for chem in chem_data],f_names, 3 + len(invert_data), locations, new_conc_fishs)
                                new_concentrations = [new_concentrations_lower[region_data[reg_index][0]], current_conc_fish_prey]
                                prior_concentrations = [prior_conc_nonfish[region_data[reg_index][0]], prior_conc_fish_prey]
                                this_fish_prior_con = prior_conc_fishs[i][j]
                                calc_cons = solve_single_fish_single_region(regions[reg_index], r_chems[reg_index], r_fish[reg_index][i], time_per_step,
                                                             prior_concentrations, new_concentrations, this_fish_prior_con, 0)

                         else:
                             reg_index = np.random.randint(0, len(regions))
                             this_fish_prior_con = prior_conc_fishs[i][j]
                             calc_cons = solve_single_fish_single_region(regions[reg_index], r_chems[reg_index], r_fish[reg_index][i], time_per_step,
                                                             None, None, this_fish_prior_con, 1)

                         new_conc_fishs[i][j] = calc_cons


                uv_results[0].append(new_concentrations_lower)
                uv_results[1].append(new_conc_fishs)
                v_count += 1
                inner_count += 1

            if model_para[8] == 'YES':

                temp = all_data[1][0]
                regions = init_region(region_data, temp, u_count, v_count)
                chemicals = init_chems(chem_data, r_con_data, regions[0], 0, u_count, v_count)
                
                phytos = init_phyto(all_data[3], chemicals, u_count, v_count)
                zoops = init_zoop(all_data[4], regions[0], chemicals, phytos[0], u_count, v_count,)
                inverts, tempinverts = init_fish_pre_region(invert_data, regions[0], chemicals, phytos[0], zoops,
                                                            diet_data, u_count, v_count)

                inverts = init_fish_post_region(inverts, tempinverts, regions[0], chemicals)
                fishs, tempfishs = init_fish_pre_region(fish_data, regions[0], chemicals, phytos[0], zoops,
                                                        diet_data, u_count, v_count)
                tempfishs = tempinverts + tempfishs[2:]
                fishs = init_fish_post_region(fishs, tempfishs, regions[0], chemicals)
                conc_log = solve_steady(regions[0], chemicals, phytos, zoops, inverts, fishs)
                steady_state.append(conc_log)
                v_count += 1
                inner_count += 1

    if model_para[8] == 'YES':
        
        print(steady_state)
        return steady_state

    return uv_results


def pretty(d, indent=0):
    for key, value in d.items():
        print('\t' * indent + str(key))
        if isinstance(value, dict):
            pretty(value, indent+1)
        else:
            print('\t' * (indent+1) + str(value))


def result_print(results_dic):
    for region in results_dic.values():
        for animal, animals in region.items():
            for chemical, cont in animals.items():
                print(cont.bestparam()[0])
                cont.show()


def filter_stat_case(dictionaries):

    # if no statistical simulation
    if len(dictionaries) == 1:


        return dictionaries[0]

    # if we ran bio monte carlo
    else:
        results_dic = pr.make_result_dist(dictionaries)
        return results_dic


def single_setup_bottom_web(region_data, temperatures, chem_data, r_con_data, phyto_data, zoop_data, invert_data, diet_data, u_count, v_count, days_per_step):
    # return regions
    regions = init_region(region_data, temperatures, u_count, v_count)

    # Dimensions of r_chems: regions by chemical
    # Dimension of r_zoop and r_invert: regions by type of invert
    r_chems = []
    r_zoops = []
    r_phytos = []
    r_inverts = []


    # sets up all chemicals in all regions
    for i in range(len(regions)):
        print('hello')
        chemicals = init_chems(chem_data, r_con_data ,regions[i], i, u_count, v_count)
        r_chems.append(chemicals)
        print(r_chems)

    for i in range(len(regions)):
        phyto = init_phyto(phyto_data, r_chems[i], u_count, v_count, per_step=days_per_step)
        r_phytos.append(phyto)
        zoops = init_zoop(zoop_data, regions[i], r_chems[i], phyto[0], u_count, v_count, per_step=days_per_step)
        r_zoops.append(zoops)

        inverts, tempprey = init_fish_pre_region(invert_data, regions[i], r_chems[i], phyto[0], zoops, diet_data,
                                                 u_count, v_count, per_step=days_per_step)
        inverts = init_fish_post_region(inverts, tempprey, regions[i], r_chems[i])
        r_inverts.append(inverts)

    return regions, r_chems, r_phytos, r_zoops, r_inverts
