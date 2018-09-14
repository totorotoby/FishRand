# TODO for monday or tuesday...fix zooplankton steady state, and nonsteady state first thing...then try to do fish differential equation

# Where we run the bioaccumlation model
import matplotlib
matplotlib.use("TkAgg")
import Classes as obj
import FR_Input_Output
import prob as pr
from copy import deepcopy

def set_all_h_and_s(model_para, all_data):

    for l in all_data:
        for item in l:
            for i in range(len(item)):
                if type(item[i]) == pr.Var:
                    pr.set_hyper_samp_cube(model_para, item[i])
                if type(item[i]) == list:
                    list_item = item[i]
                    for j in range(len(list_item)):
                        if (item[i][j]) == pr.Var:
                            pr.set_hyper_samp_cube(model_para, item[i][j])


def check_inst_non_st(inst, u_count, v_count):
    u_count = int(u_count)
    v_count = int(v_count)
    new_entry = []
    for j in range(len(inst)):
        if type(inst[j]) == pr.Var:
            if inst[j].type == 'U':
                new_entry.append(inst[j].values[u_count])
            if inst[j].type == 'V':
                new_entry.append(inst[j].values[v_count])
        else:
            new_entry.append(inst[j])

    return new_entry


# initiates all regions
def init_region(reg_data, temp, u_count, v_count):

    regions = []
    for i in range(len(reg_data)):
        region = reg_data[i]
        region = check_inst_non_st(region, u_count, v_count)
        toadd = obj.Region(region[0], temp, region[1], region[2], region[3], region[4], region[5])

        # dealing with cox
        if region[6] != '':
            toadd.set_cox(region[7])
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
            exit(0)

    return regions

# initiates chemicals for a specific single region


def init_chems(chem_data, region, region_index, u_count, v_count):

    chemicals = []
    for i in range(len(chem_data)):

        chemical = chem_data[i]
        chemical = check_inst_non_st(chemical, u_count, v_count)
        toadd = obj.Chemical(chemical[0], chemical[1], chemical[2][region_index])

        # dealing with ddoc and dpoc
        if chemical[5] and chemical[6] != '':
            toadd.set_ddoc_dpoc(chemical[5], chemical[6])

        # dealing with cwto and cwdo
        if chemical[3][i] != '':
            toadd.set_cwto(chemical[3][region_index])
        if chemical[4][i] != '':
            toadd.set_cwdo(chemical[4][region_index])

    # calculating Phi
        toadd.calc_phi_and_cwdo(region)
        toadd.calc_pore_water(region.Ocs)

        chemicals.append(toadd)

    return chemicals


def init_phyto(phyto_data, chemicals, u_count, v_count, per_step=0):

    phytos = []
    phyto = phyto_data[0]
    phyto = check_inst_non_st(phyto, u_count, v_count)

    # Potential loop in the future

    toadd = obj.Pplank(phyto[0], phyto[1], len(chemicals), per_step)

    if phyto[2] != '':
        toadd.set_vlp(phyto[2])
    if phyto[3] != '':
        toadd.set_vnp(phyto[3])

    toadd.calc_vwb()

    for i in range(len(chemicals)):

        kow = chemicals[i].Kow
        toadd.calc_k1(kow, i)
        toadd.calc_k2(kow, i)

    if toadd.init_check():
        phytos.append(toadd)
    else:
        print("Something is wrong with your ", 1,  "st Phyto Entry")
        exit(0)

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
            toadd.calc_k1(ew, j)
            toadd.calc_k2(kow, j)
            toadd.calc_kgb(kow, j)
            toadd.calc_ke(ed, j)
            toadd.calc_kd(ed, j)

        zoops.append(toadd)

    return zoops


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

    # TODO Need to redo once spatial is ready specifically regions[0] will not necessarily be regions 0
    for i in range(len(fishs)):

        fishs[i].calc_diet_per(tempfishs, region)
        fishs[i].calc_gut_per()
        fishs[i].calc_gf()

        for j in range(len(chemicals)):
            kow = chemicals[j].Kow
            ed = chemicals[j].Ed
            ew = chemicals[j].Ew
            fishs[i].calc_k1(ew, j)
            fishs[i].calc_k2(kow, j)
            fishs[i].calc_kgb(kow, j)
            fishs[i].calc_ke(ed, j)
            fishs[i].calc_kd(ed, j)

    return fishs


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


def solve_steady(region, chemicals, phytos, zoops, fishs):

    # nested dictonary where can look up first by region then chemical then animal to find concentration
    conc_log = {}
    conc_log[region.name] = {}

    phytolog = conc_log[region.name]

    # assuming one phyto
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
            Cwto = chemicals[i].Cwto
            phyto_con = conc_log[region.name]['Phytoplankton'][chemicals[i].name]

            conc_in_zoops = zoo.solve_steady_state(phi, i, Cwto, Cwdo, phyto_con)
            Cb.append(conc_in_zoops)
            zooplog[zoo.name][chemicals[i].name] = conc_in_zoops

        zoo.Cb = Cb

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


def init_prior_con_dic(regions, chems, phyto, zoops, fishs):

    prior_con_dic = {}
    for k in range(len(regions)):
        prior_con_dic[regions[k][0]] = {}
        for i in range(len(fishs)):
            prior_con_dic[regions[k][0]][fishs[i][0]] = {}
            for j in range(len(chems)):
                prior_con_dic[regions[k][0]][fishs[i][0]][chems[j][0]] = 0
        for i in range(len(zoops)):
            prior_con_dic[regions[k][0]][zoops[i][0]] = {}
            for j in range(len(chems)):
                prior_con_dic[regions[k][0]][zoops[i][0]][chems[j][0]] = 0
        prior_con_dic[regions[k][0]][phyto[0][0]] = {}
        for j in range(len(chems)):
            prior_con_dic[regions[k][0]][phyto[0][0]][chems[j][0]] = 0
    return prior_con_dic


def solve_zoop_phyto_time_period(region, chemicals, phytos, zoops, fishs, prior_chem_amounts, t):

    new_chem_amounts = deepcopy(prior_chem_amounts)

    phyto = phytos[0]
    for j in range(len(chemicals)):

        # Phytos and zoops
        prior_C = prior_chem_amounts[region.name][phyto.name][chemicals[j].name]
        new_chem_amounts[region.name][phyto.name][chemicals[j].name] = phyto.solve_next_time_step(chemicals[j].Cwdo, j,
                                                                                                  prior_C, t)
        for k in range(len(zoops)):
            prior_C = prior_chem_amounts[region.name][zoops[k].name][chemicals[j].name]
            phi = chemicals[j].phi
            Cwdo = chemicals[j].Cwdo
            Cwto = chemicals[j].Cwto

            phyto_con = (new_chem_amounts[region.name][phyto.name][chemicals[j].name] +
                         prior_chem_amounts[region.name][phyto.name][chemicals[j].name])/2
            new_chem_amounts[region.name][zoops[k].name][chemicals[j].name] =\
                zoops[k].solve_next_time_step(phi, j, Cwto, Cwdo, phyto_con, prior_C)

        # Fishes!

        fishlog_new = new_chem_amounts[region.name]
        fishlog_old = prior_chem_amounts[region.name]
        for i in range(len(fishs)):
            prior_C = prior_chem_amounts[region.name][fishs[i].name][chemicals[j].name]
            phi = chemicals[j].phi
            Cwdo = chemicals[j].Cwdo
            Cwp = chemicals[i].Cwp
            new_chem_amounts[region.name][fishs[i].name][chemicals[j].name] =\
                fishs[i].solve_next_time_step(phi, j, Cwp, Cwdo, fishlog_new, fishlog_old, chemicals[j], prior_C)

    # print(new_chem_amounts['The Pond']['Pumpkinseed PC']['PCB52'], ', ')
    return new_chem_amounts


def bio_monte_carlo_loop(model_para, all_data, temp, reg_index, u_iter, v_iter, t, per_step=0, p_dic=None):

    p_dic_list = 0
    # check if case is temporal and monte carlo
    if type(p_dic) is list:
        p_dic_list = 1
    else:
        prior_conc = p_dic

    dictionaries = []
    set_all_h_and_s(model_para, all_data)

    inner_count = 0
    u_count = 0
    while u_count < u_iter:
        u_count += 1
        v_count = 0
        while v_count < v_iter:

            if p_dic_list == 1:
                prior_conc = p_dic[inner_count]

            concentration_log = single_bio_deterministic_iter(model_para, all_data, temp, reg_index, u_count,
                                                              inner_count, t, per_step=per_step, p_dic=prior_conc)
            v_count += 1
            inner_count += 1
            dictionaries.append(concentration_log)

    return dictionaries


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

        print(dictionaries[0])
        return dictionaries[0]

    # if we ran bio monte carlo
    else:
        results_dic = pr.make_result_dist(dictionaries)
        result_print(results_dic)
        return results_dic


def single_bio_deterministic_iter(model_para, all_data, temp, reg_index, u_count, v_count, t, per_step=0, p_dic=None):

    reg_data = all_data[0]
    chem_data = all_data[2]

    regions = init_region(reg_data, temp, u_count, v_count)
    chemicals = init_chems(chem_data, regions[reg_index], reg_index, u_count, v_count)
    phytos = init_phyto(all_data[3], chemicals, u_count, v_count, per_step=per_step)
    zoops = init_zoop(all_data[4], regions[reg_index], chemicals, phytos[0], u_count, v_count, per_step=per_step)
    fishs, tempfishs = init_fish_pre_region(all_data[5], regions[reg_index], chemicals, phytos[0], zoops, all_data[6],
                                            u_count, v_count, per_step=per_step)
    fishs = reorder_fish(fishs)
    fishs = init_fish_post_region(fishs, tempfishs, regions[reg_index], chemicals)

    if len(all_data[0]) == 1 and model_para[8] == 'YES':

        conc_log = solve_steady(regions[0], chemicals, phytos, zoops, fishs)
        return conc_log

    if model_para[8] == 'NO':

        new_concentrations = solve_zoop_phyto_time_period(regions[reg_index], chemicals, phytos, zoops, fishs, p_dic, t)
        return new_concentrations

