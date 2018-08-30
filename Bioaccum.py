# Where we run the bioaccumlation model
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import Classes as obj
import FR_Input_Output
import prob as pr

def set_all_h_and_s(model_para, all_data):

    for l in all_data:
        for item in l:
            for i in range(len(item)):
                if type(item[i]) == pr.Var:
                    pr.set_hyper_samp_cube(model_para, item[i])
                if type(item[i]) == list:
                    list_item = item[i]
                    for j in range (len(list_item)):
                        if (item[i][j]) == pr.Var:
                            pr.set_hyper_samp_cube(model_para,item[i][j])


def check_inst_non_st(inst,u_iter, v_iter):
    u_iter = int(u_iter)
    v_iter = int(v_iter)
    new_entry = []
    for j in range (len(inst)):
        if type(inst[j]) == pr.Var:
            if inst[j].type == 'U':
                new_entry.append(inst[j].values[u_iter])
            if inst[j].type == 'V':
                new_entry.append(inst[j].values[v_iter])
        else:
            new_entry.append(inst[j])

    return new_entry


# initiates all regions
def init_region(time_steps, reg_data, temp_data, v_iter, u_iter):

    regions = []
    for i in range(len(reg_data)):
        region = reg_data[i]
        region = check_inst_non_st(region,u_iter, v_iter)
        temps = check_inst_non_st(temp_data[i], u_iter, v_iter)
        toadd = obj.Region(region[0],temps,region[1],region[2],region[3],region[4],region[5], time_steps)

        # dealing with cox
        if region[6] != '':
            for i in range(time_steps):
                toadd.set_cox(region[7], i)
        else:
            for i in range (time_steps):
                toadd.calc_cox(i)

        # dealing with adoc apoc
        if region[7] and region[8] != '':
            toadd.set_adoc_apoc(region[7],region[8])

        # checking to see if anything is missing
        if toadd.init_check():
            regions.append(toadd)
        else:
            print("Something is wrong with your ", i, " Regional Entry")
            exit(0)

    return regions

# initiates chemicals for a specific single region
def init_chems_unconditional(chem_data, len_regions,u_iter, v_iter):

    chemicals = []
    for i in range(len(chem_data)):

        chemical = chem_data[i]
        chemical = check_inst_non_st(chemical,u_iter, v_iter)

        toadd = obj.Chemical(chemical[0],chemical[1],chemical[2], len_regions)

        # dealing with cwto and cwdo
        if chemical[3] != '':
            toadd.set_cwto(chemical[3])
        if chemical[4] != '':
            toadd.set_cwdo(chemical[4])

        # dealing with ddoc and dpoc
        if chemical[5] and chemical[6] != '':
            toadd.set_ddoc_dpoc(chemical[5],chemical[6])

        chemicals.append(toadd)

    return chemicals

def init_chems_conditional(chems, regions):

    count = 1
    for toadd in chems:
        for i in range(len(regions)):
        # calculating Phi
            toadd.calc_phi_and_cwdo(regions[i],  i)
            toadd.calc_pore_water(regions[i], i)

        # checking to see if anything is missing
        if toadd.init_check():
            continue
        else:
            print("Something is wrong with your ",count ," Chemical Entry")
            exit(0)
        count += 1

    return chems


def init_phyto(phyto_data, chemicals, u_iter, v_iter):

    phytos = []
    phyto = phyto_data[0]
    phyto = check_inst_non_st(phyto,u_iter, v_iter)

    # Potential loop in the future

    toadd = obj.Pplank(phyto[0])

    if phyto[1] != '':
        toadd.set_kg(phyto[1])
    if phyto[2] != '':
        toadd.set_vlp(phyto[2])
    if phyto[3] != '':
        toadd.set_vnp(phyto[3])

    toadd.calc_vwp()

    for chemical in chemicals:

        kow = chemical.Kow
        k_1 = toadd.calc_k1(kow)
        toadd.k_1.append(k_1)
        k_2 = toadd.calc_k2(kow,k_1)
        toadd.k_2.append(k_2)

    if toadd.init_check():
        phytos.append(toadd)
    else:
        print("Something is wrong with your ", 1,  "st Phyto Entry")
        exit(0)


    return phytos


def init_zoop(zoo_data, regions, chemicals, phyto, time_steps, u_iter, v_iter):

    zoops = []

    for i in range (len(zoo_data)):
        zoop = zoo_data[i]
        zoop = check_inst_non_st(zoop,u_iter, v_iter)
        # Potential loop in the future

        toadd = obj.Zooplank(zoop[0],zoop[1],zoop[2], zoop[10], time_steps, len(regions), len(chemicals))

        # inital set
        if zoop[3] != '':
            toadd.set_vnb(zoop[3])
        if zoop[4] != '':
            toadd.set_mp(zoop[4])
        if zoop[5] != '' and zoop[6] != '' and zoop[7] != '':
            toadd.set_el_en_ew(zoop[5],zoop[6],zoop[7])
        if zoop[8] != '':
            toadd.set_gd(zoop[8], time_steps, len(regions))
        if zoop[9] != '':
            toadd.set_kg(zoop[9], time_steps, len(regions))

        toadd.calc_mo()
        toadd.calc_vwb()
        toadd.calc_diet_per(phyto)
        toadd.calc_gut_per()

        for j in range (len(regions)):
            for k in range (time_steps):
                toadd.calc_gv(regions[j].Cox, k, j)
                if toadd.gd_set == 0:
                    if toadd.flag == 0:
                        toadd.calc_gd_no_filter(regions[j].T, k, j)
                    if toadd.flag == 1:
                        toadd.calc_gd_filter(regions[j].Css, k, j)
                    if toadd.kg_set == 0:
                        toadd.calc_kg(regions[j].T, k, j)
                    toadd.calc_gf(k, j)

                for p in range(len(chemicals)):

                    kow = chemicals[p].Kow
                    ed = chemicals[p].Ed
                    ew = chemicals[p].Ew
                    toadd.calc_k1(ew, k, j, p)
                    toadd.calc_k2(kow, toadd.k_1, k, j ,p)
                    toadd.calc_kgb(kow, p)
                    toadd.calc_ke(ed, k, j, p)
                    toadd.calc_kd(ed, k, j ,p)

        zoops.append(toadd)

    return zoops



def init_fish_pre_region(fish_data, regions, chemicals, phyto, zoops, diet_data, time_steps, u_iter, v_iter):

    tempfishs = []
    for zoop in zoops:
        tempfishs.append(zoop)
    tempfishs.append(phyto)

    fishs = []

    for fish in fish_data:

        fish = check_inst_non_st(fish,u_iter, v_iter)

        # setting up diets
        diet = diet_data[fish[0]]

        toadd = obj.Fish(fish[0],fish[1],fish[2],fish[4], diet, fish[10], time_steps, len(regions), len(chemicals))
        tempadd = obj.Fish(fish[0],fish[1],fish[2],fish[4], diet, fish[10], time_steps, len(regions), len(chemicals))

        if fish[3] != '':
            toadd.set_vnb(fish[3])
            tempadd.set_vnb(fish[3])

        tempadd.calc_vwb()
        toadd.calc_vwb()
        tempfishs.append(tempadd)


        #inital param
        if fish[6] != '' and fish[7] != '' and fish[8] != '':
            toadd.set_el_en_ew(fish[6],fish[7],fish[8])
        if fish[8] != '':
            toadd.set_gd(fish[8], time_steps, len(regions))
        if fish[9] != '':
            toadd.set_kg(fish[9], time_steps, len(regions))

        toadd.calc_mo()

        for j in range (len(regions)):
            for k in range (time_steps):
                toadd.calc_gv(regions[j].Cox, k, j)
                if toadd.gd_set == 0:
                    if toadd.flag == 0:
                        toadd.calc_gd_no_filter(regions[j].T, k, j)
                    if toadd.flag == 1:
                        toadd.calc_gd_filter(regions[j].Css, k, j)
                    if toadd.kg_set == 0:
                        toadd.calc_kg(regions[j].T, k, j)
        #print(toadd.Vnd)
        fishs.append(toadd)

    return fishs, tempfishs


def init_fish_post_region(fishs, tempfishs, regions, chemicals, time_steps):

    # TODO Need to redo once spatial is ready specifically regions[0] will not necessarily be regions 0
    for i in range(len(fishs)):

        fishs[i].calc_diet_per(tempfishs, regions[0].Ocs)
        fishs[i].calc_gut_per()


        # ^^^ From here on is Fine only ^^^
        for j in range (len(regions)):
            for k in range (time_steps):
                fishs[i].calc_gf(k,j)


                for p in range (len(chemicals)):
                    kow = chemicals[p].Kow
                    ed = chemicals[p].Ed
                    ew = chemicals[p].Ew
                    fishs[i].calc_k1(ew, k, j, p)
                    fishs[i].calc_k2(kow, fishs[i].k_1, k, j ,p)
                    fishs[i].calc_kgb(kow, p)
                    fishs[i].calc_ke(ed, k, j, p)
                    fishs[i].calc_kd(ed, k, j ,p)


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
        if fish == None:
            break

        below.append('Phytoplankton')
        below.append('Zooplankton')
        below.append('Sediment/Detritus')

    while True:
        nextfish = find_next_fish(fishs, new_order, below)
        if type(nextfish) == obj.Fish:
            new_order.append(nextfish)
            below.append(nextfish.name)
            fishs.remove(nextfish)
        if nextfish == None:
            break

    return new_order


def find_base_fish(fishs):

    basefish = None

    for fish in fishs:
        # counter for non 0 non phyto zoop diet
        count = 0
        for entry in fish.diet_frac:

            if (entry[0] != 'Sediment/Detritus' and entry[0] != 'Phytoplankton' and entry[0] != 'Zooplankton') and entry[1] > 0:
                count += 1

        if count == 0:
            basefish = fish

    return basefish


def find_next_fish(possible_fish, current_fish, below):

    nextfish = None

    for fish in possible_fish:
        count = 0
        for entry in fish.diet_frac:
            if (entry[0] not in below) and entry[1] > 0:

                count += 1

        if count == 0:
            nextfish = fish


    return nextfish


def solve_steady_state(region, chemicals, phytos, zoops, fishs):

    # nested dictonary where can look up first by region then chemical then animal to find concentration
    conc_log = {}
    conc_log[region.name] = {}

    phytolog = conc_log[region.name]
    Cb = []
    #assuming one phyto
    phyto = phytos[0]
    phytolog[phyto.name] = {}
    for i in range (len(chemicals)):
        cwd = chemicals[i].Cwdo[0]
        conc_in_phyto = phyto.solve_steady_state(cwd, i)
        Cb.append(conc_in_phyto)
        phytolog[phyto.name][chemicals[i].name] = conc_in_phyto
    phyto.Cb = Cb

    zooplog = conc_log[region.name]
    for zoo in zoops:
        Cb = []
        zooplog[zoo.name] = {}
        for i in range (len(chemicals)):
            # since single region all chemical properties index by 0
            phi = chemicals[i].phi[0]
            Cwdp = chemicals[i].Cwp[0]
            Cwdo = chemicals[i].Cwdo[0]
            Cwto = chemicals[i].Cwto[0]
            phyto_con = conc_log[region.name]['Phytoplankton'][chemicals[i].name]

            conc_in_zoops = zoo.solve_steady_state(phi, i, Cwto, Cwdp, Cwdo, phyto_con)
            Cb.append(conc_in_zoops)
            zooplog[zoo.name][chemicals[i].name] = conc_in_zoops

        zoo.Cb = Cb

    fishlog = conc_log[region.name]

    for i in range (len(fishs)):
        fishlog[fishs[i].name] = {}
        Cb = []
        for j in range (len(chemicals)):
            phi = chemicals[j].phi[0]
            Cwp = chemicals[j].Cwp[0]
            Cwdo = chemicals[j].Cwdo[0]
            con_in_i = fishs[i].solve_steady_state(phi, j, Cwp, Cwdo, fishlog, chemicals[j])
            Cb.append(con_in_i)
            fishlog[fishs[i].name][chemicals[j].name] = con_in_i
        fishs[i].Cb = Cb

    return conc_log

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
                cont.plot_info()


# are we solving steady state on single region, solving with time on single region, or time on multiple regions
def filter_cases(filename):

    model_para, all_data, time_steps = FR_Input_Output.convert_to_lists('sheets/input/tests/testy_test.xlsx')

    v_iter = int(model_para[0])
    u_iter = int(model_para[1])

    # if steady state problem wants to be solved
    if len(all_data[0]) == 1 and model_para[8] == 'YES':

        dictionaries = []
        set_all_h_and_s(model_para, all_data)
        inner_count = 0
        u_count = 0
        print('percentage done: ')
        while (u_count < u_iter):
            u_count += 1
            v_count = 0
            while (v_count < v_iter):
                concentration_log, bundle = single_bio_iter(model_para,all_data,time_steps, u_count, v_count)
                v_count += 1
                inner_count += 1
                dictionaries.append(concentration_log)
                print('\r' + str(100 * (inner_count / (v_iter * u_iter))), end='')

        # if no statistical simulation
        if len(dictionaries) == 1:

            FR_Input_Output.deter_write_output(bundle,filename)

        # if we ran bio monte carlo
        else:
            results_dic = pr.make_result_dist(dictionaries)
            return results_dic

    #if len(all_data[0]) > 1 and model_para == 'NO':






def single_bio_iter(model_para, all_data, time_steps, u_count=0, v_count=0):

    sample_data = model_para[0:3]
    reg_data = all_data[0]
    temp_data = all_data[1]
    chem_data = all_data[2]
    set_all_h_and_s(sample_data, all_data)
    regions = init_region(time_steps, reg_data,temp_data,10, 6)
    chemicals = init_chems_unconditional(chem_data, len(regions), 10, 6)
    chemicals = init_chems_conditional(chemicals,regions)
    phytos = init_phyto(all_data[3],chemicals,10,6)
    zoops = init_zoop(all_data[4],regions,chemicals,phytos[0], time_steps, 10, 6)
    fishs, tempfishs = init_fish_pre_region(all_data[5], regions, chemicals, phytos[0], zoops, all_data[6], time_steps, 10, 6)

    ##
    ## figure out what region we are in somewhere in here (insert spatial)
    ##

    fishs = init_fish_post_region(fishs,tempfishs, regions, chemicals, time_steps)

    if model_para[8] == 'YES':
        conc_log = solve_steady_state(regions[0], chemicals, phytos, zoops, fishs)
        return conc_log, [regions, fishs, zoops, phytos, chemicals]

    #else solve non-steady state


def test_changing_time():
    model_para, all_data, time_steps = FR_Input_Output.convert_to_lists('sheets/input/tests/testy_test.xlsx')
    sample_data = model_para[0:3]
    reg_data = all_data[0]
    temp_data = all_data[1]
    chem_data = all_data[2]
    set_all_h_and_s(sample_data, all_data)
    regions = init_region(time_steps, reg_data,temp_data,10, 6)
    chemicals = init_chems_unconditional(chem_data, len(regions), 10, 6)
    chemicals = init_chems_conditional(chemicals,regions)
    phytos = init_phyto(all_data[3],chemicals,10,6)
    zoops = init_zoop(all_data[4],regions,chemicals,phytos[0], time_steps, 10, 6)
    fishs, tempfishs = init_fish_pre_region(all_data[5], regions, chemicals, phytos[0], zoops, all_data[6], time_steps, 10, 6)
    ##
    ## figure out what region we are in somewhere in here (insert spatial)
    ##
    fishs = init_fish_post_region(fishs,tempfishs, regions, chemicals, time_steps)
    print(fishs)



#test_changing_time()
filter_cases('sheets/input/tests/testy_test.xlsx')