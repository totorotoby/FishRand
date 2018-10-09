import Bioaccum
import FR_Input_Output
import prob as pr
import spatial
import numpy as np
import copy
from matplotlib import pyplot as plt


def loc_setup(fishs, boundary, regions, hotspots, draw_num):

    fl = []
    for fish in fishs:
        fl.append(fish[0])

    setups = []
    for f in fl:
        setups.append(spatial.location_step(boundary,regions, hotspots, f, draw_num))

    return setups, fl


def get_locs_matrix(loc_setups, draws, mig_data, timestep):

    locations = []
    for info in loc_setups:
        loc = spatial.new_draw(info[0], info[1], info[2], info[3], draws)
        locations.append(loc)
    locations = np.asarray(locations)


    #going back and adding in probabilities of migrating out of the pond
    indices = [k for k in range(draws)]
    for i in range(len(locations)):

        fish_name = loc_setups[i][0][0].fish
        abund = mig_data[fish_name][timestep]
        percentage_fish_out = 1 - abund
        number_out = int(draws * percentage_fish_out)
        out_indices = np.random.choice(indices, size=number_out, replace=False)
        for j in range(len(out_indices)):
            locations[i][out_indices[j]] = 1000


    return locations


def get_fish_dic(fish_cons, lower_cons, chem_name, fish_name, region_areas):

    fish_dic = {}
    lower_avg_dic = copy.deepcopy(list(lower_cons[0].items())[0][1])
    if len(fish_cons) == 1:
        fish_cons = fish_cons[0]
        for i in range(len(fish_cons)):
            single_fish_avg_cons = []
            fish_dic[fish_name[i]] = {}
            for k in range(len(chem_name)):
                chem_con_over_pops = [pop_con[k] for pop_con in fish_cons[i]]
                avg_chem_over_pops = np.mean(chem_con_over_pops)
                fish_dic[fish_name[i]][chem_name[k]] = avg_chem_over_pops

        for regions in lower_cons[0].values():
            for animal, chemicals in regions.items():
                for chem_name, con in  chemicals.items():
                    if type(lower_avg_dic[animal][chem_name]) == np.float64 or type(lower_avg_dic[animal][chem_name]) == float:
                        lower_avg_dic[animal][chem_name] = []
                        lower_avg_dic[animal][chem_name].append(con)
                    else:
                        lower_avg_dic[animal][chem_name].append(con)

        for animal, chemicals in lower_avg_dic.items():
            for chem_name, chem in chemicals.items():

                lower_avg_dic[animal][chem_name] = np.average(chem, weights = region_areas)

        return lower_avg_dic, fish_dic

    else:

        # making lower distrubtion results
        lower_result = pr.make_result_dist(lower_cons)


        # make upper results dictionary
        upper_result = {}
        chem_values = [[[] for _ in range(len(chem_name))] for _ in range(len(fish_name))]
        for i in range(len(fish_cons)):
            for j in range(len(fish_cons[i])):
                if fish_name[j] not in upper_result:
                    upper_result[fish_name[j]] = {}
                for k in range(len(fish_cons[i][j])):
                    for p in range(len(fish_cons[i][j][k])):
                        chem_values[j][p].append(fish_cons[i][j][k][p])

        for i in range(len(chem_values)):
            for j in range(len(chem_values[i])):
                turn_result = np.unique(chem_values[i][j])
                upper_result[fish_name[i]][chem_name[j]] = pr.ResultDist(turn_result, chem_name[j], fish_name[i])
        return lower_result, upper_result




def graph_by_time(data, fish_name, chem_name, time_interval):

    time_steps = len(data)
    values = []
    times = [i for i in range(time_steps)]
    for i in range(len(data)):
        con_at_t = data[i][fish_name][chem_name]
        values.append(con_at_t)

    plt.plot(times, values, 'ro', )
    plt.xlabel('Timesteps (' + time_interval + ')')
    plt.ylabel('Concentration of ' + chem_name + ' in ' + fish_name + ' (ng/g)')
    plt.show()


def get_fish_in_region(locations, fish_names, reg_len):

    f_by_region = [[] for _ in range(reg_len)]

    for fish_index in range (len(locations)):
        for pop_index in range (len(locations[fish_index])):
            region_index = locations[fish_index][pop_index]
            if region_index != 1000:
                if fish_names[fish_index] not in f_by_region[region_index]:
                    f_by_region[region_index].append(fish_names[fish_index])

    return f_by_region


def niceprint(prior_total_cons, locations):



    count = 0
    for uv_dic in prior_total_cons[0]:
        print('uv count:', count)
        count += 1
        for region, animals in uv_dic.items():
            print('\t', region)
            for animal, chemicals in animals.items():
                print('\t\t', animal)
                for chem, conc in chemicals.items():
                    print('\t\t\t', chem, ' concentation: ', conc)

    print('\n\n')
    print('locations: \n')
    for fish in locations:
        print(fish)
    print()
    count = 0
    for uv_sample in prior_total_cons[1]:

        print('uv count:', count)
        count += 1
        count1 = 0
        for fishs in uv_sample:
            print('\tfish number: ', count1)
            count1 += 1
            count2 = 0
            for population in fishs:
                print('\t\tfish population:', count2)
                count2 += 1
                count3 = 0
                for chemical in population:
                    print('\t\t\tconcentation of chemical ', count3, ' : ', chemical)
                    count3 += 1


# are we solving steady state on single region, solving with time on single region, or time on multiple regions
def filter_cases(filename, stops):

    model_para, all_data, time_steps, time_per_step, site_data, foodweb_graph = FR_Input_Output.convert_to_lists(filename)
    # set up for
    u_iter = int(model_para[0])
    v_iter = int(model_para[1])
    writing_info = []
    sample_data = model_para[0:3]
    Bioaccum.set_all_h_and_s(sample_data, all_data)
    stops_count = 0
    t = 0
    # if steady state problem wants to be solved
    if model_para[8] == 'YES':

        # since we are doing steady state
        old_fish_by_region = None
        fish_by_region = None
        prior_locations = None
        locations = None

        stat_check = Bioaccum.set_all_h_and_s(model_para, all_data)
        total_cons = Bioaccum.bio_monte_carlo_loop(model_para, all_data, 0, time_per_step, old_fish_by_region ,fish_by_region, prior_locations, locations, stat_check, u_iter, v_iter)

        if stat_check == True:
            total_cons = pr.make_result_dist(total_cons)


    elif model_para[8] == 'NO':

        #assigning mig_data for later adding to locations
        mig_data = all_data[8]
        # The setup before time iterations start includes:
            # • Setting put the regions as Polygon, objects
            # • finding the probabilities assosiated with regions and hotspot polygons
            # • getting number of draws per fish per timestep from the input
            # • Making the prior concentrations dictionary for all the animals

        boundary, regions, hotspots= spatial.setup(site_data)
        region_areas = [region[1].area for region in regions]
        loc_setups, f_names = loc_setup(all_data[6], boundary, regions, hotspots, site_data[3])
        draws = site_data[3]
        graph_data = []
        # turns all_data that is distributions into array of samples which can be iterated through with u_count or v_count
        # Will return true if there is at least one distrubtion input
        stat_check = Bioaccum.set_all_h_and_s(model_para, all_data)

        #Setup for how concentrations will be returned
        # - None Fish (prestat processed:
        # { Region { animal { chemical { concentration: x }}}
        # - Fish (prestat Processed:
        # [num montecarlo iterations: region {fish:{ {chemical: x, other chemical: y}}} ]]]

        prior_concentrations_non_fish, prior_concentrations_fish = Bioaccum.init_prior_con_dic(stat_check, u_iter*v_iter,
                                                                                               all_data[0], all_data[2],
                                                                                               all_data[3], all_data[4],
                                                                                               all_data[5], all_data[6],
                                                                                               draws)
        total_cons = [prior_concentrations_non_fish, prior_concentrations_fish]
        #################################################################################

        # The time loop
        locations = get_locs_matrix(loc_setups, draws, mig_data, t)
        fish_by_region = get_fish_in_region(locations, f_names, len(all_data[0]))
        while t < time_steps:

            # get all locations for each fish... locations matrix is fish number by list of region numbers draws long
            prior_locations = locations
            locations = get_locs_matrix(loc_setups, draws, mig_data, t)
            print(locations)
            # fish_by loc is for each region there is list of names of fish that are in that region at this time step
            old_fish_by_region = fish_by_region
            fish_by_region = get_fish_in_region(locations, f_names, len(all_data[0]))

            # do entire simulation for a time step, and return concentration dictionaries
            uv_single = Bioaccum.bio_monte_carlo_loop(model_para, all_data, t, time_per_step, old_fish_by_region ,fish_by_region, prior_locations, locations,stat_check, u_iter, v_iter,
                                                         p_dic=total_cons)

            total_cons = uv_single

            lower_cons, fish_dic = get_fish_dic(total_cons[1], total_cons[0], [chem[0] for chem in all_data[2]], [fish[0] for fish in all_data[6]], region_areas)
            graph_data.append(fish_dic)
            if t in stops:
                if stat_check == 0:
                    writing_info.append([total_cons[0], lower_cons, fish_dic])
                else:
                    writing_info.append([lower_cons, fish_dic])

            t += 1
    else:

        print('Need to set steady state in the Sample and Time Input tab to YES or NO.')


    if model_para[8] == 'YES':

        return [model_para[8], total_cons, stat_check, foodweb_graph]

    if model_para[8] == 'NO':

        return [model_para[8], writing_info, stat_check, region_areas, graph_data, model_para[6], [boundary,regions,hotspots], foodweb_graph]


    #################### STEADY STATE EXCEL WRITING ########################


def steady_state_output(total_cons, stat_check, output_name, dist_type):

        FR_Input_Output.write_output_steady(total_cons[0], output_name, 'Steady State Simulation', dist_type)

def temporal_output(stat_check, to_write, output_name, stops, region_areas, dist_type):

        #Non-statistical
        if stat_check == False:
            FR_Input_Output.write_temporal_excel(to_write, output_name, stops, 0, region_areas, dist_type)
        #statistical
        else:
            FR_Input_Output.write_temporal_excel(to_write, output_name, stops, 1, region_areas, dist_type)

