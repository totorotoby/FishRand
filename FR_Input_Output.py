import xlrd, xlsxwriter
import networkx as nx
import prob as pr
import Time_parser
import numpy as np
from spatial import HotSpot


# Reads in all data from spread sheets


def convert_to_lists(filename):

    all_sheets = xlrd.open_workbook(filename)

    reg_len = 9  # number of inputs per region

    
    ## getting model param from first excel sheet ##
    
    model_para = []
    model_sheet = all_sheets.sheet_by_index(0)
    para_col = model_sheet.col(1)

    get_model_para(para_col, model_para)
    
    ## getting regional parameters from excel ##

    reg_sheet = all_sheets.sheet_by_index(1)

    entry_col = reg_sheet.col(1)

    region_data = []
    get_data(entry_col, reg_len, region_data)


    ## temp data from excel ##
    
    temp_data = []
    temp_sheet = all_sheets.sheet_by_index(2)

    num_timestep, time_per_step = Time_parser.num_steps(model_para[4], model_para[5], model_para[6])
    
    get_temp_data(temp_data, temp_sheet, len(region_data), num_timestep)

    ## chemical properties from excel ##
    

    chem_len = 9
    chem_sheet = all_sheets.sheet_by_index(3)
    chem_entry_col = chem_sheet.col(1)
    chem_data = []
    
    get_data(chem_entry_col, chem_len, chem_data)

    ## chemical concentrations from excel ##
    
    con_sheet = all_sheets.sheet_by_index(4)
    con_data = []
    
    get_con_data(con_sheet, con_data, len(region_data), len(chem_data), num_timestep)

    
    ## organism properties from excel ##
    
    org_sheet = all_sheets.sheet_by_index(5)

    f_len = 11  # number of inputs per fish
    zo_len = 11  # number of inputs per zoo
    ph_len = 4  # number of inputs per phyto

    fish_entry_col = org_sheet.col(1)
    invert_entry_col = org_sheet.col(3)
    zoop_entry_col = org_sheet.col(5)
    phyto_entry_col = org_sheet.col(7)

    fish_data = []
    invert_data = []
    zoop_data = []
    phyto_data = []

    get_data(fish_entry_col, f_len, fish_data)
    get_data(invert_entry_col, f_len, invert_data)
    get_data(zoop_entry_col, zo_len, zoop_data)
    get_data(phyto_entry_col, ph_len, phyto_data)

    ## organisms diets from excel ##
    
    diet_data = {}
    entrysize = len(fish_data) + len(invert_data) + 5
    diet_sheet = all_sheets.sheet_by_index(6)
    get_diet_data(diet_sheet, diet_data, entrysize)

    #Used to visualize foodweb
    foodweb_graph = foodweb_to_network_struc(diet_data)

    ## Fish migration data from excel ##
    
    mig_data = {}
    mig_sheet = all_sheets.sheet_by_index(7)
    get_mig_data(mig_data, mig_sheet)

    
    ## geometric site, and hotspot data from excel ##
    
    sites_sheet = all_sheets.sheet_by_index(8)
    sites_data = get_sites_data(sites_sheet)

    ## pass to program ##

    total = [region_data, temp_data, chem_data, phyto_data, zoop_data, invert_data, fish_data, diet_data, mig_data, con_data]

    return model_para, total, num_timestep, time_per_step, sites_data, foodweb_graph


def get_model_para(para_col, model_para):

    for i in range(1, len(para_col)):
        model_para.append(para_col[i].value)


def get_temp_data(temp_data, temp_sheet, reg_len, timesteps):

    #loop over regions by row
    for i in range (1, reg_len + 1):
        r_temps = []
        row = temp_sheet.row(i)
        #loop over timesteps by column
        for j in range (1, timesteps + 1):
            data_get_helper(row[j], r_temps)
            
        temp_data.append(r_temps)

    return temp_data
    
def get_data(entry_col, instance_len, new_list):

    new_entry = []
    for i in range(len(entry_col)):
        if i % (instance_len + 1) == 0:
            entry = entry_col[i].value
            if len(new_entry) != 0:
                new_list.append(new_entry)
            new_entry = []
            if entry == 'END':
                break
        else:
            if i % (instance_len + 1) == 1:
                if entry_col[i].value == '':
                    break
                else:
                    new_entry.append(entry_col[i].value)
            else:
                data_get_helper(entry_col[i], new_entry)
    if len(new_entry) != 0:
        new_list.append(new_entry)


def data_get_helper(preentry, new_entry):

    entry = preentry.value
    
    if type(entry) == float:
        new_entry.append(entry)
    elif type(entry) == str and entry != '':
        split = entry.split(" ", 2)
        ty = split[0].replace(',', '')
        name = split[1]
        params = split[2].replace('(', '')
        params = params.replace(')', '')
        params = params.split(', ')
        try:
            params = [float(i) for i in params]
        except ValueError:
            print(entry)
            print('In ' + str(new_entry[0]) +  ' something is wrong with the format of' + entry)
            exit(0)

        to_add = pr.Var(ty, name, params)
        new_entry.append(to_add)
    else:
        new_entry.append(entry)


def get_con_data(con_sheet, con_data, num_reg, num_chem, timesteps):

    
    for i in range (1, timesteps+1):
        t_cons = []
        column = con_sheet.col(1)
        for j in range (num_reg):
            r_cons = []
            for k in range (0, num_chem):
                c_cons = []
                for p in range (1,5):
                    
                    # This is a mess so...
                    # num_chem * 4 * j counts all red spaces in all of the  prior regional sections
                    # 1 + j counts the number of yellow regional labels
                    # k*4 counts the number of red spaces from this region to the partiular chemical
                    # 2*j + 1 + k counts the number of green labels for chemicals
                    # p counts the remaining red spaces down to the type of concentration to append
                    #print(column[(((num_chem * 4 * j) + (1 + j) ) + ( (k * 4) + ( (2 * j) + (1 + k) ) ) + p)])
                    #print((((num_chem * 4 * j) + (1 + j) ) + ( (k * 4) + ( (2 * j) + (1 + k) ) ) + p))
                    data_get_helper(column[(1+j) + ((j * num_chem) + (1+k)) + ((j * num_chem * 4) + (k * 4) + p)], c_cons)
                r_cons.append(c_cons)
            t_cons.append(r_cons)
        con_data.append(t_cons)

    
            

def get_diet_data(diet_sheet, diet_data, entrysize):

    new_name = 0
    rows = diet_sheet.nrows
    
    for i in range(rows):
        if i % entrysize == 0:
            continue
        if i % entrysize == 1:
            new_name = diet_sheet.row(i)[0].value
            diet_data[new_name] = []
        else:
            try:
                new_prey_data = [diet_sheet.row(i)[0].value, float(diet_sheet.row(i)[1].value)]
            except ValueError:
                print("you are missing an organism in a diet list or the entirety of a diet for a organism.")
                exit(0)
            diet_data[new_name].append(new_prey_data)


def get_mig_data(mig_data, mig_sheet):

    for i in range(1, mig_sheet.nrows):
        row = mig_sheet.row(i)
        mig_data[row[0].value] = [row[i].value for i in range(1, len(row))]


def foodweb_to_network_struc(foodweb):

    network_struc = {}
    nodes = [entry for entry in list(foodweb.values())[0]]

    for node in nodes:
            network_struc[node[0]] = []
    
    no_eat = ['Sediment/Detritus', 'Zooplankton', 'Phytoplankton']
    network_struc['Zooplankton'].append('Phytoplankton')
    for predator, prey in network_struc.items():
        if predator not in no_eat:
            lookup = foodweb[predator]
            for entry in lookup:
                if entry[1] != 0:
                    prey.append(entry[0])


    foodweb_graph = nx.DiGraph(network_struc)

    return foodweb_graph

def get_sites_data(sites_sheet):

    boundary = []
    bound_row = sites_sheet.row(1)
    col = 1
    while bound_row[col].value != '':
        coord = [float(i) for i in (bound_row[col].value.replace('', '').split(','))]
        boundary.append(coord)
        col += 1

    sites = []
    site_col_name = sites_sheet.col(1)
    site_col_coord = sites_sheet.col(2)
    row = 3
    while site_col_coord[row].value != '' and site_col_name[row].value != '':
        site = (site_col_name[row].value, [float(i) for i in (site_col_coord[row].value.replace('', '').split(', '))])
        sites.append(site)
        row += 1

    hotspots = []
    row_num = 40
    row = sites_sheet.row(row_num)
    deftype = sites_sheet.row(38)[3].value
    while row[0].value != '':
        col = 3
        row = sites_sheet.row(row_num)
        defin = []
        while row[col].value != '':
            defin.append(sites_sheet.row(row_num)[col].value)
            col += 1
        hotspot = HotSpot(sites_sheet.row(row_num)[0].value,deftype,
                          sites_sheet.row(row_num)[1].value, sites_sheet.row(row_num)[2].value, defin)
        if len(hotspot.definition) != 0:
            hotspots.append(hotspot)
        row_num += 1

    draws = int(sites_sheet.row(38)[1].value)
    return [boundary, sites, hotspots, draws]


def write_output_steady(total_cons, output_name, stop, dist_type):


    types = ["Normal", 'Lognormal', 'Uniform', 'Gamma', "KS Best"]
    index = types.index(dist_type)

    workbook = xlsxwriter.Workbook(output_name)
    bold = workbook.add_format({'bold': True})
    sheet = workbook.add_worksheet(stop)

    org_list = []
    for region, animals in total_cons.items():

        for animal, chems in animals.items():
            org_list.append(animal)
    org_len = len(org_list)
    chem_list = list(list(list(total_cons.values())[0].values())[0].items())
    chem_len = len(chem_list)
    #write down region which is dictionaries first key
    sheet.write(0,0, 'Region: ' + str(list(total_cons.keys())[0]) + ' Concentrations (ng/g ww)', bold)

    #write chemicals at the top
    for j in range(chem_len):
        sheet.write(0, j+1, chem_list[j][0])
    #write orgs on the column
    for i in range(org_len):
        sheet.write(i+1, 0, org_list[i])


    # if statstical
    if type(chem_list[0][1]) == pr.ResultDist:
        for i in range(org_len):
            for j in range(chem_len):
                if index != 4:
                    total_cons[region][org_list[i]][chem_list[j][0]].index = index
                    to_print = total_cons[region][org_list[i]][chem_list[j][0]].bestparam()
                    sheet.write(i + 1, j + 1, to_print[0])
                else:
                    sheet.write(i + 1, j + 1, total_cons[region][org_list[i]][chem_list[j][0]].best_para[0])

        sheet.write(org_len + 1, 0, 'Unfitted Mean and Std of Simulations', bold)
        for k in range(org_len):
            for j in range(chem_len):
                sheet.write(org_len+(k+2),0,org_list[k])
                sheet.write(org_len+(k+2), j+1, str(total_cons[region][org_list[k]][chem_list[j][0]].v_mean_std))

    # if nonstatstical
    else:
        for i in range(org_len):
            for j in range(chem_len):
                sheet.write(i+1,j+1, round(total_cons[region][org_list[i]][chem_list[j][0]],6))

    workbook.close()



def write_temporal_excel(array, output_name, stops, stat_flag, regional_areas, dist_type):

    types = ["Normal", 'Lognormal', 'Uniform', 'Gamma', "KS Best"]
    index = types.index(dist_type)


    workbook = xlsxwriter.Workbook(output_name)
    bold = workbook.add_format({'bold': True})
    big = workbook.add_format({'font_size': 14, 'bold': True})

    if stat_flag == 0:

        lower_org_list = []
        for animal, chems in array[0][1].items():
                lower_org_list.append(animal)

        upper_org_list = []
        for animal, chems in array[0][2].items():
            upper_org_list.append(animal)
        chem_list = list(list(array[0][1].values())[0].keys())

        reg_list =list(array[0][0][0].keys())
        #write down region which is dic

        for i in range(len(stops)):
            sheet = workbook.add_worksheet('Timestep ' + str(stops[i]))
            data_at_time = array[i]
            lower_non_avg = data_at_time[0][0]
            lower_avg = data_at_time[1]
            upper_avg = data_at_time[2]
            sheet.write(0,0, 'Lower Food Web Concentrations by Region (ng/g ww)', big)
            # write chemicals at the top
            for j in range(len(chem_list)):
                sheet.write(0, j + 1, chem_list[j])
            # write orgs on the column
            count = 0
            for pair in lower_non_avg.items():
                sheet.write(count + 1, 0, pair[0], bold)
                count += 1
                for animal in pair[1].keys():
                    sheet.write(count + 1, 0, animal)
                    count +=1

            for i in range(len(reg_list)):
                for j in range(len(lower_org_list)):
                    for k in range(len(chem_list)):
                        # row first term everything before, then over organisms
                        sheet.write((1 + i * (1 + len(lower_org_list))) + (1+j) , k+1, round(lower_non_avg[reg_list[i]][lower_org_list[j]][chem_list[k]],6))



            sheet.write(0, len(chem_list) + 1, 'Lower Food Web Concentrations Average Concentrations Weighted by Regional Area (ng/g ww)', big)
            for j in range(len(chem_list)):
                sheet.write(0, len(chem_list) + 1 + (j + 1), chem_list[j])
            for i in range(len(lower_org_list)):
                for j in range(len(chem_list)):
                    sheet.write(2 + i, 2 + len(chem_list) + j, round(lower_avg[lower_org_list[i]][chem_list[j]],6))

            sheet.write(count + 1, 0, 'Upper Food Web Concentrations Averaged over Populations (ng/g ww)', big)
            for i in range(len(upper_org_list)):
                sheet.write(count + 2 + i, 0, upper_org_list[i])
                for j in range (len(chem_list)):
                    sheet.write(count + 2 + i, j+1, round(upper_avg[upper_org_list[i]][chem_list[j]],6))


    if stat_flag == 1:

        normed_reg_areas = np.asarray(regional_areas) / sum(regional_areas)

        for i in range(len(stops)):

            sheet = workbook.add_worksheet('Timestep ' + str(stops[i]))

            lower_dists = array[i][0]
            upper_dists = array[i][1]

            reg_list = list(lower_dists.keys())
            lower_org_list = list(list(lower_dists.values())[0].keys())
            upper_org_list = list(upper_dists.keys())
            chem_list = list(list(upper_dists.values())[0].keys())

            sheet.write(0, 0, 'Lower Food Web Concentrations by Region (ng/g ww)', big)

            for j in range(len(chem_list)):

                sheet.write(0, j + 1, chem_list[j])

            count = 1
            for reg in reg_list:
                sheet.write(count, 0, reg, bold)
                count += 1
                for animal in lower_org_list:
                    sheet.write(count, 0, animal)
                    count += 1

            for i in range(len(reg_list)):
                for j in range(len(lower_org_list)):
                    for k in range(len(chem_list)):
                        if index != 4:
                            lower_dists[reg_list[i]][lower_org_list[j]][chem_list[k]].index = index
                            toprint = lower_dists[reg_list[i]][lower_org_list[j]][chem_list[k]].bestparam()
                            sheet.write((1 + i * (1 + len(lower_org_list))) + (1+j), k + 1, toprint[0])
                        else:
                            sheet.write((1 + i * (1 + len(lower_org_list))) + (1+j), k + 1,
                                        lower_dists[reg_list[i]][lower_org_list[j]][chem_list[k]].best_para[0])

            sheet.write(0, len(chem_list) + 1, 'Lower Food Web Mean Concentrations Weighted by Regional Area (ng/g ww)', big)
            sheet.write(1, len(chem_list) + 1, 'All Regions', bold)

            for j in range(len(chem_list)):
                sheet.write(0, len(chem_list) + 1 + (j + 1), chem_list[j])


            for i in range(len(chem_list)):
                for j in range(len(lower_org_list)):
                    weight_avg = 0
                    for k in range(len(reg_list)):
                        weight_avg =+normed_reg_areas[k]*lower_dists[reg_list[k]][lower_org_list[j]][chem_list[i]].v_mean_std[0]

                    sheet.write(2 + j, 2 + len(chem_list) + i, weight_avg)

            sheet.write(count + 1, 0, 'Upper Food Web Concentrations (ng/g ww)', big)
            count += 1

            for i in range(len(upper_org_list)):
                sheet.write(count + 1 + i, 0,upper_org_list[i])
                for j in range(len(chem_list)):
                    if index != 4:
                        upper_dists[upper_org_list[i]][chem_list[j]].index = index
                        toprint = upper_dists[upper_org_list[i]][chem_list[j]].bestparam()
                        sheet.write(count + 1 + i, j+1, toprint[0])
                    else:
                        sheet.write(count + 1 + i, j + 1, upper_dists[upper_org_list[i]][chem_list[j]].best_para[0])

            sheet.write(count, 1 + len(chem_list), 'Upper Food Web Concentrations (ng/g ww) (Mean and standard Deviation of Samples)', big)

            for i in range(len(upper_org_list)):
                for j in range(len(chem_list)):
                    sheet.write(count + 1 + i, j+ len(chem_list) + 2, str(upper_dists[upper_org_list[i]][chem_list[j]].v_mean_std))


    workbook.close()
