import xlrd, xlsxwriter
import datetime
#import Classes as cs
import prob as pr
import Time_parser
import numpy as np
from spatial import HotSpot

#TODO Fix output parameters

# Reads in all data from spread sheets


def convert_to_lists(filename):

    all_sheets = xlrd.open_workbook(filename)

    f_len = 11  # number of inputs per fish
    zo_len = 11  # number of inputs per zoo
    ph_len = 4  # number of inputs per phyto
    reg_len = 9  # number of inputs per region

    model_para = []
    model_sheet = all_sheets.sheet_by_index(0)
    para_col = model_sheet.col(1)

    get_model_para(para_col, model_para)
    #print(model_para)
    reg_sheet = all_sheets.sheet_by_index(1)

    entry_col = reg_sheet.col(1)
    dist_col = reg_sheet.col(2)

    region_data = []

    get_data(entry_col, dist_col, reg_len, region_data)

    temp_data = []
    temp_sheet = all_sheets.sheet_by_index(2)

    num_timestep, time_per_step = Time_parser.num_steps(model_para[4], model_para[5], model_para[6])
    num_timestep = num_timestep + 1
    get_temp_data(temp_data, temp_sheet, len(region_data))

    chem_len = 5 + (len(region_data)*3)

    chem_sheet = all_sheets.sheet_by_index(3)

    entry_col = chem_sheet.col(1)
    dist_col = chem_sheet.col(2)

    chem_data = []

    get_chem_data(entry_col, dist_col, chem_len, chem_data, len(region_data))

    org_sheet = all_sheets.sheet_by_index(4)

    fish_entry_col = org_sheet.col(1)
    fish_dist_col = org_sheet.col(2)
    invert_entry_col = org_sheet.col(4)
    invert_dist_col = org_sheet.col(5)
    zoop_entry_col = org_sheet.col(7)
    zoop_dist_col = org_sheet.col(8)
    phyto_entry_col = org_sheet.col(10)
    phyto_dist_col = org_sheet.col(11)

    fish_data = []
    invert_data = []
    zoop_data = []
    phyto_data = []

    get_data(fish_entry_col, fish_dist_col, f_len, fish_data)
    get_data(invert_entry_col, invert_dist_col, f_len, invert_data)
    get_data(zoop_entry_col, zoop_dist_col, zo_len, zoop_data)
    get_data(phyto_entry_col, phyto_dist_col, ph_len, phyto_data)
    diet_data = {}
    entrysize = len(fish_data) + len(invert_data) + 5
    diet_sheet = all_sheets.sheet_by_index(5)
    get_diet_data(diet_sheet, diet_data, entrysize)

    mig_data = {}
    mig_sheet = all_sheets.sheet_by_index(6)
    get_mig_data(mig_data, mig_sheet)

    sites_sheet = all_sheets.sheet_by_index(7)
    sites_data = get_sites_data(sites_sheet)

    total = [region_data, temp_data, chem_data, phyto_data, zoop_data, invert_data, fish_data, diet_data, mig_data]

    return model_para, total, num_timestep, time_per_step, sites_data


def get_model_para(para_col, model_para):

    for i in range(1, len(para_col)):
        model_para.append(para_col[i].value)


def get_temp_data(temp_data, temp_sheet, reg_len):
    for i in range(2, (reg_len * 2) + 1, 2):
        try:
            row = temp_sheet.row(i)

        except:
            break
        regional_temps = []
        for j in range(1, len(row), 2):
            entry = row[j].value
            if type(entry) == float:
                regional_temps.append(entry)
            if type(entry) == str:
                entry = entry.split(', ')
                param = row[j+1].value.split(', ')
                param = [float(i) for i in param]
                toadd = pr.Var(entry[0], entry[1], param)
                regional_temps.append(toadd)

        temp_data.append(regional_temps)


def get_data(entry_col, dist_col, instance_len, new_list):

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
                new_entry.append(entry_col[i].value)

            else:
                data_get_helper(entry_col[i], dist_col[i], new_entry)
    if len(new_entry) != 0:
        new_list.append(new_entry)


def data_get_helper(preentry, dist, new_entry):

    entry = preentry.value
    dist_par = dist.value
    if entry == '' and type(dist_par) == str and dist_par != '':
        print('You must specify a distribution where you provide ' + dist_par)
        exit(0)
    elif type(entry) == float and dist_par == '':
        new_entry.append(entry)
    elif type(dist_par) == str and entry != '':
        entry = entry.split(',')
        ty = entry[0].strip()
        dist_name = entry[1].strip()

        dist_par = dist_par.split(',')
        dist_par = map(str.strip, dist_par)
        try:
            dist_par = [float(i) for i in dist_par]
        except ValueError:
            print('In ' + str(new_entry[0]) + ', ' + str(list(dist_par)) + ' something is wrong with format.')
            exit(0)

        to_add = pr.Var(ty, dist_name, dist_par)
        new_entry.append(to_add)
    else:
        new_entry.append(dist_par)


def get_diet_data(diet_sheet, diet_data, entrysize):
    # getting Diet data #
    new_name = 0
    rows = diet_sheet.nrows

    for i in range(rows):
        if i % entrysize == 0:
            continue
        if i % entrysize == 1:
            new_name = diet_sheet.row(i)[1].value
            diet_data[new_name] = []
        else:
            new_prey_data = [diet_sheet.row(i)[1].value, float(diet_sheet.row(i)[2].value)]
            diet_data[new_name].append(new_prey_data)


def get_chem_data(entry_col, dist_col, chem_len, chem_data, num_regions):

    for i in range(len(entry_col)):
        if i % chem_len == 0:
            new_chem = []
            sed_con = []
            total_con = []
            dis_con = []
        if i % chem_len == 1:
            new_chem.append(entry_col[i].value)
            data_get_helper(entry_col[i+1], dist_col[i+1], new_chem)
        if i % chem_len == 3:
            for j in range(i, i+num_regions*3):
                if (j-i) % 3 == 0:
                    data_get_helper(entry_col[j], dist_col[j], sed_con)
                if (j - i) % 3 == 1:
                    data_get_helper(entry_col[j], dist_col[j], total_con)
                if (j - i) % 3 == 2:
                    data_get_helper(entry_col[j], dist_col[j], dis_con)
            new_chem.append(sed_con)
            new_chem.append(total_con)
            new_chem.append(dis_con)
        if i % chem_len == chem_len-2:
            data_get_helper(entry_col[i], dist_col[i], new_chem)
            data_get_helper(entry_col[i+1], dist_col[i+1], new_chem)
            chem_data.append(new_chem)


def get_mig_data(mig_data, mig_sheet):

    for i in range(1, mig_sheet.nrows):
        row = mig_sheet.row(i)
        mig_data[row[0].value] = [row[i].value for i in range(1, len(row))]


def get_sites_data(sites_sheet):

    boundary = []
    bound_row = sites_sheet.row(1)
    col = 1
    while bound_row[col].value != '':
        coord = [float(i) for i in (bound_row[col].value.replace(' ', '').split(','))]
        boundary.append(coord)
        col += 1

    sites = []
    site_col_name = sites_sheet.col(1)
    site_col_coord = sites_sheet.col(2)
    row = 3
    while site_col_coord[row].value != '' and site_col_name[row].value != '':
        site = (site_col_name[row].value, [float(i) for i in (site_col_coord[row].value.replace(' ', '').split(','))])
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


def write_output_steady(total_cons, output_name, stop):
    #print(total_cons, output_name, stop)
    workbook = xlsxwriter.Workbook(output_name)
    bold = workbook.add_format({'bold': True})
    sheet = workbook.add_worksheet(stop)

    # becasue non-statistical
    total_cons = total_cons[0]
    org_list = []
    for region, animals in total_cons.items():

        for animal, chems in animals.items():
            org_list.append(animal)
    org_len = len(org_list)
    chem_list = list(list(list(total_cons.values())[0].values())[0].items())
    chem_len = len(chem_list)
    #write down region which is dictionaries first key
    sheet.write(0,0, 'Region: ' + str(list(total_cons.keys())[0]), bold)

    #write chemicals at the top
    for j in range(chem_len):
        sheet.write(0, j+1, chem_list[j][0])
    #write orgs on the column
    for i in range(org_len):
        sheet.write(i+1, 0, org_list[i])

    # if nonstatstical
    if type(chem_list[0][1]) == np.float64:
        for i in range(org_len):
            for j in range(chem_len):
                sheet.write(i+1,j+1, total_cons[region][org_list[i]][chem_list[j][0]])

    # if statstical
    if type(chem_list[0][1]) == pr.ResultDist:
        for i in range(org_len):
            for j in range(chem_len):
                sheet.write(i + 1, j + 1, total_cons[region][org_list[i]][chem_list[j][0]].best_para[0])

        sheet.write(org_len + 1, 0, 'Unfitted Mean and Std of Simulations', bold)
        for k in range(org_len):
            for j in range(chem_len):
                sheet.write(org_len+(k+2),0,org_list[k])
                sheet.write(org_len+(k+2), j+1, str(total_cons[region][org_list[k]][chem_list[j][0]].v_mean_std))



def write_temporal_excel(array, output_name, stops, stat_flag, regional_areas):

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
        #print(array[0][0])
        reg_list =list(array[0][0][0].keys())
        #write down region which is dic

        for i in range(len(stops)):
            sheet = workbook.add_worksheet('Timestep ' + str(stops[i]))
            data_at_time = array[i]
            lower_non_avg = data_at_time[0][0]
            lower_avg = data_at_time[1]
            upper_avg = data_at_time[2]
            sheet.write(0,0, 'Lower Food Web Concentrations by Region', big)
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

            count1 = -1
            for i in range(len(reg_list)):
                count1 += 1
                for j in range(len(lower_org_list)):
                    for k in range(len(chem_list)):
                        #if
                        sheet.write(count1 + 2 + (i + j) + (len(lower_non_avg) * i), k+1, lower_non_avg[reg_list[i]][lower_org_list[j]][chem_list[k]])
                        #print()


            sheet.write(0, len(chem_list) + 1, 'Lower Food Web Concentrations Average Concentrations Weighted by Regional Area', big)
            for j in range(len(chem_list)):
                sheet.write(0, len(chem_list) + 1 + (j + 1), chem_list[j])


            for i in range(len(lower_org_list)):
                for j in range(len(chem_list)):
                    sheet.write(1 + i, 2 + len(chem_list) + j, lower_avg[lower_org_list[i]][chem_list[j]])

            sheet.write(count + 1, 0, 'Upper Food Web Concentrations Averaged over Populations', big)
            for i in range(len(upper_org_list)):
                sheet.write(count + 2 + i, 0, upper_org_list[i])
                for j in range (len(chem_list)):
                    sheet.write(count + 2 + i, j+1, upper_avg[upper_org_list[i]][chem_list[j]])


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

            sheet.write(0, 0, 'Lower Food Web Concentrations by Region', big)

            for j in range(len(chem_list)):
                sheet.write(0, j + 1, chem_list[j])

            count = 1
            for reg in reg_list:
                sheet.write(count, 0, reg, bold)
                count += 1
                for animal in lower_org_list:
                    sheet.write(count, 0, animal)
                    count += 1

            count1 = -1
            for i in range(len(reg_list)):
                count1 += 1
                for j in range(len(lower_org_list)):
                    for k in range(len(chem_list)):
                        # if
                        sheet.write(count1 + 2 + (i + j) + (len(lower_dists) * i), k + 1,
                                    lower_dists[reg_list[i]][lower_org_list[j]][chem_list[k]].best_para[0])

            sheet.write(0, len(chem_list) + 1, 'Lower Food Web Mean Concentrations Weighted by Regional Area', big)
            sheet.write(1, len(chem_list) + 1, 'All Regions', bold)

            for j in range(len(chem_list)):
                sheet.write(0, len(chem_list) + 1 + (j + 1), chem_list[j])


            for i in range(len(chem_list)):
                for j in range(len(lower_org_list)):
                    weight_avg = 0
                    for k in range(len(reg_list)):
                        weight_avg =+normed_reg_areas[k]*lower_dists[reg_list[k]][lower_org_list[j]][chem_list[i]].v_mean_std[0]

                    sheet.write(2 + j, 2 + len(chem_list) + i, weight_avg)

            sheet.write(count + 1, 0, 'Upper Food Web Concentrations', big)
            count += 1

            for i in range(len(upper_org_list)):
                sheet.write(count + 1 + i, 0,upper_org_list[i])
                for j in range(len(chem_list)):
                    sheet.write(count + 1 + i, j+1, upper_dists[upper_org_list[i]][chem_list[j]].best_para[0])

            sheet.write(count, 1 + len(chem_list), 'Upper Food Web Concentrations (Mean and standard Deviation of Samples)', big)

            for i in range(len(upper_org_list)):
                for j in range(len(chem_list)):
                    sheet.write(count + 1 + i, j+ len(chem_list) + 2, str(upper_dists[upper_org_list[i]][chem_list[j]].v_mean_std))

            workbook.close()























































