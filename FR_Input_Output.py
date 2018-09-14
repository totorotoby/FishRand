import xlrd, xlsxwriter
import datetime
#import Classes as cs
import prob as pr
import Time_parser
from spatial import HotSpot



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
    zoop_entry_col = org_sheet.col(4)
    zoop_dist_col = org_sheet.col(5)
    phyto_entry_col = org_sheet.col(8)
    phyto_dist_col = org_sheet.col(9)

    fish_data = []
    zoop_data = []
    phyto_data = []

    get_data(fish_entry_col, fish_dist_col, f_len, fish_data)
    get_data(zoop_entry_col, zoop_dist_col, zo_len, zoop_data)
    get_data(phyto_entry_col, phyto_dist_col, ph_len, phyto_data)

    diet_data = {}
    entrysize = len(fish_data) + 5
    diet_sheet = all_sheets.sheet_by_index(5)
    get_diet_data(diet_sheet, diet_data, entrysize)

    mig_data = {}
    mig_sheet = all_sheets.sheet_by_index(6)
    get_mig_data(mig_data, mig_sheet)

    sites_sheet = all_sheets.sheet_by_index(7)
    sites_data = get_sites_data(sites_sheet)

    total = [region_data, temp_data, chem_data, phyto_data, zoop_data, fish_data, diet_data, mig_data]

    return model_para, total, num_timestep, time_per_step, sites_data


def get_model_para(para_col, model_para):

    for i in range(1, len(para_col)):
        model_para.append(para_col[i].value)


def get_temp_data(temp_data, temp_sheet, reg_len):

    for i in range(2, reg_len + 2):
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
            new_prey_data = [diet_sheet.row(i)[1].value, diet_sheet.row(i)[2].value]
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
    #print(sites_sheet.row(row_num + 1))
    while row[0].value != '':
        col = 4
        row = sites_sheet.row(row_num)
        defin = []
        while row[col].value != '':
            defin.append(sites_sheet.row(row_num)[col].value)
            col += 1
        hotspot = HotSpot(sites_sheet.row(row_num)[0].value,sites_sheet.row(row_num)[1].value,
                          sites_sheet.row(row_num)[2].value, sites_sheet.row(row_num)[3].value, defin)
        if len(hotspot.defintion) != 0:
            hotspots.append(hotspot)
        row_num += 1

    return [boundary, sites, hotspots]


def deter_write_output(*args):
    output_name = "sheets/output/FR_Model_" + '{:%Y-%m-%d %H:%M}'.format(datetime.datetime.now()) +\
                   '_from_' + str(args[5]) + '.xls'
    workbook = xlsxwriter.Workbook(output_name)
    worksheet = workbook.add_worksheet()
    num_org = len(args[1]) + len(args[3]) + len(args[2])

    for i in range(len(args[4])):
        worksheet.write(i + 1, 0, args[4][i].name)
    for i in range(len(args[0])):
        write_region = (i * (num_org + 1))
        worksheet.write(0, write_region, args[0][i].name)
        for j in range(len(args[3])):
            write_phyto = (write_region + 1) + j
            worksheet.write(0, write_phyto, args[3][j].name)
            for p in range(len(args[3][j].Cb)):
                chem_write = p + 1
                worksheet.write(chem_write, write_phyto, args[3][j].Cb[p])
        for j in range(len(args[2])):
            write_zoop = ((write_region + write_phyto + 1) + j)
            worksheet.write(0, write_zoop, args[2][j].name)
            for p in range(len(args[2][j].Cb)):
                chem_write = p + 1
                worksheet.write(chem_write, write_zoop, args[2][j].Cb[p])
        for j in range(len(args[1])):
            write_fish = ((write_region + write_phyto + write_zoop) + j)
            worksheet.write(0, write_fish, args[1][j].name)
            for p in range(len(args[1][j].Cb)):
                chem_write = p + 1
                worksheet.write(chem_write, write_fish, args[1][j].Cb[p])

    workbook.close()

