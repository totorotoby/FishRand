import xlrd, xlsxwriter
import datetime
#import Classes as cs
import prob as pr

# Reads in all data from spread sheets

def convert_to_lists(filename):

    all_sheets = xlrd.open_workbook(filename)

    f_len = 11  # number of inputs per fish
    zo_len = 11  # number of inputs per zoo
    ph_len = 4  # number of inputs per phyto
    reg_len = 10  # number of inputs per region
    chem_len = 7

    model_para = []
    model_sheet = all_sheets.sheet_by_index(0)
    para_col = model_sheet.col(1)

    get_model_para(para_col, model_para)

    reg_sheet = all_sheets.sheet_by_index(1)

    entry_col = reg_sheet.col(1)
    dist_col = reg_sheet.col(2)

    region_data = []

    get_data(entry_col,dist_col,reg_len,region_data)


    chem_sheet = all_sheets.sheet_by_index(2)

    entry_col = chem_sheet.col(1)
    dist_col = chem_sheet.col(2)

    chem_data = []

    get_data(entry_col,dist_col,chem_len,chem_data)


    org_sheet = all_sheets.sheet_by_index(3)

    fish_entry_col = org_sheet.col(1)
    fish_dist_col = org_sheet.col(2)
    zoop_entry_col = org_sheet.col(4)
    zoop_dist_col = org_sheet.col(5)
    phyto_entry_col = org_sheet.col(8)
    phyto_dist_col = org_sheet.col(9)

    fish_data = []
    zoop_data = []
    phyto_data = []

    get_data(fish_entry_col,fish_dist_col,f_len, fish_data)
    get_data(zoop_entry_col,zoop_dist_col,zo_len,zoop_data)
    get_data(phyto_entry_col, phyto_dist_col, ph_len, phyto_data)

    diet_data = {}
    entrysize = len(fish_data) + 5
    diet_sheet = all_sheets.sheet_by_index(4)
    get_diet_data(fish_data,diet_sheet,diet_data,entrysize)

    total = [region_data, chem_data, fish_data, zoop_data, phyto_data, diet_data]

    return model_para, total


def get_model_para(para_col, model_para):

    for i in range (3):
        model_para.append(para_col[i].value)

def get_data(entry_col, dist_col, instance_len, new_list):


    new_entry = []
    for i in range (len(entry_col)):
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
                entry = entry_col[i].value
                dist_par = dist_col[i].value
                if type(entry) == float and dist_par == '':
                     new_entry.append(entry)
                elif type(dist_par) == str and entry != '':
                    entry = entry.split(', ')
                    ty = entry[0]
                    dist_name = entry[1]
                    dist_par = dist_par.split(', ')
                    try:
                        dist_par = [float(i) for i in dist_par]
                    except ValueError:
                        print('In ' + new_entry[0] + ', ' + str(dist_par) + ' something is wrong with format.')
                        exit(0)
                    to_add = pr.Var(ty, dist_name, dist_par)
                    new_entry.append(to_add)
                else:
                    new_entry.append(dist_par)
    if len(new_entry) != 0:
        new_list.append(new_entry)

def get_diet_data(fish_data, diet_sheet, diet_data, entrysize):
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



def deter_write_output(regions, fish, chemicals, phyto, zoop, inputfilename):
    output_name  = "sheets/output/FR_Model_" + '{:%Y-%m-%d %H:%M}'.format(datetime.datetime.now()) + '_from_' + str(inputfilename)
    workbook = xlsxwriter.Workbook(output_name)
    worksheet = workbook.add_worksheet()
    num_org = len(fish) + len(phyto) + len(zoop)

    for i in range(len(chemicals)):
        worksheet.write(i+1, 0, chemicals[i].name)
    for i in range (len(regions)):
        write_region = (i * (num_org + 1))
        worksheet.write(0, write_region , regions[i].name)
        for j in range (len(phyto)):
            write_phyto = (write_region + 1) + j
            worksheet.write(0, write_phyto, phyto[j].name)
            for p in range (len(phyto[j].Cb)):
                chem_write = p + 1
                worksheet.write(chem_write, write_phyto, phyto[j].Cb[p])
        for j in range (len(zoop)):
            write_zoop = ((write_region + write_phyto + 1) + j)
            worksheet.write(0, write_zoop, zoop[j].name)
            for p in range(len(zoop[j].Cb)):
                chem_write = p + 1
                worksheet.write(chem_write, write_zoop, zoop[j].Cb[p])
        for j in range (len(fish)):
            write_fish = ((write_region + write_phyto + write_zoop) + j)
            worksheet.write(0, write_fish, fish[j].name)
            for p in range(len(fish[j].Cb)):
                chem_write = p + 1
                worksheet.write(chem_write, write_fish, fish[j].Cb[p])

    workbook.close()
