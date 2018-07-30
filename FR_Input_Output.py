import xlrd
import xlsxwriter

# Reads in all data from spread sheets

#def stat_convert_to_lists(filename):




def determ_convert_to_lists(filename):

    fish_data_raw = []  # messy output from xlrd
    fish_data = []  # clean fish array
    zoo_data = []  # clean zoo array
    phyto_data = []  # clean phyto array
    f_len = 11  # number of inputs per fish
    zo_len = 11  # number of inputs per zoo
    ph_len = 4  # number of inputs per phyto
    reg_chem_len = [10,7] # number of inputs per region

    # getting columns #
    all_sheets = xlrd.open_workbook(filename)

    org_sheet = all_sheets.sheet_by_index(2)

    fish_column = org_sheet.col(1)
    zoop_column = org_sheet.col(3)
    phy_column = org_sheet.col(5)

    # getting region and chemical data #
    chem_reg_data=[]
    for k in range (2):
        data = []
        sheet = all_sheets.sheet_by_index(k)
        reg_column = sheet.col(1)
        data_raw = []
        for cell in reg_column:
            data_raw.append(cell.value)

        new_reg = []
        for i in range(len(data_raw)):
            # if we aren't on a label
            if i % (reg_chem_len[k] + 1) == 0:
                # if the list has something in it
                if len(new_reg) != 0:
                    data.append(new_reg)
                new_reg = []
            else:
                new_reg.append(data_raw[i])
        if len(new_reg) != 0:
            data.append(new_reg)
        chem_reg_data.append(data)

#!!!!!!!!!!! might need to change so we can have multiple phyto !!!!!!!!!!!

    # getting phyto and zoop data #
    for i in range (zo_len+1):
        if i > 0:
            zoo_data.append(zoop_column[i].value)
    for i in range (ph_len+1):
        if i > 0:
            phyto_data.append(phy_column[i].value)


# !!!!!!!!!!! might need to change so we can have multiple phyto !!!!!!!!!!!

    # getting fish_data #

    for cell in fish_column:
        fish_data_raw.append(cell.value)

    new_fish = []
    for i in range(len(fish_data_raw)):
        # if we aren't on a label
        if i%(f_len+1) == 0:
            # if the list has something in it
            if len(new_fish) != 0:
                fish_data.append(new_fish)
            new_fish = []
        else:

            new_fish.append(fish_data_raw[i])
    fish_data.append(new_fish)


    # getting Diet data #

    diet_data = {}
    entrysize = len(fish_data) + 5
    diet_sheet = all_sheets.sheet_by_index(3)
    rows = diet_sheet.nrows


    for i in range(rows):
        if i%entrysize == 0:
            continue
        if i%entrysize == 1:
            new_name = diet_sheet.row(i)[1].value
            diet_data[new_name] = []
        else:
            new_prey_data = [diet_sheet.row(i)[1].value,diet_sheet.row(i)[2].value]
            diet_data[new_name].append(new_prey_data)


    reg_data = chem_reg_data[0]
    chem_data = chem_reg_data[1]

    return reg_data, chem_data, fish_data, zoo_data, phyto_data, diet_data


def deter_write_output(regions, fish, chemicals, phyto, zoop):

    workbook = xlsxwriter.Workbook('hello.xlsx')
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

                print('hello')
                chem_write = p + 1

                worksheet.write(chem_write, write_fish, fish[j].Cb[p])
