import xlrd

# Reads in all data from spread sheets

def convert_to_lists(filename):

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
