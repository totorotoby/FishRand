import xlrd, xlsxwriter
import networkx as nx
import prob as pr
import Time_parser
import numpy as np
import math
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
            try:
                data_get_helper(row[j], r_temps)
            except IndexError:
                print("Not enough tempature data for the number of timesteps you have input")
                exit(0)
            
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

    try:
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
                        
                        
                        data_get_helper(column[(1+j) + ((j * num_chem) + (1+k)) + ((j * num_chem * 4) + (k * 4) + p)], c_cons)
                    r_cons.append(c_cons)
                t_cons.append(r_cons)
            con_data.append(t_cons)
    except:
        print("Something is wrong with your chemical concentration entries.\nCheck the format on Page 13 section 3.7 of the user manual.")
        exit(0)
            

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
        try:
            site = (site_col_name[row].value, [float(i) for i in (site_col_coord[row].value.replace('', '').split(', '))])
        except ValueError:
            print('Some coordinate entries are not formatted correctly in Sample Sites.\nRefer to page 14 section 3.9 of the user manual for help.')
            exit(0)
        sites.append(site)
        row += 1

    hotspots = []
    row_num = 40
    row = sites_sheet.row(row_num)
    while row[0].value != '':
        col = 3
        row = sites_sheet.row(row_num)
        defin = []
        while row[col].value != '':
            if row[col].value == 'END':
                break
            defin.append(sites_sheet.row(row_num)[col].value)
            col += 1

        hotspot = HotSpot(sites_sheet.row(row_num)[0].value, 'Polygon',
                          sites_sheet.row(row_num)[1].value, sites_sheet.row(row_num)[2].value, defin)
        if len(hotspot.definition) != 0:
            hotspots.append(hotspot)
        row_num += 1

    draws = int(sites_sheet.row(38)[1].value)
    print(draws)
    return [boundary, sites, hotspots, draws]


def write_output_steady(total_cons, output_name, stop, dist_type, inputs, data, tofit):


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
    if tofit == 1:
        sheet.write(0,0, 'Region: ' + str(list(total_cons.keys())[0]) + ' Concentrations (ng/g ww)', bold)

    #write chemicals at the top
    for j in range(chem_len):
        sheet.write(0, j+1, chem_list[j][0])
    #write orgs on the column
    for i in range(org_len):
        sheet.write(i+1, 0, org_list[i])


    # if statstical
    if type(chem_list[0][1]) == pr.ResultDist:
        if tofit == 1:
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
        else:
            sheet.write(0, 0, 'Unfitted Mean and Std of Simulations (ng/g ww)', bold)
            for i in range(org_len):
                for j in range(chem_len):
                    to_print = total_cons[region][org_list[i]][chem_list[j][0]].v_mean_stdString
                    sheet.write(i + 1, j + 1, to_print)

    # if nonstatstical
    else:
        for i in range(org_len):
            for j in range(chem_len):
                if type(total_cons[region][org_list[i]][chem_list[j][0]]) == float:
                    sheet.write(i+1,j+1, round(total_cons[region][org_list[i]][chem_list[j][0]],6))
                else:
                    sheet.write(i+1,j+1, round(total_cons[region][org_list[i]][chem_list[j][0]][0],6))

    workbook.close()
    print('Output written to: '+ output_name)


# returns the weighted standard devations over animals, weighted by area of regions.
def get_weighted_stds(v_dic, areas):

    weights = [areas[i]/sum(areas) for i in range(len(areas))]
    
    valuesArray = []
    meansArray = []
    Arrayify(v_dic, valuesArray)

    for j in range(len(valuesArray[0])):
        animal = []
        for k in range(len(valuesArray[0][0])):
            mean = 0
            for i in range(len(valuesArray)):
                mean += valuesArray[i][j][k] * weights[i]
            animal.append(round(mean,6))
        meansArray.append(animal)
    
    stds = []
    for j in range(len(valuesArray[0])):
        animal = []
        for k in range(len(valuesArray[0][0])):
            std = 0
            for i in range(len(valuesArray)):
                std += round(math.sqrt(weights[i] *(valuesArray[i][j][k] - meansArray[j][k])**2), 6)
            animal.append(std)
        stds.append(animal)
        
    return stds, meansArray

# Function takes in a nested dictonary of arbitrary depth, drops all the keys
# and turns the values into a nested dictonary.
def Arrayify(D, L):

    if type(list(D.items())[0][1]) == float or  type(list(D.items())[0][1]) == np.float64:
        
        v = list(D.items())
        for pair in v:
            L.append(pair[1])   
    else:
        count = 0
        for subD in list(D.values()):
            L.append([])
            Arrayify(subD, L[count])
            count += 1





def makeRegionTab(workbook, rAreas, rInfo):

    rSheet = workbook.add_worksheet('Regional Info')
    big = workbook.add_format({'font_size': 14, 'bold': True})

    # list of coordnites of interior regions
    rCoords = [list(r[1].exterior.coords) for r in rInfo[1]]
    rAreas =  [r[1].area for r in rInfo[1]]
    intNames = [r[0] for r in rInfo[1]]
    # max number of sides of interior regions
    maxNCoor = max([len(rCoord) for rCoord in rCoords])
    # list of boundry coordinates
    bCoords = list(rInfo[0].exterior.coords)
    bArea = rInfo[0].area
    # max number of sides over boundary and interior
    maxLen = max([len(bCoords), maxNCoor])
    
    rSheet.write(0,0, 'Coordinates and Area of Thessian Polygon Generated Regions (x,y)', big)

    # Labels
    for i in range(maxLen):
        rSheet.write(1, i+1, "Coord " + str(i))

    rSheet.write(1, maxLen + 1, "Area")

    # Writing Boundary
    rSheet.write(2, 0, 'Boundary')
    for j in range(len(bCoords)):
        rSheet.write(2, j + 1, str(bCoords[j]).replace('(', '').replace(')', ''))
        
    rSheet.write(2, maxLen + 1, bArea)

    # Writing Interior regions
    for i in range(len(intNames)):
        rSheet.write(3 + i, 0, str(intNames[i]))
        for j in range(len(rCoords[i])):
            rSheet.write(3 + i, j + 1, str(rCoords[i][j]).replace('(', '').replace(')', ''))
        rSheet.write(3+i, maxLen + 1, rAreas[i])
            
    
    
    
    
        

def writeInputsTab(workbook, inputs, data):

   
    regions = data[0]
    chemicals = data[2]
    phyto = data[3][0]
    zoop = data[4][0]
    nonMove = data[5]
    move = data[6]
    chemCon = data[len(data)-1]
    phytoA = inputs[2][0][0].A
    phytoB = inputs[2][0][0].B

    sheet = workbook.add_worksheet("Chemical and Biological Inputs")
    big = workbook.add_format({'font_size': 14, 'bold': True})
    rCount = 1
    sheet.write(0, 0, 'Regional inputs', big)
    sheet.write(1, 0, 'Dissolved Organic Carbon Content (kg/L)')
    sheet.write(2, 0, 'Particulate Organic Carbon Content (kg/L)')
    sheet.write(3, 0, 'Concentration Suspended Solids (g/L)')
    sheet.write(4, 0, 'Fraction Organic Carbon Content in Sediment')
    sheet.write(5, 0, 'Disolved Oxygen Saturation')
    sheet.write(6, 0, 'Disolved Oxygen Concentration (mg O2/L)')
    sheet.write(7, 0, 'DOC octanol Proportionality Constant')
    sheet.write(8, 0, 'POC octanol Proportionality Constant')
    
    for region in regions:
        for i in range(len(region)):
            if i in [1, 2, 3, 4, 6]:
                if region[i] == '':
                  sheet.write(i, rCount, 'DNE')
            if i == 7 and region[i] == '':
                sheet.write(i, rCount, .08)
            if i == 8 and region[i] == '':
                sheet.write(i, rCount, .35)
                  
            sheet.write(i, rCount, str(region[i]))
            
        rCount = rCount + 2

    lenR = 9
    
    sheet.write(9, 0, 'Chemical Inputs', big)
    sheet.write(10, 0, 'Kow')    
    sheet.write(11, 0, 'Disequilibrium factor of Disolved Organic Carbon')
    sheet.write(12, 0, 'Disequilibrium factor of Particualte Organic Carbon')
    sheet.write(13, 0, 'Beta Lipid Omn')
    sheet.write(14, 0, 'Beta Non-Lipid Om')
    sheet.write(15, 0, 'Beta Lipid digesta')
    sheet.write(16, 0, 'Beta Sediment OC')
    sheet.write(17, 0, 'Beta Non-lipid digesta')

    cCount = 1

    for chemical in chemicals:
        for i in range(len(chemical)):
            if i in [2,3] and chemical[i] == '':
                sheet.write(i+lenR, cCount, 'DNE')
            if i == 4 and chemical[i] == '':
                sheet.write(i+lenR, cCount, 1)
            if i == 5 and chemical[i] == '':
                sheet.write(i+lenR, cCount, .035)
            if i == 6 and chemical[i] == '':
                sheet.write(i+lenR, cCount, 1)
            if i == 7 and chemical[i] == '':
                sheet.write(i+lenR, cCount, .35)
            if i == 8 and chemical[i] == '':
                sheet.write(i+lenR, cCount, .035)
            sheet.write(i+lenR, cCount, str(chemical[i]))
        cCount = cCount + 2

        
    sheet.write(18, 0, 'Phyto Inputs', big)
    sheet.write(19, 0, 'Growth Rate (1/days)')
    sheet.write(20, 0, 'Lipid Content (kg/kg)')
    sheet.write(21, 0, 'Non-Lipid Content (kg/kg)')
    sheet.write(22, 0, 'Constant A (see equation 10 Arnot 2004) (Time)')
    sheet.write(23, 0, 'Constant B (see equation 10 Arnot 2004) (Time)')

    RClen = 18
    
    for i in range(len(phyto)):
        if i == 1 and phyto[i] == '':
            sheet.write(RClen + i, 1, .08)
        if i == 2 and phyto[i] == '':    
            sheet.write(RClen + i, 1, .005)
        if i == 3 and phyto[i] == '':
            sheet.write(RClen + i, 1, .065)
        sheet.write(RClen + i, 1, str(phyto[i]))
    sheet.write(22, 1, phytoA)
    sheet.write(23, 1, phytoB)

    sheet.write(24, 0, 'Zoops Inputs', big)
    sheet.write(25, 0, 'Weight (Kg)')
    sheet.write(26, 0, 'Lipid Content (Kg/Kg)')
    sheet.write(27, 0, 'Non-Lipid Content (Kg/Kg)')
    sheet.write(28, 0, 'Fraction Pore Water Ventilated')
    sheet.write(29, 0, 'Dietary absorption efficency of lipid OM')
    sheet.write(30, 0, 'Dietary absorption efficency of non-lipid OM')
    sheet.write(31, 0, 'Dietary absorption efficency of water')
    sheet.write(32, 0, 'Feeding Rate (Kg/day)')
    sheet.write(33, 0, 'Growth Rate (1/Day)')
    sheet.write(34, 0, 'Filter Feeder')

    RCPlen = 24
    
    for i in range(len(zoop)):
        if i == 3 and zoop[i] == '':
            sheet.write(RCPlen + i, 1, .2)
        if i == 4 and zoop[i] == '':
            sheet.write(RCPlen + i, 1, 0)
        if i == 5 and zoop[i] == '':
            sheet.write(RCPlen + i, 1, .72)
        if i == 6 and zoop[i] == '':
            sheet.write(RCPlen + i, 1, .72)
        if i == 7 and zoop[i] == '':    
            sheet.write(RCPlen + i, 1, .25)
        if i in [8,9] and zoop[i] == '':
            sheet.write(i+RCPlen, cCount, 'DNE')
        sheet.write(i+RCPlen, 1, str(zoop[i]))
        
        if i == 10:
            if zoop[i] == 1:
                sheet.write(RCPlen + i, 1, 'Yes')
            else:
                sheet.write(RCPlen + i, 1, 'No')


            
    sheet.write(35, 0, 'Non-moving Organisms Inputs', big)
    sheet.write(36, 0, 'Weight (Kg)')
    sheet.write(37, 0, 'Lipid Content (Kg/Kg)')
    sheet.write(38, 0, 'Non-Lipid Content (Kg/Kg)')
    sheet.write(39, 0, 'Fraction Pore Water Ventilated')
    sheet.write(40, 0, 'Dietary absorption efficency of lipid OM')
    sheet.write(41, 0, 'Dietary absorption efficency of non-lipid OM')
    sheet.write(42, 0, 'Dietary absorption efficency of water')
    sheet.write(43, 0, 'Feeding Rate (Kg/day)')
    sheet.write(44, 0, 'Growth Rate (1/day)')
    sheet.write(45, 0, 'Filter Feeder')

    s = 35
    nmCount = 1
    
    for nonmove in nonMove:
        for i in range(len(nonmove)):
            if i == 3 and nonmove[i] == '':
                sheet.write(s + i, nmCount, .2)
            if i == 4 and nonmove[i] == '':
                sheet.write(s + i, nmCount, 0)
            if i == 5 and nonmove[i] == '':
                sheet.write(s + i, nmCount, .72)
            if i == 6 and nonmove[i] == '':
                sheet.write(s + i, nmCount, .72)
            if i == 7 and nonmove[i] == '':    
                sheet.write(s + i, nmCount, .25)
            if i in [8,9] and nonmove[i] == '':
                sheet.write(i+s, nmCount, 'DNE')
                
            sheet.write(i+s, nmCount, str(nonmove[i]))
            
            if i == 10:
                if nonmove[i] == 1:
                    sheet.write(s + i, nmCount, 'Yes')
                else:
                    sheet.write(s + i, nmCount, 'No')
                
        nmCount = nmCount + 2


    sheet.write(46, 0, 'Moving Organisms Inputs', big)
    sheet.write(47, 0, 'Weight (Kg)')
    sheet.write(48, 0, 'Lipid Content (Kg/Kg)')
    sheet.write(49, 0, 'Non-Lipid Content (Kg/Kg)')
    sheet.write(50, 0, 'Fraction Pore Water Ventilated')
    sheet.write(51, 0, 'Dietary absorption efficency of lipid OM')
    sheet.write(52, 0, 'Dietary absorption efficency of non-lipid OM')
    sheet.write(53, 0, 'Dietary absorption efficency of water')
    sheet.write(54, 0, 'Feeding Rate (Kg/Day)')
    sheet.write(55, 0, 'Growth Rate (1/Day)')
    sheet.write(56, 0, 'Filter Feeder')

    st = 46
    
    mCount = 1
    for m in move:
        for i in range(len(m)):
            if i == 3 and m[i] == '':
                sheet.write(st + i, mCount, .2)
            if i == 4 and m[i] == '':
                sheet.write(st + i, mCount, 0)
            if i == 5 and m[i] == '':
                sheet.write(st + i, mCount, .75)
            if i == 6 and m[i] == '':
                sheet.write(st + i, mCount, .75)
            if i == 7 and m[i] == '':    
                sheet.write(st + i, mCount, .5)
            if i in [8,9] and m[i] == '':
                sheet.write(i + st, mCount, 'DNE')
                
            sheet.write(i + st, mCount, str(m[i]))

            if i == 10:
                if m[i] == 1:
                    sheet.write(st + i, mCount, 'Yes')
                else:
                    sheet.write(st + i, mCount, 'No')
        mCount = mCount + 2
        
                     
    workbook.close()
    
def write_temporal_excel(array, output_name, stops, stat_flag, regional_info, dist_type, inputs, data, tofit=None):

    
    types = ["Normal", 'Lognormal', 'Uniform', 'Gamma', "KS Best"]
    index = types.index(dist_type)


    workbook = xlsxwriter.Workbook(output_name)
    bold = workbook.add_format({'bold': True})
    big = workbook.add_format({'font_size': 14, 'bold': True})
    regional_areas = [region[1].area for region in regional_info[1]]
    
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
            sheet = workbook.add_worksheet('Timestep ' + str(stops[i]+1))
            data_at_time = array[i]
            lower_non_avg = data_at_time[0][0]
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

            stds_lower, mean_lower = get_weighted_stds(lower_non_avg, regional_areas)
            sheet.write(0, len(chem_list) + 1, 'Lower Food Web Mean Concentrations and Standard Deviations Weighted by Regional Area (ng/g ww)', big)
            for j in range(len(chem_list)):
                sheet.write(0, len(chem_list) + 1 + (j + 1), chem_list[j])
            for i in range(len(lower_org_list)):
                for j in range(len(chem_list)):
                    sheet.write(2 + i, 2 + len(chem_list) + j, str(mean_lower[i][j]) + ', ' + str(stds_lower[i][j]) )


            sheet.write(count + 1, 0, 'Upper Food Web Concentrations Averaged over Populations (ng/g ww)', big)
            for i in range(len(upper_org_list)):
                sheet.write(count + 2 + i, 0, upper_org_list[i])
                for j in range (len(chem_list)):
                    sheet.write(count + 2 + i, j+1, str(round(upper_avg[upper_org_list[i]][chem_list[j]][0],6)) + ', ' + str(round(upper_avg[upper_org_list[i]][chem_list[j]][1],6)))


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

            if tofit == 1:
                sheet.write(0, 0, 'Lower Food Web Concentrations by Region (ng/g ww)', big)
            else:
                sheet.write(0, 0, 'Lower Food Web Concentrations by Region (ng/g ww) (Mean, Standard Deviation)', big)

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
                        if tofit == 1:
                            if index != 4:
                                lower_dists[reg_list[i]][lower_org_list[j]][chem_list[k]].index = index
                                toprint = lower_dists[reg_list[i]][lower_org_list[j]][chem_list[k]].bestparam()
                                sheet.write((1 + i * (1 + len(lower_org_list))) + (1+j), k + 1, toprint[0])
                            else:
                                sheet.write((1 + i * (1 + len(lower_org_list))) + (1+j), k + 1,
                                            lower_dists[reg_list[i]][lower_org_list[j]][chem_list[k]].best_para[0])
                        else:
                            toprint = lower_dists[reg_list[i]][lower_org_list[j]][chem_list[k]].v_mean_stdString
                            
                            sheet.write((1 + i * (1 + len(lower_org_list))) + (1+j), k + 1, toprint)


            sheet.write(0, len(chem_list) + 1, 'Lower Food Web Mean Concentrations and Standard Deviations Weighted by Regional Area (ng/g ww)', big)
            sheet.write(1, len(chem_list) + 1, 'All Regions', bold)

            for j in range(len(chem_list)):
                sheet.write(0, len(chem_list) + 1 + (j + 1), chem_list[j])


            for i in range(len(chem_list)):
                for j in range(len(lower_org_list)):
                    weight_avg = 0
                    weight_std = 0
                    rStds = [lower_dists[reg_list[k]][lower_org_list[j]][chem_list[i]].v_mean_std[1] for k in range(len(reg_list))]
                    for k in range(len(reg_list)):
                        weight_avg += normed_reg_areas[k]*lower_dists[reg_list[k]][lower_org_list[j]][chem_list[i]].v_mean_std[0]
                       #weight_std += normed_reg_areas[k]*lower_dists[reg_list[k]][lower_org_list[j]][chem_list[i]].v_mean_std[1]
                    stdAllRegions(rStds, weight_avg, weights) 
                       
                    sheet.write(2 + j, 2 + len(chem_list) + i, str(round(weight_avg, 6)) + ', ' + str(round(weight_std, 6)))



            if tofit == 1:
                sheet.write(count + 1, 0, 'Upper Food Web Concentrations (ng/g ww)', big)
            else:
                sheet.write(count + 1, 0, 'Upper Food Web Concentrations (ng/g ww) (Mean, Standard Deviation)', big)
            count += 1

            for i in range(len(upper_org_list)):
                sheet.write(count + 1 + i, 0,upper_org_list[i])
                for j in range(len(chem_list)):
                    if tofit == 1:
                        if index != 4:
                            upper_dists[upper_org_list[i]][chem_list[j]].index = index
                            toprint = upper_dists[upper_org_list[i]][chem_list[j]].bestparam()
                            sheet.write(count + 1 + i, j+1, toprint[0])
                        else:
                            sheet.write(count + 1 + i, j + 1, upper_dists[upper_org_list[i]][chem_list[j]].best_para[0])
                    else:
                        toprint = upper_dists[upper_org_list[i]][chem_list[j]].v_mean_stdString
                        sheet.write(count + 1 + i, j+1, toprint)

            if tofit == 1:
                sheet.write(count, 1 + len(chem_list), 'Upper Food Web Concentrations (ng/g ww) (Mean and Standard Deviation of Samples)', big)

                for i in range(len(upper_org_list)):
                    for j in range(len(chem_list)):
                        sheet.write(count + 1 + i, j+ len(chem_list) + 2, str(upper_dists[upper_org_list[i]][chem_list[j]].v_mean_std))


    makeRegionTab(workbook, regional_areas, regional_info)
    writeInputsTab(workbook, inputs, data)
    
        

    print('Output written to: '+ output_name)
   

def stdAllRegions(rStds, wMeanGroup, weights):

    print(rStds)
    print(wMeanGroup)
    print(weigths)







    
