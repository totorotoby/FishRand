import sys
sys.path.append("../")
from Bioaccum import *
from numpy import linspace
import openpyxl
import xlrd
import matplotlib.pyplot as plt
import numpy


def reprint(write_work, sheet, write_sheet, dist_col_num, instance_len, scale):

    dist_col = sheet.col(dist_col_num)
    dist_name_col = sheet.col(dist_col_num-1)

    for i in range(len(dist_col)):
            if i % (instance_len + 1) == 0:
                name = dist_col[i].value
                if name == 'END':
                    break
            else:
                if i % (instance_len + 1) == 1:
                    continue
                else:
                    dist_par = dist_col[i].value
                    if dist_par != '':
                        try:
                            mean = dist_par.split(', ')[0]
                        except:
                            continue
                        mean = float(mean)
                        new_std = mean * scale
                        new_print = str(mean) + ', ' + str(new_std)
                        write_sheet.cell(row=i+1, column=dist_col_num+1).value = new_print
    write_work.save('../sheets/input/tests/FR_Input_st_small_Var.xlsx')


def single_round(it, steps):
    f_len = 11  # number of inputs per fish
    zo_len = 11  # number of inputs per zoo
    ph_len = 4  # number of inputs per phyto
    reg_len = 10  # number of inputs per region
    chem_len = 7

    all_sheets = xlrd.open_workbook('../sheets/input/tests/approach_mean_read.xlsx')
    write_work = openpyxl.load_workbook('../sheets/input/tests/FR_Input_st_small_Var.xlsx')
    reg_sheet_read = all_sheets.sheet_by_index(1)
    reg_sheet_write = write_work.worksheets[1]

    chem_sheet_read = all_sheets.sheet_by_index(2)
    chem_sheet_write = write_work.worksheets[2]

    an_sheet_read = all_sheets.sheet_by_index(3)
    an_sheet_write = write_work.worksheets[3]

    reprint(write_work, reg_sheet_read, reg_sheet_write,2,reg_len, steps[it])
    reprint(write_work, chem_sheet_read, chem_sheet_write, 2, chem_len, steps[it])
    reprint(write_work, an_sheet_read, an_sheet_write, 2, f_len, steps[it])
    reprint(write_work, an_sheet_read, an_sheet_write, 5, zo_len, steps[it])
    reprint(write_work, an_sheet_read, an_sheet_write, 9, ph_len, steps[it])


    single_dic = run_bio(1,'../sheets/input/tests/FR_Input_st_small_Var.xlsx', 'approach_mean_write')

    return single_dic

def loop():

    iter = 1
    det_results = run_bio(0,'../sheets/input/FR_Input_det.xls', 'FR_Input_det')
    steps = linspace(.01, .7, iter)
    result_dic = {}
    print('Fraction done:')
    for i in range (iter):
        print('\r' + str(i/iter), end='')
        single_dic = single_round(i,steps)

        for reg, values in single_dic.items():
            if i == 0:
                result_dic[reg] = {}
            for animal, values1 in values.items():
                if i == 0:
                    result_dic[reg][animal] = {}
                for chem in values1.keys():
                    conver_dic(result_dic,single_dic,reg,animal,chem, i)


    print_approach_g(steps,result_dic,det_results)


def print_approach_g(steps,result_dic, det_results):

    for reg, values in result_dic.items():
        for animal, values1 in values.items():
            for chemical, list in values1.items():
                y_label = '(ng/g) of ' + chemical + ' in ' + animal
                x_label = 'coefficient of variation'
                plt.xlabel(x_label)
                plt.ylabel(y_label)
                plt.plot(steps,list, 'ro', color='b')
                plt.axhline(det_results[reg][animal][chemical], color='r')
                plt.show()


def conver_dic(result_dic, single_dic,reg,animal,chem, i):

    if i == 0:
        result_dic[reg][animal][chem] = []
    values = single_dic[reg][animal][chem].values
    mean = numpy.mean(values)
    result_dic[reg][animal][chem].append(mean)

    return result_dic

loop()
