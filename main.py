import Bioaccum
import FR_Input_Output
import spatial

# are we solving steady state on single region, solving with time on single region, or time on multiple regions
def filter_cases(filename):

    model_para, all_data, time_steps, time_per_step, site_data = FR_Input_Output.convert_to_lists(filename)

    # # set up for
    # u_iter = int(model_para[0])
    # v_iter = int(model_para[1])
    # sample_data = model_para[0:3]
    # Bioaccum.set_all_h_and_s(sample_data, all_data)
    b, r, h =spatial.setup(site_data)
    spatial.get_location(b,r,h, 'Pumpkinseed', 100)
    # # if steady state problem wants to be solved
    # if len(all_data[0]) == 1 and model_para[8] == 'YES':
    #     dictionaries = Bioaccum.bio_monte_carlo_loop(model_para, all_data, all_data[1][0][0], 0, u_iter, v_iter, 0)
    #     results_dic= Bioaccum.filter_stat_case(dictionaries, filename)
    #
    # if model_para[8] == 'NO':
    #
    #     prior_concentrations = Bioaccum.init_prior_con_dic(all_data[0], all_data[2], all_data[3], all_data[4], all_data[5])
    #     t = 0
    #     steps_dic = []
    #     # print(time_steps, time_per_step)
    #     while t < time_steps:
    #         region_index = 0
    #         temp = all_data[1][region_index][t]
    #         days = t * time_per_step
    #         dictionaries = Bioaccum.bio_monte_carlo_loop(model_para, all_data, temp, region_index, u_iter, v_iter,
    #                                             days, per_step=time_per_step, p_dic=prior_concentrations)
    #
    #         steps_dic.append(dictionaries)
    #         # Need to get prior concentrations, in none statistical case this is the only dictionary
    #         if len(dictionaries) == 1:
    #             prior_concentrations = dictionaries[0]
    #         # we pass all of the dictonaries, and have to iterate through them according to u_count, and v_count
    #         else:
    #             prior_concentrations = dictionaries
    #
    #         t += 1
    #
    #     results_dic = Bioaccum.filter_stat_case(dictionaries)


filter_cases('sheets/input/tests/testy_test.xlsx')