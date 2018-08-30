from datetime import datetime

def num_steps(start_string, end_string, step):

    start = datetime.strptime(start_string, '%m,%d,%Y').date()
    end = datetime.strptime(end_string, '%m,%d,%Y').date()

    if step == 'Week':
        difference = (end-start).days//7
    if step == 'Month':
        difference = (end-start).days//30
    if step == 'Quarter':
        difference = (end - start).days // 91

    return difference

