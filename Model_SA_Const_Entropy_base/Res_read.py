import os
import subprocess
import pandas as pd
import numpy as np
import math
import re

import sys
import pickle
#TARGET_DIR = int(sys.argv[1])
start = int(sys.argv[1])
end = int(sys.argv[2])

spec_ans = []
reac_ans = []
power_ans = []
for _ in range(start, end):
    print('I am {}.'.format(_))
    with open('./saved_res_{}'.format(_), 'rb') as fp:
        spec_ans_tmp, reac_ans_tmp, power_ans_tmp = pickle.load(fp)
    spec_ans += spec_ans_tmp
    reac_ans += reac_ans_tmp
    power_ans += power_ans_tmp
with open('./total_rea', 'wb') as fp:
    pickle.dump((spec_ans, reac_ans, power_ans), fp)
