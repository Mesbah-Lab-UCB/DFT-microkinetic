import pickle
import os
import subprocess
import pandas as pd
import numpy as np
import math
import re
from scipy.interpolate import interp1d
import sys
sys.path.insert(0, '/expanse/lustre/projects/ucb316/mesbah/qtplaskin')
from qtplaskin import FastDirData
 
start = int(sys.argv[1])
SLURM_SUBMIT_DIR = sys.argv[2]
ARRAY_TASK_DIR = sys.argv[3]
def auto_parse_to_kin(id):
    spec_ans = []
    reac_ans = []
    power_ans = []
    
    param_values = np.loadtxt(SLURM_SUBMIT_DIR+'/sobol_vals.txt')
    with open(SLURM_SUBMIT_DIR+'/pred_collect.pkl', 'rb') as fp:
        pred_collect = pickle.load(fp)
    with open(SLURM_SUBMIT_DIR+'/basis.pkl', 'rb') as fp:
        basis = pickle.load(fp)
    with open(SLURM_SUBMIT_DIR+'/entropyPara.pkl', 'rb') as fp:
        entropyPara = pickle.load(fp)
    with open(SLURM_SUBMIT_DIR+'/EleFieldforEntropy.pkl', 'rb') as fp:
        EleFieldforEntropy = pickle.load(fp)
    for i in range(id*10,(id+1)*10):
        curPara = param_values[i,:]
        Temp = curPara[3]
        ED = curPara[0]
        EF = curPara[1]
        N2_frac = curPara[2]
        
        ele = [-1.393915286, -0.407617857, 0., 0.405267196, 0.814412714, 0.081528]
        basis_col = np.zeros((basis.shape[0], 4))
        ele_final_col_500 = np.zeros((1, 4)) # only 4 reactions using reaction energies
        ele_final_col_1000 = np.zeros((1, 4))
        ele_final_col_5000 = np.zeros((1, 4))
        ele_final_col_10000 = np.zeros((1, 4))
        
        for ele_id, cur_ele in enumerate(pred_collect):
            mymodel = np.poly1d(np.polyfit(ele[1:], cur_ele[1:], 3))
            #mymodel(myline)
            if ele_id < 12:
                continue
            else:
                basis_col[:, ele_id-12] = mymodel(basis.reshape(-1))
                # TODO: are they really negative?
                ele_final_col_500[:, ele_id-12] = mymodel(-EF/2)
                ele_final_col_1000[:, ele_id-12] = mymodel(-EF)
                ele_final_col_5000[:, ele_id-12] = mymodel(-EF*5)
                ele_final_col_10000[:, ele_id-12] = mymodel(-EF*10)
                
        basis_entropy = np.array([0.081528/2, 0.081528, 0.081528*5, 0.081528*10]).reshape(-1,1)
        basis_col_entropy = np.zeros((basis.shape[0], entropyPara.shape[0]))
        ele_final_col_500_entropy = np.zeros((1, entropyPara.shape[0]))
        ele_final_col_1000_entropy = np.zeros((1, entropyPara.shape[0]))
        ele_final_col_5000_entropy = np.zeros((1, entropyPara.shape[0]))
        ele_final_col_10000_entropy = np.zeros((1, entropyPara.shape[0]))
        for ele_id, cur_ele in enumerate(entropyPara):
            mymodel = np.poly1d(np.polyfit(EleFieldforEntropy, cur_ele, 3))
            basis_col_entropy[:, ele_id] = mymodel(basis_entropy.reshape(-1))
            ele_final_col_500_entropy[:, ele_id] = mymodel(EF/2)
            ele_final_col_1000_entropy[:, ele_id] = mymodel(EF)
            ele_final_col_5000_entropy[:, ele_id] = mymodel(EF*5)
            ele_final_col_10000_entropy[:, ele_id] = mymodel(EF*10)
        
        # necessary, mapping the magnitude back for ZDPlaskin
        EF = EF*1e10/1e3
        
        subprocess.run('cd ' + ARRAY_TASK_DIR +' && mkdir {}_datastore'.format(i),shell=True)
        
        
        fake_pulse = pd.DataFrame(np.vstack((EF,ED)).T)
        fake_pulse.columns = ["EField", "Electrons_cm-3"]
        fake_pulse.to_csv(ARRAY_TASK_DIR + '/{}_datastore/Ele_{}.dat'.format(i,i), sep=' ',index=False)

        fake_gas_T = pd.DataFrame(np.array([Temp]).reshape(1,-1))
        fake_gas_T.columns = ["K"]
        fake_gas_T.to_csv(ARRAY_TASK_DIR + '/{}_datastore/Tgas_{}.dat'.format(i,i), sep=' ',index=False)

        other_para = pd.DataFrame(np.array([N2_frac]).reshape(1,-1))
        other_para.columns = ["N2_frac"]
        other_para.to_csv(ARRAY_TASK_DIR + '/{}_datastore/other_para_{}.dat'.format(i,i),sep=' ',index=False)

        fake_ac_E = pd.DataFrame(ele_final_col_10000[0, :].reshape(1,-1))
        fake_ac_E.columns = ['H2_NHs_ace', 'Ns_Hs_ACT_E', 'NHs_Hs_ACT_E', 'NH2s_Hs_ACT_E']
        fake_ac_E.to_csv(ARRAY_TASK_DIR + '/{}_datastore/REACTION_E_IN_{}.DAT'.format(i,i), sep=' ',index=False)

        fake_ac_E = pd.DataFrame(basis_col[3, :].reshape(1,-1)) 
        #0 is 500 enhancement, 1 is 1000 enhancement, 2 is 5000, 3 is 10000
        fake_ac_E.columns = ['H2_NHs_ace', 'Ns_Hs_ACT_E', 'NHs_Hs_ACT_E', 'NH2s_Hs_ACT_E']
        fake_ac_E.to_csv(ARRAY_TASK_DIR + '/{}_datastore/REACTION_E_BASIS_{}.DAT'.format(i,i), sep=' ',index=False)

        fake_ac_E = pd.DataFrame(ele_final_col_10000_entropy[0, :].reshape(1,-1))
        fake_ac_E.columns = ['H2_vib_1', 'N2_vib_1', 'NH_vib_1', 'NH2_vib_1', 'NH2_vib_2', 'NH2_vib_3', 'NH3_vib_1', 'NH3_vib_2', 'NH3_vib_3', 'NH3_vib_4', 'NH3_vib_5', 'NH3_vib_6', 'H2_inertia_1', 'N2_inertia_1', 'NH_inertia_1', 'NH2_inertia_1', 'NH2_inertia_2', 'NH2_inertia_3', 'NH3_inertia_1', 'NH3_inertia_2', 'NH3_inertia_3']
        fake_ac_E.to_csv(ARRAY_TASK_DIR + '/{}_datastore/ENTROPY_PARA_IN_{}.DAT'.format(i,i), sep=' ',index=False)

        fake_ac_E = pd.DataFrame(basis_col_entropy[3, :].reshape(1,-1)) 
        #0 is 500 enhancement, 1 is 1000 enhancement, 2 is 5000, 3 is 10000
        fake_ac_E.columns = ['H2_vib_basis_1', 'N2_vib_basis_1', 'NH_vib_basis_1', 'NH2_vib_basis_1', 'NH2_vib_basis_2', 'NH2_vib_basis_3', 'NH3_vib_basis_1', 'NH3_vib_basis_2', 'NH3_vib_basis_3', 'NH3_vib_basis_4', 'NH3_vib_basis_5', 'NH3_vib_basis_6', 'H2_inertia_basis_1', 'N2_inertia_basis_1', 'NH_inertia_basis_1', 'NH2_inertia_basis_1', 'NH2_inertia_basis_2', 'NH2_inertia_basis_3', 'NH3_inertia_basis_1', 'NH3_inertia_basis_2', 'NH3_inertia_basis_3']
        fake_ac_E.to_csv(ARRAY_TASK_DIR + '/{}_datastore/ENTROPY_INFO_BASIS_{}.DAT'.format(i,i), sep=' ',index=False)
        
        with open(SLURM_SUBMIT_DIR + '/bolsigdb.dat') as f:
            lines = f.readlines()

        with open(ARRAY_TASK_DIR + '/{}_datastore/bolsigdb.dat'.format(i), 'w') as f:
            for line in lines:
                f.write(line)
        
        with open(SLURM_SUBMIT_DIR + '/Const_E.F90') as f:
            lines = f.readlines()

        newlines = []
        for line in lines:
            temp_line = re.sub('Ele.dat', 'Ele_{}.dat'.format(i), line)
            temp_line = re.sub('other_para.dat', 'other_para_{}.dat'.format(i), temp_line)
            newlines.append(re.sub('Tgas.dat', 'Tgas_{}.dat'.format(i), temp_line))

        with open(ARRAY_TASK_DIR + '/{}_datastore/Const_E_{}.F90'.format(i,i), 'w') as f:
            for line in newlines:
                f.write(line)
        
        
        with open(SLURM_SUBMIT_DIR + '/zdplaskin_m_DFT_in_entropy_varying_basis.F90') as f:
            lines = f.readlines()

        newlines = []
        for line in lines:
            temp_line = re.sub('REACTION_E_IN.DAT', 'REACTION_E_IN_{}.DAT'.format(i), line)
            temp_line = re.sub('REACTION_E_BASIS.DAT', 'REACTION_E_BASIS_{}.DAT'.format(i), temp_line)
            temp_line = re.sub('ENTROPY_PARA_IN.DAT', 'ENTROPY_PARA_IN_{}.DAT'.format(i), temp_line)
            newlines.append(re.sub('ENTROPY_INFO_BASIS.DAT', 'ENTROPY_INFO_BASIS_{}.DAT'.format(i), temp_line))

        with open(ARRAY_TASK_DIR + '/{}_datastore/zdplaskin_m_DFT_in_entropy_varying_basis_{}.F90'.format(i,i), 'w') as f:
            for line in newlines:
                f.write(line)
        
        subprocess.run('cd ' + ARRAY_TASK_DIR + '/{}_datastore && gfortran -O3 -ffree-line-length-none -o entropyV_extract_{}.exe {}/dvode_f90_m.f90 zdplaskin_m_DFT_in_entropy_varying_basis_{}.F90 Const_E_{}.F90 {}/bolsig_x86_64.so'.format(i,i,SLURM_SUBMIT_DIR,i,i,SLURM_SUBMIT_DIR), shell=True)
        try:
            subprocess.run('cd ' + ARRAY_TASK_DIR + '/{}_datastore && export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH && ./entropyV_extract_{}.exe'.format(i,i), shell=True, timeout=60*60, stdout=open(os.devnull, "w"))
        except:
            pass
        
        
        d = FastDirData(ARRAY_TASK_DIR + '/{}_datastore'.format(i))
        if d.t[-1]<4000.:
            spec_ans.append([0]*len(d.species))
            reac_ans.append([0]*len(d.reactions))
            power_ans.append(0.)
        else:
            tmp_spec = []
            for s in d.species:
                approx = interp1d(d.t, d.get_spec(s))
                tmp_spec.append(approx(4000))
            tmp_reac = []
            for r in d.reactions:
                approx = interp1d(d.t, d.get_rate(r))
                tmp_reac.append(approx(4000))
            spec_ans.append(tmp_spec[:])
            reac_ans.append(tmp_reac[:])
            approx = interp1d(d.t, d.get_cond('power'))
            power_ans.append(approx(4000))
    with open(SLURM_SUBMIT_DIR+'/saved_res_{}'.format(id), 'wb') as fp:
        pickle.dump((spec_ans, reac_ans, power_ans), fp)
        
# Need a wrapper function because map() only operates on one argument
def wrapper(i):
    return(auto_parse_to_kin(i))

#n = len(range(start, end))
import time
time.time()
wrapper(start) # Run calculation in parallel
time.time()