# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 14:44:09 2017

@author: Issun
"""

import math as m
import numpy as np
import matplotlib.pyplot as plt
import datetime
import os
import json
import scipy.stats as stats

def TheoClineFcn(s, sigma):
    return lambda x: (-2 + 3*(m.tanh(x*m.sqrt(s/2)/sigma + m.atanh(m.sqrt(2/3)))**2))/2 + .5

def PlotAndSave(xarrays,yarrays,xlabel,ylabel,title,savename,labels):
   
    for i in range(len(yarrays)):
        if i in [2, 3]:
            plt.plot(xarrays[i], yarrays[i][0:len(xarrays[i])], '-.',\
                     c = plt.cm.rainbow(1/len(yarrays)), label = labels[i], alpha = .8)
            
        else:
            plt.plot(xarrays[i], yarrays[i][0:len(xarrays[i])],\
                     c = plt.cm.rainbow(i/len(yarrays)), label = labels[i])
        
        
        #legend_line.append(mlines.Line2D([], [], color=plt.cm.rainbow(i/len(yarrays)), label=r'$S_%d$' % i))
    
    # plt.legend(handles=legend_line)
    plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.subplots_adjust(left=0.1, right=.75, top=0.9, bottom=0.1)
    plt.savefig(savename, format='eps', dpi=900)
    plt.show()
    plt.clf()
    
    
# Constants    
computername = os.environ['COMPUTERNAME']
filename = '0119_output_' + computername + '.txt'

# Parameters
K = K0 = 150
N = N0 = 100

#s_v = [10**(-i/6) for i in range(7, 10)]
s_v = [1e-3]
sigma = sigma0 = .5

realizations_v = [20 for i in range(len(s_v))]  # 20
n_w_v = [10 for i in range(len(s_v))] # 2000
delay = 0

len_H_v = 1 # !!!!!!!!!!!!!!!!!!!!!!!!
len_H_v_by_2 = int(len_H_v/2)
resolution = 20
bootstrap_n = 100

plot = 0
savedata = 0
deterministic_start = 0
tm = 1

rescale = 0
rescale_factor = 4

mu = 1e-5

tol_equil = 1e-1 # 1e-5?

# Init

if rescale:
    
    sigma = sigma0/rescale_factor
    N = int(rescale_factor * N0)
    K = int(K0/rescale_factor)

with open(filename, 'a') as f:
    f.write('\n')
    f.write(str(datetime.datetime.now()))
    f.write('\n')
    f.write('s,K,N,sigma,<w(p)>,w(<p>),w_theo,n_real,<F>,iter\n')

for s_i, s in enumerate(s_v):
    
    n_w = n_w_v[s_i]
    realizations = realizations_v[s_i]
    
    print('s = %.4f, sigma = %.3f, K = %d, N = %d' % (s, sigma, K, N))
    
    env_change = int(K/2)

    lamb = K * N * 2 * mu
    
    fit_v = [(1-s)**2, 1-s, 1]
    
    # n_w = int(40 * 2e-2 / s)
    
    w_pq_vv = []
    freqs_vv = []
    freq_means_v = []
    w_pq_mean_v = []

    F_v = []
    
    
    for r_i in range(realizations):
        print('realization no. %d' % (r_i+1))

        # Initialization of realization
        if not deterministic_start:
            alleles = [[[0, 0] for i in range(N)] if j < env_change else [[1, 1] for i in range(N)] for j in range(K)]
        else:
            alleles = [[] for deme in range(K)]
            cline = TheoClineFcn(s, sigma)
            theo_cline_values = [((env_change <= j)*2 - 1)*cline(abs(env_change -.5 - j)) + (env_change > j) for j in range(K)]
            for deme in range(K):
                for n in range(N):
                    if n < N * theo_cline_values[deme]:
                        alleles[deme].append([1, 1])
                    else:
                        alleles[deme].append([0, 0])
        w_pq_v = []
        w_iter = 0

        H_v = []

        freqs_v = []
        # F_v = []

        iteration = 0
        equil = done = False

        while not done:
            
            # Migration
            
            moves = np.random.normal(0, sigma, K * N)
            
            alleles_new = [[] for j in range(K)]
            
            for deme in range(K):
                
                for i, ind in enumerate(alleles[deme]):
                    
                    # ind = alleles[deme][i]
                    move = int(round(moves[deme*N + i]))
                    # move = moves[deme*N + i]
                    # move = floor(move) + (move%1 > np.random.rand())
                    
                    if deme + move < 0:
                        new_loc = 0
                    elif deme + move >= K:
                        new_loc = K - 1
                    else:
                        new_loc = deme + move
                    
                    alleles_new[new_loc].append(ind)
        
            # Mating and selection
            
            for deme, inds_in_deme in enumerate(alleles_new):
                
                # inds_in_deme = alleles_new[deme]
        
                n_ind = len(inds_in_deme)
                
                if n_ind:
            
                    if deme < env_change: # <
                        fitness_v = (1-s)**np.sum(inds_in_deme, 1)
                    else:
                        fitness_v = (1-s)**(2-np.sum(inds_in_deme, 1))
                    
                    fitness_v /= np.sum(fitness_v)
                    
                    chosen_indices = np.random.choice(list(range(n_ind)), N, True, fitness_v)
                    
                    gametes = np.array(inds_in_deme)[chosen_indices].reshape(N*2)
                    
                    new_inds = gametes[np.concatenate((np.arange(1, 2*N), [0]), 0)].reshape(N,2)
                    
                    alleles[deme] = new_inds
        
            # Mutation
            
            n_mut = np.random.poisson(lamb)
            
            for i in range(n_mut):
                deme = np.random.randint(K)
                ind = np.random.randint(N)
                allele = np.random.randint(2)
                
                alleles[deme][ind][allele] = 1 - alleles[deme][ind][allele]
            
            # Iteration and calculation of frequencies
            
            if (not iteration % resolution) and (iteration >= delay):
            
                freqs = [0 for i in range(K)]
                
                for deme in range(K):
                    
                    freqs[deme] = np.sum(alleles[deme]) / (2*N)
            
                H = 1 - (sum(freqs)/K)**2 - (sum(1 - np.array(freqs))/K)**2

                H_v.append(H)
                
                w_pq = 4 * sum([freqs[i] * (1-freqs[i]) for i in range(K)])
                w_pq_v.append(w_pq)
                    
                freqs_v.append(freqs)
                
                if not equil:
                    
                    if len(H_v) >= len_H_v:
                        
                        # Criteria for equilibrium
                        #if abs(sum(H_v[-len_H_v:-len_H_v_by_2]) - \
                        #       sum(H_v[-len_H_v_by_2:]))/len_H_v_by_2 < tol_equil:
                        if True:    
                            equil = True
                            print('equil')
                        
                else:
                    
                    w_iter += 1
                    
                    if w_iter == n_w:
                        done = True
                        
            iteration += 1
            
        # Calculations each realization
        
        w_pq_mean = np.mean(w_pq_v[-n_w:])
        w_max_mean = 3/4 * w_pq_mean

        print('w_pq = %.5E' % w_pq_mean)
        print('w_max = %.5E' % w_max_mean)
        
        w_pq_vv.append(w_pq_v)
        w_pq_mean_v.append(w_pq_mean)
        freq_means = [sum([freqs_v[i][deme] for i in range(len(freqs_v))])/len(freqs_v) for deme in range(K)]
        
        freqs_time_var = [np.var([freqs_v[i][deme] for i in range(len(freqs_v))]) for deme in range(K)]
        F = sum(freqs_time_var)/sum([p*(1-p) for p in freq_means])
        F_v.append(F)         
             
        freq_means_v.append(freq_means)
        
        freqs_vv.append(freqs_v)
        
        print('%d iterations' % iteration)
                        
        #plt.cla()
        #plt.axis([0, K, 0, 1])
        #plt.plot(list(range(K)), freqs, 'r-*')
        #plt.xlabel('Deme')
        #plt.ylabel('f(A)')
        #plt.grid(True)
        #plt.title(r'f(A) for each deme in iteration %d. $\sigma _D$ = %.1E, s = %.1E' \
        #          % (iteration, sigma, s))
        #plt.savefig('0119_sigma_%.1E_s_%.1E.eps' % (sigma, s), format='eps', dpi=900)
        
        
        #w_max = 3 / 4 * w_2_mean
        #w_theo = sqrt(3*sigma**2 / s)
        
    # Final calculations and plots
        
    w_theo = m.sqrt(3*sigma**2 / s) # !!!!!!!!! 
    
    min_len = min([len(i) for i in w_pq_vv])
    
    F_mean = np.mean(F_v)
    
    av_cline = [[np.mean([freqs_vv[r_i][t][deme] for r_i in range(realizations)])\
                for deme in range(K)] for t in range(min_len)]

    w_pq_of_av_cline_v = [4 * sum([av_cline[t][deme] * (1-av_cline[t][deme]) for deme in range(K)]) \
                          for t in range(min_len)]
    w_max_of_av_cline_v = [3/4 * i for i in w_pq_of_av_cline_v]
    if plot:
        '''
        av_cline_conf_int = [[stats.bayes_mvs([freqs_vv[r_i][t][deme] for r_i in range(realizations)])[0][1] \
                    for deme in range(K)] for t in range(min_len)]
                    
        av_cline_low = np.array([[av_cline_conf_int[t][deme][0] for deme in range(K)]\
                                 for t in range(min_len)])
        col_mean = np.nanmean(av_cline_low, axis=0)
        indices = np.where(np.isnan(av_cline_low))
        av_cline_low[indices] = np.take(col_mean, indices[1])
        
        av_cline_high = np.array([[av_cline_conf_int[t][deme][1] for deme in range(K)]\
                                 for t in range(min_len)])
        col_mean = np.nanmean(av_cline_high, axis=0)
        indices = np.where(np.isnan(av_cline_high))
        av_cline_high[indices] = np.take(col_mean, indices[1])
        '''
        
        w_max_of_av_cline_v_i_v = []
        for i in range(bootstrap_n):
            indices = np.random.randint(0, realizations, realizations)
            av_cline_i = [[np.mean([freqs_vv[r_i][t][deme] for r_i in indices])\
                                   for deme in range(K)] for t in range(min_len)]
            w_max_of_av_cline_v_i = [3 * sum([av_cline_i[t][deme] * (1-av_cline_i[t][deme]) for deme in range(K)]) \
                                             for t in range(min_len)]
            w_max_of_av_cline_v_i_v.append(w_max_of_av_cline_v_i)
            
#        w_max_of_av_cline_v_sds = np.sqrt(np.var(w_max_of_av_cline_v_i_v, 0))
#        
#        w_max_of_av_cline_v_low = np.array(w_max_of_av_cline_v) - 1.96*w_max_of_av_cline_v_sds
#        w_max_of_av_cline_v_high = np.array(w_max_of_av_cline_v) + 1.96*w_max_of_av_cline_v_sds
        w_max_of_av_cline_v_low = np.percentile(w_max_of_av_cline_v_i_v, 2.5, 0)
        w_max_of_av_cline_v_high = np.percentile(w_max_of_av_cline_v_i_v, 97.5, 0)
#        w_max_of_av_cline_v_high = [3 * np.nansum([av_cline_high[t][deme] * (1-av_cline_high[t][deme]) for deme in range(K)]) \
#                          for t in range(min_len)]  
        

    
    av_w_pq_of_cline = sum(w_pq_mean_v)/realizations

    av_w_max_of_cline = av_w_pq_of_cline * 3/4
    
    av_w_pq_of_cline_v = [sum([w_pq_vv[r_i][t] for r_i in range(realizations)])/realizations \
                           for t in range(min_len)]
    percentile_low_w_max_v = 3/4 * np.percentile(w_pq_vv, 2.5, 0)
    percentile_high_w_max_v = 3/4 * np.percentile(w_pq_vv, 97.5, 0)
    
    
    av_w_max_of_cline_v = [3/4 * i for i in av_w_pq_of_cline_v]

    w_max_of_av_cline = sum(w_max_of_av_cline_v[-n_w:])/n_w

    w_pq_of_av_cline = w_max_of_av_cline * 4/3

    w_pq_of_av_cline_v_modif = [w_pq_of_av_cline_v[t] * (1 - F_mean) for t in range(min_len)]
    w_max_of_av_cline_v_modif = [w_max_of_av_cline_v[t] * (1 - F_mean) for t in range(min_len)]

    #w_max_of_av_cline_v = 
    x = list(range(min_len))
    
    labels = [r'$w(\overline{p})$', r'$\overline{w(p)}$', '95% CI',\
              '95% CI', r'$(1-\overline{F} ) w(\overline{p})$', r'$\sqrt{\frac{3}{s}} \sigma$']
              
    args = [[x,x,x,x,x,x], [w_max_of_av_cline_v, av_w_max_of_cline_v, percentile_low_w_max_v, \
            percentile_high_w_max_v, w_max_of_av_cline_v_modif, [w_theo for i in range(min_len)]],\
                'iteration / %d' % resolution, 'w', r'w(time), $\sigma =$ %.1f, s = %.3f' % (sigma, s),\
                '0119_w_s_%.3f_sigma_%.1f_rscl_%.2f.eps' % (s, sigma, rescale * rescale_factor), labels]
    
    if rescale:
        args[-3] = r'w(time), $\sigma =$ %.1f, s = %.3f, rescale = %.1f' % (sigma, s, rescale_factor)
        
    if plot:
        PlotAndSave(*args)
        
    if tm:
        for j, v in enumerate([w_max_of_av_cline_v, av_w_max_of_cline_v, percentile_low_w_max_v,\
                               percentile_high_w_max_v, w_max_of_av_cline_v_modif]):
            v_temp = []
            for i in range(min_len):
                mean = np.mean(v[:i+1])
                v_temp.append(mean)
            args[1][j] = v_temp
        args[-2] = '0119_w_s_%.3f_sigma_%.1f_rscl_%.2f_tm_%d.eps' % (s, sigma, rescale * rescale_factor, tm) 
        
        if plot:
            PlotAndSave(*args)
    
    
    with open(filename, 'a') as f:
        f.write('%f,%d,%d,%f,%f,%f,%f,%d,%f,%d\n' % \
                (s, K, N, sigma, av_w_pq_of_cline, \
                 w_pq_of_av_cline, w_theo, realizations, \
                 F_mean, n_w*resolution))
    
    if savedata:
        filenamedata = 'freqs_vv_s_%.1e_sigma_%.1e_rscl_%.2f_' % (s, sigma, rescale * rescale_factor)
        
        i = 1
        
        while os.path.isfile('D:\\' + filenamedata + str(i) + '.txt'):
            i += 1
        
        with open('D:\\' + filenamedata + str(i) + '.txt', 'w') as f:
            json.dump(freqs_vv, f)
            
        
        
        
        
                
        
                        
                        
                        
                        
                        
                        
                
                
                
                
                
                
                
                
                
                




