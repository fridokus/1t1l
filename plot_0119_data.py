# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 08:46:31 2017

@author: rax
"""
from math import *
import numpy as np
import matplotlib.pyplot as plt
import datetime
import os
import csv
import pandas as pd

def w_theo(s, sigma):
    return sqrt(3/s) * sigma

def PlotAndSave(xarrays,yarrays,xlabel,ylabel,title,savename,labels):
   
    for i in range(len(yarrays)):
        if i == 3:
            marker = '-'
        else:
            marker = '.'
        
        plt.plot(xarrays[i], yarrays[i][0:len(xarrays[i])],marker,\
                 c = plt.cm.gnuplot(i/len(yarrays)), label = labels[i])
        
        #legend_line.append(mlines.Line2D([], [], color=plt.cm.rainbow(i/len(yarrays)), label=r'$S_%d$' % i))
    
    # plt.legend(handles=legend_line)
    plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.subplots_adjust(left=0.12, right=.7, top=0.9, bottom=0.1)
    plt.savefig(savename, format='eps', dpi=900)
    plt.show()
    plt.clf()

def PlotAndSave2(xarrays,yarrays,xlabel,ylabel,title,savename):
   
    for i in range(len(yarrays)):
        if i == 1:
            marker = '-'
        else:
            marker = '.'
        
        plt.plot(xarrays[i], yarrays[i][0:len(xarrays[i])],marker,\
                 c = plt.cm.gnuplot(i/len(yarrays)))
        
        #legend_line.append(mlines.Line2D([], [], color=plt.cm.rainbow(i/len(yarrays)), label=r'$S_%d$' % i))
    
    # plt.legend(handles=legend_line)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    axes = plt.gca()
    axes.set_ylim([0, .4])
    plt.title(title)
    plt.grid(True)
    plt.subplots_adjust(left=0.1, right=.9, top=0.9, bottom=0.1)
    plt.savefig(savename, format='eps', dpi=900)
    plt.show()
    plt.clf()
            
df = pd.read_csv('0119_data_long.txt')

sigma_vv = [[.5, 1]]
K_v = [150]


for K_i, sigma_v in enumerate(sigma_vv):

    s_vv = []
    av_w_of_p_vv = []
    w_of_av_p_vv = []
    N_v = []
    #w_theo_vv = []
    F_vv = []
    w_of_av_p_modif_vv = []
    rel_diff_vv = []

    K = K_v[K_i]
    
    for i, sigma in enumerate(sigma_v):
        s_vv.append(df['s'][df['sigma'] == sigma][df['K'] == K])
        av_w_of_p_vv.append(df['<w(p)>'][df['sigma'] == sigma][df['K'] == K]*3/4)
        w_of_av_p_vv.append(df['w(<p>)'][df['sigma'] == sigma][df['K'] == K]*3/4)
        N_v.append(df['N'][df['sigma'] == sigma][df['K'] == K])
        #w_theo_vv.append(df['w_theo'][df['sigma'] == sigma][df['K'] == K])
        w_theo_here = df['w_theo'][df['sigma'] == sigma][df['K'] == K]
        F_vv.append(df['<F>'][df['sigma'] == sigma][df['K'] == K])
        w_of_av_p_modif_vv.append((1 - F_vv[i]) * w_of_av_p_vv[i])
        rel_diff_vv.append(np.abs(w_of_av_p_vv[-1] - w_theo_here)/w_theo_here)
    
    
    s_v_theo = np.arange(np.min(s_vv[0])/4*3.5, np.max(s_vv[0]) + np.min(s_vv[0]), 1e-4)
    w_theo_vv = [[w_theo(s, sigma) for s in s_v_theo] for sigma in sigma_v]
        
    labels = [r'$\overline{w(p)}$', r'$w(\overline{p})$', r'$(1-\overline{F} ) w(\overline{p})$', r'$\sqrt{\frac{3}{s}} \sigma$']
              
    for i, sigma in enumerate(sigma_v):
        N = N_v[i]
        N = N[N.index[0]]
    
        args = [[1/np.sqrt(s_vv[i]) for j in range(3)] + [1/np.sqrt(s_v_theo)],\
                [av_w_of_p_vv[i], w_of_av_p_vv[i], w_of_av_p_modif_vv[i], w_theo_vv[i]],\
                    r'$1/\sqrt{s}$', r'$w$', r'$w(1/\sqrt{s})$, $\sigma$ = %.2f, K = %d, N = %d' % (sigma, K, N),\
                    'w_of_s_K_%d_nr_%d_0.eps' % (K, i), labels]
                
        PlotAndSave(*args)

#sigma_v = np.squeeze(sigma_vv)    
    labels = [r'$\sigma$ = ' + str(sigma_v[i]) for i in range(len(sigma_v))]
    args = [[1/np.sqrt(s_vv[i]) for i in range(len(sigma_v))],\
                [rel_diff_vv[i] for i in range(len(sigma_v))],\
                    r'$1/\sqrt{s}$', r'$\left| w(\overline{p}) - \sqrt{\frac{3}{s}} \sigma \right| /\sqrt{\frac{3}{s}} \sigma$',\
                    r'Simulated and theoretical $w$, relative difference',\
                    'rel_diff_K_%d_0.eps' % K, labels]
                
    PlotAndSave(*args)
        
'''
s = round(float(10**(-7/4)), 6)
sigma_vv = []

av_w_of_p_vv = []
w_of_av_p_vv = []
w_theo_vv = []
F_vv = []
w_of_av_p_modif_vv = []

K_v = [150, 50]

for i, K in enumerate(K_v):
    sigma_vv.append(df['sigma'][df['s'] == s][df['K'] == K])
    #sigma_vv[-1].sort_values(inplace=True)
    av_w_of_p_vv.append(df['<w(p)>'][df['s'] == s][df['K'] == K]*4/3)
    w_of_av_p_vv.append(df['w(<p>)'][df['s'] == s][df['K'] == K]*4/3)
    w_theo_vv.append(df['w_theo'][df['s'] == s][df['K'] == K])
    F_vv.append(df['<F>'][df['s'] == s][df['K'] == K])
    w_of_av_p_modif_vv.append((1 - F_vv[i]) * w_of_av_p_vv[i])
	
	
for i, K in enumerate(K_v):

    args = [[np.log(sigma_vv[i]) for j in range(4)], [av_w_of_p_vv[i], w_of_av_p_vv[i], w_theo_vv[i], w_of_av_p_modif_vv[i]],\
                'ln($\sigma _D$)', 'w', r'w($\sigma _D$), s = %.3f' % s,\
                'w_of_sigma_%d.png' % i, labels]
            
    if K != 150:
        rescale_factor = int(150/K)
        args[-3] = r'w($\sigma$), s = %.3f, rescale factor %d' % (s, rescale_factor)

    PlotAndSave(*args)



'''








             
             
             
             
             

             
             
            
            
