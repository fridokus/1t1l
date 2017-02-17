# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 11:34:06 2017

@author: Issun
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt
import datetime


def PlotAndSave(xarrays,yarrays,xlabel,ylabel,title,savename):
    
    for i in range(len(yarrays)):
        if i ==0:
            plt.plot(xarrays[i], yarrays[i], Colors(i))
        if i ==1:
            plt.plot(xarrays[i], yarrays[i], Colors(i)+'--')
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.savefig(savename, format='eps', dpi=900)
    plt.show()
    plt.clf()


# Parameters

K = 150
N = 100

s_v = [5e-2, 1e-2, 5e-3, 1e-3, 5e-4]
sigma_v = [.5, 1, 1.5, 2, 3, 5]

#s = 1e-3
#sigma = 1.5
f
mu = 1e-5

tol_equil = 1e-3

env_change = int(K/2)

lamb = K * N * 2 * mu

w_calc_distance = 1
sum_distance = 10
n_w = 100

#n_traits = len(env_change)

f = open('0118_output.txt', 'a')
f.write('\n')
f.write(str(datetime.datetime.now()))
f.write('\n')
f.write('s, sigma, w_mle, w_pq \n')
f.close()


for s in s_v:
    print('s = %f' % s)
    
    for sigma in sigma_v:
        print('sigma = %f' % sigma)
        
        # Initialization
        
        alleles = [[[0, 0] for i in range(N)] if j < env_change else [[1, 1] for i in range(N)] for j in range(K)]
        
        done = equil = False
        
        iteration = 0
        
        freq_means_r = []
        freq_means_l = []
        
        last_mean_recent_freq_r = 0
        last_mean_recent_freq_l = 0
        
        w_v = []
        w_2_v = []
        
        w_iter = 0 
        
        while not done:
            
            # Migration
            
            alleles_new = [[] for j in range(K)]
            
            for deme in range(K):
                
                for ind in alleles[deme]:
                    
                    move = int(np.random.normal(0, sigma, 1))
                    
                    if deme + move < 0:
                        new_loc = 0
                    elif deme + move >= K:
                        new_loc = K - 1
                    else:
                        new_loc = deme + move
                    
                    alleles_new[new_loc].append(ind)
        
            # Mating and selection
            
            for deme in range(K):
                
                inds_in_deme = alleles_new[deme]
        
                n_ind = len(inds_in_deme)
        
                if deme < env_change: # <
                    fitness_v = (1-s)**np.sum(np.array(inds_in_deme), 1)
                else:
                    fitness_v = (1-s)**(2-np.sum(np.array(inds_in_deme), 1))
        
                sum_fitness = sum(fitness_v)
                fitness_v = fitness_v/sum_fitness
                
                chosen_indices = np.random.choice(list(range(n_ind)), N, True, fitness_v)
                chosen_inds = [inds_in_deme[i] for i in chosen_indices]
            
                offspring = [[chosen_inds[i][0], chosen_inds[i - 1][1]] for i in range(N)]
                             
                alleles[deme] = [i for i in offspring]
        
            # Mutation
            
            n_mut = np.random.poisson(lamb)
            
            for i in range(n_mut):
                deme = np.random.randint(K)
                ind = np.random.randint(N)
                allele = np.random.randint(2)
                
                alleles[deme][ind][allele] = 1 - alleles[deme][ind][allele]
            
            # Iteration and calc freqs
            
            iteration += 1
            
            '''
            if iteration > 500:
                done = True
            '''
            
            if not iteration % 20:
            
                freqs = [0 for i in range(K)]
                
                for deme in range(K):
                    
                    freqs[deme] = np.sum(alleles[deme]) / (2*N)
                    
                freqs_r = [freqs[i] for i in range(env_change, K)]
                freqs_l = [freqs[i] for i in range(env_change)]
                    
                freqs_mean_r = sum(freqs_r)/env_change
                freqs_mean_l = sum(freqs_l)/env_change
        
                freq_means_r.append(freqs_mean_r)
                freq_means_l.append(freqs_mean_l)
                
                index_list = list(range(int(max(0, iteration/20 - sum_distance)), int(iteration/20)))
                n_indices = len(index_list)
                
                mean_recent_freq_r = sum([freq_means_r[i] for i in index_list])/n_indices
                mean_recent_freq_l = sum([freq_means_l[i] for i in index_list])/n_indices
            
                # H = 1 - (sum(freqs)/K)**2 - (sum(1 - np.array(freqs))/K)**2
                
                if not equil:
                    if abs(mean_recent_freq_l - last_mean_recent_freq_l) < tol_equil and \
                         abs(mean_recent_freq_r - last_mean_recent_freq_r) < tol_equil:
                        equil = True
                        print('equil')
                    
                    last_mean_recent_freq_r = mean_recent_freq_r
                    last_mean_recent_freq_l = mean_recent_freq_l
                        
                else:
                    
                    w = w_calc_distance/(freqs[int(env_change + (w_calc_distance-1)/2)] -\
                                     freqs[int(env_change - (w_calc_distance+1)/2)])
                    
                    w_iter += 1
                    w_v.append(w)
                    
                    w_2 = 4 * sum([freqs[i] * (1-freqs[i]) for i in range(K)])
                    w_2_v.append(w_2)
                    
                    if w_iter == n_w:
                        
                        done = True
                        
                        while w_v.count(inf):
                            w_v.remove(inf)
                        
                        w_mean = sum(w_v)/n_w
                        w_2_mean = sum(w_2_v)/n_w
        
                        print('w_mle = %.5E' % w_mean)
                        print('w_pq = %.5E' % w_2_mean)
                        
        plt.cla()
        plt.axis([0, K, 0, 1])
        plt.plot(list(range(K)), freqs, 'r-*')
        plt.xlabel('Deme')
        plt.ylabel('f(A)')
        plt.grid(True)
        plt.title(r'f(A) for each deme in iteration %d. $\sigma _D$ = %.1E, s = %.1E' \
                  % (iteration, sigma, s))
        plt.savefig('0118_sigma_%.1E_s_%.1E.eps' % (sigma, s), format='eps', dpi=900)
        
        f = open('0118_output.txt', 'a')
        f.write('%f, %f, %f, %f \n' % (s, sigma, w_mean, w_2_mean))
        f.close()
        
        
        
        
        
                
        
                        
                        
                        
                        
                        
                        
                
                
                
                
                
                
                
                
                
                




