# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 13:55:04 2017

@author: oskarfridell
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt

def PlotAndSave(xarrays,yarrays,xlabel,ylabel,title,savename, labels):
   
    for i in range(len(yarrays)):
        plt.plot(xarrays[i], yarrays[i][0:len(xarrays[i])],'o',\
                 c = plt.cm.gnuplot(i/len(yarrays)), label = labels[i])
        
        #legend_line.append(mlines.Line2D([], [], color=plt.cm.rainbow(i/len(yarrays)), label=r'$S_%d$' % i))
    
    # plt.legend(handles=legend_line)
    plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.subplots_adjust(left=0.1, right=.7, top=0.9, bottom=0.1)
    plt.savefig(savename, format='png', dpi=900)
    plt.show()
    plt.clf()
 

w_pq_mean_mean_vv = [[0.64306952182952182, 1.1808741978021979, 1.8779971190476188, 2.7661961846153851, 3.8710425502489398, 5.3683226834381559, 7.6427288752142086, 11.397883451851852], [2.161908752598753, 3.2800129523809516, 4.6947305357142852, 6.4264942666666673, 8.690016383920339, 11.848306004192874, 16.092516677052501, 21.808327044444447], [3.0722429729729726, 4.5524918388278381, 6.385186428571429, 8.7582634974358982, 11.862566992439612, 16.161324276729562, 21.612953777847014, 29.490592088888892], [7.6392538461538466, 11.059353728937733, 15.361452964285714, 20.791433415384621, 28.192859708648349, 38.04132685115303, 50.874360250817858, 68.275453266666645]]
var_vv = [[0.0014158396438742928, 0.0029187292671989697, 0.0066523027445700026, 0.019489929588834297, 0.06323109549957176, 0.25466016004800329, 1.1280919094847213, 5.2833776447877456], [0.0030227010987024612, 0.0056836756060087224, 0.013936669020663241, 0.024864006170540753, 0.12306773748124372, 0.47336146272884527, 0.88335946398413756, 6.0147578234834267], [0.00442741884405799, 0.0056615984598128687, 0.015940098537967771, 0.034539459542496168, 0.14185257915407579, 0.3989218050010353, 2.1800874335729525, 5.6438875318465147], [0.0069758578155437109, 0.01197662025746493, 0.045833537943540857, 0.11139604964782845, 0.34862677108798318, 0.88322370308029119, 2.5001080710808834, 12.050666790427812]]
var_vv_theo = [[0.09924465050361478, 0.2353461546318528, 0.5580936828225578, 1.3234486847412896, 3.1383914117879956, 7.442298871988413, 17.648471854708717, 41.85112210136604], [0.15879144080578364, 0.37655384741096454, 0.8929498925160926, 2.1175178955860634, 5.021426258860793, 11.907678195181461, 28.23755496753395, 66.96179536218567], [0.19848930100722956, 0.4706923092637056, 1.1161873656451156, 2.6468973694825793, 6.276782823575991, 14.884597743976826, 35.296943709417434, 83.70224420273207], [0.3969786020144591, 0.9413846185274112, 2.232374731290231, 5.2937947389651585, 12.553565647151983, 29.769195487953652, 70.59388741883487, 167.40448840546415]]
F_mean_vv = [[0.02132431651224876, 0.031120966227700279, 0.044243776359073214, 0.06068353458107005, 0.079732519108155625, 0.10433420085727055, 0.12643606294999074, 0.15654397808450193], [0.018241612068083304, 0.022696673509434877, 0.027228814038036439, 0.033609284388282673, 0.042232860455555098, 0.052683359259064706, 0.06403070725783129, 0.08111818527936536], [0.017073253557231514, 0.019386398521341403, 0.023058781873327047, 0.027655722914572821, 0.033697819325995075, 0.04189891076366057, 0.051094542254717046, 0.063967615930828906], [0.013674837972889703, 0.01455546242141679, 0.01603600137058531, 0.017759192296851925, 0.020468018130696278, 0.023983307667550715, 0.028854397289142634, 0.036931400050347625]]

sigma_v = [.5, .8, 1, 2]
s_vv = [[10**(-i/4) for i in range(3,11)] for j in range(4)]
sigma_vv = [[i for j in range(3,11)] for i in sigma_v]

#var_vv = np.array(var_vv).reshape(1,32)
#var_vv = var_vv.squeeze()
#var_vv_theo = np.array(var_vv_theo).reshape(1,32)
#var_vv_theo = var_vv_theo.squeeze()

neigh_vv = [[2*sqrt(pi)*sigma[0]*100 for i in sigma] for sigma in sigma_vv]

x = [[1/(sqrt(s_vv[j][i]) * neigh_vv[j][i]) for i in range(8)] for j in range(4)]
#x = np.array(x).reshape(1,32)
#x = x.squeeze()

for i, sigma in enumerate(sigma_v):
    
    plt.plot(np.log2(x[i]), np.log10(var_vv[i]), c = plt.cm.rainbow(i/len(sigma_v)), label = r'$\sigma = %.1f$' % sigma)
    plt.plot(np.log2(x[i]), np.log10(var_vv_theo[i]), c = plt.cm.rainbow((i+.5)/len(sigma_v)), label = 'Theory')
    
    
    
#plt.plot(w_pq_mean_mean_vv[i], var_vv_theo[i], c = plt.cm.rainbow((i+.5)/len(sigma_v)), label = 'var theo')
    
    #legend_line.append(mlines.Line2D([], [], color=plt.cm.rainbow(i/len(yarrays)), label=r'$S_%d$' % i))

# plt.legend(handles=legend_line)
plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
plt.xlabel(r'$\log_2 (1/\mathcal{N} \sqrt{s})$')
plt.ylabel('$\log_{10}$ (var(w))')
plt.title(r'w as a function of $1/\mathcal{N} \sqrt{s}$')
plt.grid(True)
plt.subplots_adjust(left=0.1, right=.75, top=0.9, bottom=0.1)
plt.savefig('0202_var_of_w_neigh_all_time.png', format='png', dpi=900)
plt.show()
plt.clf()

'''

for i, sigma in enumerate(sigma_v):
    
    plt.plot(w_pq_mean_mean_vv[i], F_mean_vv[i], c = plt.cm.rainbow(i/len(sigma_v)), label = r'$\sigma = %.1f$' % sigma)
        
        #legend_line.append(mlines.Line2D([], [], color=plt.cm.rainbow(i/len(yarrays)), label=r'$S_%d$' % i))

# plt.legend(handles=legend_line)
plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
plt.xlabel('w')
plt.ylabel('<F>')
plt.title('<F> for different w')
plt.grid(True)
plt.subplots_adjust(left=0.1, right=.75, top=0.9, bottom=0.1)
plt.savefig('0202_F.png', format='png', dpi=900)
plt.show()
plt.clf()



'''














