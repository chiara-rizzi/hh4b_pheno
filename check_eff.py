import argparse
import os
import numpy as np
import pandas as pd
import lorentz 
import math

parser=argparse.ArgumentParser()
parser.add_argument("-f","--file",help="file name", type=str, default="ggHH_smallR-jets.lhco")
parser.add_argument("-4bex","--exactly4b", help="select only events with exactly 4b (as opposed to >= 4b)", type=int, default=0)
args=parser.parse_args()

print("Ciao! :) ")

def set_lv(pt, eta, phi, m):
     x = pt*math.cos(phi)
     y = pt*math.sin(phi)
     z = pt*math.sinh(eta)
     t = math.sqrt(x*x + y*y + z*z + m*m)
     return t,x,y,z

data = pd.read_csv(args.file, sep="\s+")

#print(data)
# structude of the input file: 
# #  typ      eta      phi      pt    jmas   ntrk   btag  had/em   dum1   dum2

n_events = {}
labels = ['tot','4b', 'dR', 'pT', 'eta', 'Xhh', 'XWt']
for l in labels:
     n_events[l] = 0
eff={}
eff['atlas']={}
eff['us']={}
eff['atlas']['tot']=1
eff['atlas']['4b']=0.49
eff['atlas']['dR']=0.445
eff['atlas']['pT']=0.42
eff['atlas']['eta']=0.37
eff['atlas']['Xhh']=0.19
eff['atlas']['XWt']=0.175



#n_events=0
n_bjets=0 # count the number of b-tagged jets, to be set to zero at each event
#n_events_4b=0 # number of events with 4b
#n_events_dR=0 # number of events that pass dR selection
#n_events_pT=0 # number of events that pass pT selection
#n_events_eta=0 # number of events that pass eta selection
bjets_event = [] # list of lorentz vectors of the b-jets in the event

def pair_jets(bjets):
    #print("bjets[0].pt",bjets[0].pt)
    i_h1 = (0,0)
    i_h2 = (0,0)
    if len(bjets) < 4: 
        return i_h1, i_h2
    i_h1_j1 = 0
    diff = 999999
    for i_h1_j2_appo in range(1,len(bjets)):
        h1_appo = bjets[i_h1_j1] + bjets[i_h1_j2_appo]
        i_h2_list = [x for x in range(1,len(bjets)) if x != i_h1_j2_appo]
        i_h2_j1_appo = i_h2_list[0]
        for i_h2_j2_appo in i_h2_list[1:]:             
             i_h2_appo = [i_h2_j1_appo, i_h2_j2_appo]
             h2_appo = bjets[i_h2_appo[0]] + bjets[i_h2_appo[1]]
             diff_appo = math.fabs(h1_appo.m() - h2_appo.m() )
             # the dR selection makes a difference between dR1 and dR2
             # so I have to treat them differently
             # pretty ugly in the code since then I do it again later, will clean it up
             if h1_appo.pt > h2_appo.pt:
                  dR_h1_appo = dR(bjets[i_h1_j1], bjets[i_h1_j2_appo])
                  dR_h2_appo = dR(bjets[i_h2_appo[0]], bjets[i_h2_appo[1]])
             else:
                  dR_h2_appo = dR(bjets[i_h1_j1], bjets[i_h1_j2_appo])
                  dR_h1_appo = dR(bjets[i_h2_appo[0]], bjets[i_h2_appo[1]])                       
             m4j_appo = (h1_appo+h2_appo).m()
             # consider only the pairings that satisfy the dR selecion
             if  pass_dR(dR_h1_appo, dR_h2_appo, m4j_appo):
                  # choose the one with the smallest mass difference
                  if diff_appo < diff:
                       diff = diff_appo
                       i_h1 = 0, i_h1_j2_appo
                       i_h2 = i_h2_appo[0], i_h2_appo[1]
    h1 = bjets[i_h1[0]] + bjets[i_h1[1]]
    h2 = bjets[i_h2[0]] + bjets[i_h2[1]]
    #if diff > 999998:
    #     print("no good match")
    #     no_good+=1
    if h1.pt > h2.pt:
         return i_h1, i_h2
    else:
         return i_h2, i_h1

def dR(j1, j2):
    deta = j1.eta - j2.eta
    dphi = j1.phi - j2.phi
    return math.sqrt(deta*deta+dphi*dphi)

def pass_dR(dR_h1, dR_h2, m4j, do_print=False):
     if do_print:
          print("dR_h1",dR_h1,"   dR_h2",dR_h2)
          print("m4j",m4j)
     if m4j < 1250:
          if (360./m4j) - 0.5 < dR_h1 and dR_h1 < (653./m4j) + 0.475:
               if (235./m4j) < dR_h2 and dR_h2 < (875./m4j) + 0.35: 
                    return True
     else: 
          if 0 <  dR_h1 and dR_h1 < 1:
               if 0 <  dR_h2 and dR_h2 < 1:
                    return True
     return False


def pass_pT(h1, h2):
     m4j = (h1+h2).m()
     if h1.pt > 0.5*m4j - 103:
          if h2.pt > 0.33*m4j - 73:
               return True
     return False

def pass_dEta(h1, h2):
     if math.fabs(h1.eta - h2.eta) < 1.5:
          return True
     return False

# loop on all the particles
for index, row in data.iterrows():
    #if n_events>500: break
     if row['#']==0:
          if n_events['tot']%5000==0:
               print("Looking at event: "+str(n_events['tot']))
          # just finished previous event: check if it had >= 4 b-tagged jets, and compute quantities 
          if (not args.exactly4b and n_bjets >=4) or n_bjets==4:
               n_events['4b']+=1
               i_h1, i_h2 = pair_jets(bjets_event) 
               h1 = bjets_event[i_h1[0]] + bjets_event[i_h1[1]]
               h2 = bjets_event[i_h2[0]] + bjets_event[i_h2[1]]
               m4j = (h1+h2).m()
               dR_h1 =  dR(bjets_event[i_h1[0]], bjets_event[i_h1[1]])
               dR_h2 =  dR(bjets_event[i_h2[0]], bjets_event[i_h2[1]])
               if pass_dR(dR_h1, dR_h2, m4j):
                    n_events['dR'] +=1
                    if pass_pT(h1, h2):
                         n_events['pT'] +=1
                         if pass_dEta(h1,h2):
                              n_events['eta'] +=1
          # it's a new event! set event-by-event counters to zero 
          n_bjets=0 
          bjets_event = []        
          # increase the number of events        
          n_events['tot']+=1 
     if row["typ"]==4 and row["btag"]>0 and row["pt"]>40 and row["eta"]<2.5:
          n_bjets+=1
          n = set_lv(row["pt"], row["eta"], row["phi"], row["jmas"])
          bjets_event.append(lorentz.FourMomentum(n[0],n[1],n[2], n[3]))
          
# For the way my code is setup, the counter n_events_4b is increased on the following event. 
# This means that for the last event we need to make the check (and in case increase the counter) outside the loop
# Will put this back later on
#if n_bjets >=4:
#     n_events_4b+=1

# for the way the code is set up, I'm disregardin the last event        
n_events['tot'] = n_events['tot']-1.0


print('\nNumber of events')
for l in labels:
     print("N events",l,":", n_events[l])

print('\nAcceptance x Efficiency')
for l in labels:
     if l=='tot': 
          eff['us'][l]=1
     else:
          eff['us'][l]=1.0*n_events[l]/n_events['tot']* 1/(0.58*0.58)
          print(l,"  atlas:",eff['atlas'][l],'   us:',eff['us'][l])

print('\nRelative to previous selection')
for i in range(1,len(labels)):
     l = labels[i]
     l_prev=labels[i-1]
     if eff['us'][l_prev] >0:          
          print(l,"  atlas:",eff['atlas'][l]/eff['atlas'][l_prev],'   us:',eff['us'][l]/eff['us'][l_prev])
