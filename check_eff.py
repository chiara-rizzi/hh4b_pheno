import argparse
import os
import numpy as np
import pandas as pd
import lorentz 
import math

parser=argparse.ArgumentParser()
parser.add_argument("-f","--file",help="file name", type=str, default="ggHH_smallR-jets.lhco")
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
# #  typ      eta      phi      pt    jmas   ntrk   btag  had/em   dum1   dum2

n_events=0
n_bjets=0 #count the number of b-tagged jets, to be set to zero at each event
n_events_4b=0 # number of events with 4b
n_events_dR=0 # number of events with 4b
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
             if diff_appo < diff:
                  #if True:
                  if h1_appo.pt > h2_appo.pt:
                       dR_h1_appo = dR(bjets[i_h1_j1], bjets[i_h1_j2_appo])
                       dR_h2_appo = dR(bjets[i_h2_appo[0]], bjets[i_h2_appo[1]])
                  else:
                       dR_h2_appo = dR(bjets[i_h1_j1], bjets[i_h1_j2_appo])
                       dR_h1_appo = dR(bjets[i_h2_appo[0]], bjets[i_h2_appo[1]])                       
                  m4j_appo = (h1_appo+h2_appo).m()
                  if  pass_dR(dR_h1_appo, dR_h2_appo, m4j_appo):
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

# loop on all the particles
for index, row in data.iterrows():
    #if n_events>500: break
     if row['#']==0:
          if n_events%5000==0:
               print("Looking at event: "+str(n_events))
          # just finished previous event: check if it had >= 4 b-tagged jets, and compute quantities 
          if n_bjets ==4:
               n_events_4b+=1
               i_h1, i_h2 = pair_jets(bjets_event) 
               h1 = bjets_event[i_h1[0]] + bjets_event[i_h1[1]]
               h2 = bjets_event[i_h2[0]] + bjets_event[i_h2[1]]
               m4j = (h1+h2).m()
               dR_h1 =  dR(bjets_event[i_h1[0]], bjets_event[i_h1[1]])
               dR_h2 =  dR(bjets_event[i_h2[0]], bjets_event[i_h2[1]])
               if pass_dR(dR_h1, dR_h2, m4j):
                    n_events_dR +=1
          # it's a new event! set event-by-event counters to zero 
          n_bjets=0 
          bjets_event = []        
          # increase the number of events        
          n_events+=1 
     if row["typ"]==4 and row["btag"]>0 and row["pt"]>40:
          n_bjets+=1
          n = set_lv(row["pt"], row["eta"], row["phi"], row["jmas"])
          bjets_event.append(lorentz.FourMomentum(n[0],n[1],n[2], n[3]))
          
# For the way my code is setup, the counter n_events_4b is increased on the following event. 
# This means that for the last event we need to make the check (and in case increase the counter) outside the loop
if n_bjets >=4:
     n_events_4b+=1
        

print("N events: "+str(n_events))
print("N events >=4 b-jets pt>40 GeV: "+str(n_events_4b))
print("N events dR: "+str(n_events_dR))
print("Acc x Eff 4b:", 1.0*n_events_4b/n_events * 1/(0.58*0.58))
print("Acc x Eff 4b:", 1.0*n_events_dR/n_events * 1/(0.58*0.58))


