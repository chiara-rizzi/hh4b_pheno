import argparse
import os
import numpy as np
import pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument("-f","--file",help="file name", type=str, default="ggHH_smallR-jets.lhco")
args=parser.parse_args()

print("Ciao! :) ")

data = pd.read_csv(args.file, sep="\s+")
#print(data)
# #  typ      eta      phi      pt    jmas   ntrk   btag  had/em   dum1   dum2

n_events=0
n_bjets=0 #count the number of b-tagged jets, to be set to zero at each event
n_events_4b=0 # number of events with 4b

# loop on all the particles
for index, row in data.iterrows():
    if row['#']==0:
        if n_events%5000==0:
            print("Looking at event: "+str(n_events))
        if n_bjets >=4:
            n_events_4b+=1 # just finished previous event: check if it had >= 4 b-tagged jets
        n_bjets=0 # it's a new event! set event-by-event counters to zero 
        n_events+=1 # increase the number of events
    if row["typ"]==4 and row["btag"]>0 and row["pt"]>40:
        n_bjets+=1

# For the way my code is setup, the counter n_events_4b is increased on the following event. 
# This means that for the last event we need to make the check (and in case increase the counter) outside the loop
if n_bjets >=4:
    n_events_4b+=1
        

print("N events: "+str(n_events))
print("N events >=4 b-jets pt>40 GeV: "+str(n_events_4b))




