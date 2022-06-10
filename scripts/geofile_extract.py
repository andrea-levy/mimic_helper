import numpy as np
import pandas as pd
import sys
import os.path


#help function accessible via python geofile_extract.py --help (or -h) 
import argparse
parser = argparse.ArgumentParser(description='Extract GEOMETRY file from TRAJECTORY for selected step.')
parser.add_argument('TRAJECTORY_PATH', type=str,
                    help='path to TRAJECTORY file')
parser.add_argument('ENERGIES_PATH', type=str, 
                    help='path to ENERGIES file')
parser.add_argument('STEP', type=int,
                    help='step at which configuration is extracted')
args = parser.parse_args()

#input files check
TRAJECTORY = args.TRAJECTORY_PATH
print(TRAJECTORY)
if not (os.path.isfile(TRAJECTORY)):
    print('\nERROR: Path to TRAJECTORY file does not exist.')
    print('Please check path: '+TRAJECTORY)
    exit()

ENERGIES = args.ENERGIES_PATH
if not (os.path.isfile(ENERGIES)):
    print('\nERROR: Path to ENERGIES file does not exist.')
    print('Please check input path: '+ENERGIES)
    exit()

step = args.STEP

#read energy file to extract temperature
#and check required step is present
EneFile = open(ENERGIES, 'r')
temperature = []
for line in EneFile:
    if(int(line.split()[0])==step):
        temperature = float(line.split()[2])
        break
if(temperature == []):
    print('\nERROR: step required ('+str(step)+') is not in TRAJECTORY file')
    print('Please check steps saved in TRAJECTORY file '+TRAJECTORY)
    exit()
EneFile.close()

#read trajectory to extract geometry at selected step
df = pd.DataFrame(columns=('step', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
dfi = 0

TrajFile = open(TRAJECTORY, 'r')
for line in TrajFile:
    if(int(line.split()[0])==step):
        df.loc[dfi] = line.split()
        dfi += 1
TrajFile.close()

#write output file
if(os.path.isabs(TRAJECTORY)):
    OUTPUT = '/'.join(TRAJECTORY.split('/')[:-1])+'/GEO_'+str(int(temperature))+'K'
else:
    OUTPUT = '.'+'/'.join(TRAJECTORY.split('/')[:-1])+'/GEO_'+str(int(temperature))+'K'

OutFile = open(OUTPUT, 'w')
for i in range(dfi):
    OutFile.write(str(df['x'].loc[i])+'\t'+str(df['y'].loc[i])+'\t'+str(df['z'].loc[i])+'\t'+str(df['vx'].loc[i])+'\t'+str(df['vy'].loc[i])+'\t'+str(df['vz'].loc[i])+'\n') 
OutFile.close()

#print output message
print('------')
print('Geometry file extracted from step '+str(step)+' saved as:')
print('\t'+OUTPUT)
print('Corresponding to temperature of '+str(temperature)+'K')
