import numpy as np
import pandas as pd
import sys
import os.path


#help function accessible via python geofile_extract.py --help (or -h) 
import argparse
parser = argparse.ArgumentParser(description='Compute temperature from velocities in TRAJECTORY file.')
parser.add_argument('TRAJECTORY_PATH', type=str,
                    help='path to TRAJECTORY file')
parser.add_argument('GEOMETRY_PATH', type=str, 
                    help='path to GEOMETRY.xyz file')
parser.add_argument('CPMDOUT_PATH', type=str, 
                    help='path to CPMD output file')
parser.add_argument('-p', '--plot', action='store_true',
                    help='plot temperatures')
args = parser.parse_args()

# function to get mass from atom name
def get_mass(element):
    elements = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np']
    masses = [1.00797,4.0026,6.941,9.01218,10.81,12.011,14.0067,15.9994,18.998403,20.179,22.98977,24.305,26.98154,28.0855,30.97376,32.06,35.453,39.948,39.0983,40.08,44.9559,47.9,50.9415,51.996,54.938,55.847,58.9332,58.7,63.546,65.38,69.72,72.59,74.9216,78.96,79.904,83.8,85.4678,87.62,88.9059,91.22,92.9064,95.94,-98,101.07,102.9055,106.4,107.868,112.41,114.82,118.69,121.75,127.6,126.9045,131.3,132.9054,137.33,138.9055,140.12,140.9077,144.24,-145,150.4,151.96,157.25,158.9254,162.5,164.9304,167.26,168.9342,173.04,174.967,178.49,180.9479,183.85,186.207,190.2,192.22,195.09,196.9665,200.59,204.37,207.2,208.9804,-209,-210,-222,-223,226.0254,227.0278,232.0381,231.0359,238.029,237.0482] #amu - source: https://www.angelo.edu/faculty/kboudrea/periodic/structure_numbers.html    
    return masses[elements.index(element)]

#input files check
TRAJECTORY = args.TRAJECTORY_PATH
print(TRAJECTORY)
if not (os.path.isfile(TRAJECTORY)):
    print('\nERROR: Path to TRAJECTORY file does not exist.')
    print('Please check path: '+TRAJECTORY)
    exit()

GEOMETRY = args.GEOMETRY_PATH
if not (os.path.isfile(GEOMETRY)):
    print('\nERROR: Path to GEOMETRY.xyz file does not exist.')
    print('Please check input path: '+GEOMETRY)
    exit()

CPMDOUT = args.CPMDOUT_PATH
if not (os.path.isfile(CPMDOUT)):
    print('\nERROR: Path to CPMD output file does not exist.')
    print('Please check input path: '+CPMDOUT)
    exit()

#extract number of QM atoms and total DOF from CPMD output
QMATOMS = None  
with open(CPMDOUT, 'r') as f:  
    for line in f:  
        if QMATOMS is None: 
            if ('NR   TYPE        X(BOHR)        Y(BOHR)        Z(BOHR)     MBL' in line):  #start
                QMATOMS = []  # begin saving ATOMS
        elif (' ****************************************************************' in line):  #stop
            break  
        else:  
            QMATOMS.append(line)  # append ATOMS
QMDOF = len(QMATOMS)*3-3 #QM degrees of freedom - no constraints

with open(CPMDOUT, 'r') as f:
    for line in f:
        if('DEGREES OF FREEDOM FOR SYSTEM:' in line):
            TOTDOF = int(line.split()[-1]) #total degrees of freedom
            break
MMDOF = TOTDOF - QMDOF #MM degrees of freedom

#read geometry file to extract element types
elements = []
GeoFile = open(GEOMETRY, 'r')
lineindex = 0
for line in GeoFile:
    if(lineindex==0):
        N = int(line.split()[0])
    if(lineindex>1):
        elements.append(line.split()[0])
    lineindex += 1  
GeoFile.close()

#read trajectory once to read step stored
TrajFile = open(TRAJECTORY, 'r')
Nlines = 1
stepNumbers = []
for line in TrajFile:
    Nlines += 1
    if(int(line.split()[0]) in stepNumbers):
        continue
    else:
        stepNumbers.append(int(line.split()[0]))
TrajFile.close()

#read trajectory a second time to compute Ek and T
df = pd.DataFrame(columns=('step', 'T_QM', 'T_MM', 'T_QMMM'))
TrajFile = open(TRAJECTORY, 'r')
Ek = []
Ek_QM = 0.0
Ek_MM = 0.0
i=0
dfi = 0

'''
Conversion factor included in ConvFactor
    T = sum(0.5 m_i v_i^2)/(0.5 N_dof kB)
    where   m_i amu --> m_i/0.00054857990943 au
            v_i au
            kB = 1.380649E−23J/K= 3.166811563E−6 Ha/K
           
    --> T = Ek/(0.00054857990943*0.5*kB*N_dof) = Ek*ConvFactor/N_dof
'''
kB = 3.166811563e-6
ConvFactor = 1/(0.00054857990943*0.5*kB)

for line in TrajFile:
    if(i%N == 0 and i>0):
        dfi += 1
        T_QM = Ek_QM*ConvFactor/QMDOF
        T_MM = Ek_MM*ConvFactor/MMDOF
        T_QMMM = (Ek_QM+Ek_MM)*ConvFactor/(QMDOF+MMDOF)
        df.loc[dfi] = [stepNumbers[dfi-1], T_QM, T_MM, T_QMMM]
        i=0
        Ek_QM = 0.0
        Ek_MM = 0.0
    if(i<len(QMATOMS)):
        v2 = float(line.split()[-3])**2+float(line.split()[-2])**2+float(line.split()[-1])**2
        Ek_QM += 0.5*get_mass(elements[i])*v2
    elif(i>=len(QMATOMS)):    
        v2 = float(line.split()[-3])**2+float(line.split()[-2])**2+float(line.split()[-1])**2
        Ek_MM += 0.5*get_mass(elements[i])*v2
    i += 1

dfi += 1
df.loc[dfi] = [stepNumbers[dfi-1], T_QM, T_MM, T_QMMM] 
TrajFile.close()

#write output file
if(os.path.isabs(TRAJECTORY)):
    OUTPUT = '/'.join(TRAJECTORY.split('/')[:-1])+'/temp_check.dat'
else:
    OUTPUT = '.'+'/'.join(TRAJECTORY.split('/')[:-1])+'/temp_check.dat'

OutFile = open(OUTPUT, 'w')
for i in range(dfi)[:]:
    OutFile.write(str(df['step'].loc[i+1])+'\t'+str(df['T_QM'].loc[i+1])+'\t'+str(df['T_MM'].loc[i+1])+'\t'+str(df['T_QMMM'].loc[i+1])+'\n') 
OutFile.close()

#print output message
print('------')
print('Temperatures extracted from trajectory saved as:')
print('\t'+OUTPUT)

#plot (optional)
if(args.plot==True):
    import matplotlib.pyplot as plt
    plt.plot(df['step'], df['T_QMMM'], 'o-', color='royalblue', label = 'QM+MM')
    plt.plot(df['step'], df['T_MM'], 'o-', color='mediumseagreen', label = 'MM')
    plt.plot(df['step'], df['T_QM'], 'o-', color='lightcoral', label = 'QM')
    plt.axhline(y = df['T_QM'].mean(), color='indianred', linestyle = '--', label = 'QM mean')
    plt.xlabel('Step')
    plt.ylabel('Temperature (K)')
    plt.plot([], [], ' ', label='(atoms: '+str(len(QMATOMS))+' QM, '+str(N-len(QMATOMS))+' MM)') # reports number of QM and MM atoms
    leg = plt.legend()
    labels = leg.get_texts()
    labels[-1].set_size('small')

    plt.show()
