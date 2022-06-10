import sys
import os.path

#help function accessible via python3 traj_xyz_convert.py --help (or -h) 
import argparse
parser = argparse.ArgumentParser(description='Generate TRAJEC_converted.xyz file from TRAJECTORY to easily visualize in VMD.')
parser.add_argument('TRAJECTORY_PATH', type=str,
                    help='path to TRAJECTORY file')
parser.add_argument('GEOMETRY_PATH', type=str, 
                    help='path to GEOMETRY.xyz file')
args = parser.parse_args()

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

#read geometry file to extract atom elements
GeoFile = open(GEOMETRY, 'r')
element_list = []
linecounter = 0
for line in GeoFile:
    if(linecounter>1):
        element_list.append(line.split()[0])
    linecounter+=1

Natoms = len(element_list)

#read configurations and velocities from trajectory  file
#and save them in output file with correct formatting 
TrajFile = open(TRAJECTORY, 'r')

if(os.path.isabs(TRAJECTORY)):
    OUTPUT = '/'.join(TRAJECTORY.split('/')[:-1])+'/TRAJECTORY_converted.xyz'
else:
    OUTPUT = '.'+'/'.join(TRAJECTORY.split('/')[:-1])+'/TRAJECTORY_converted.xyz'
OutFile = open(OUTPUT, 'w')


linecounter = 0
for line in TrajFile:
    print(linecounter, line)
    if(linecounter%Natoms==0):
        linecounter = 0
        OutFile.write('       '+str(Natoms)+'\n')
        OutFile.write(' STEP:           '+str(line.split()[0])+'\n')
    #TRAJECTORY.xyz format: X      x_coord      y_coord      z_coord')
    OutFile.write(str(element_list[linecounter])+'\t'+str(line.split()[1])+'\t'+str(line.split()[2])+'\t'+str(line.split()[3]+'\n'))
    
    linecounter += 1

TrajFile.close()
OutFile.close()

#print output message
print('------')
print('TRAJECTORY file converted in xyz format saved as:')
print('\t'+OUTPUT)
