import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import argparse
from mpl_toolkits.mplot3d import Axes3D


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-N', '--number' , type=int, nargs=1, help='number of csv files')
parser.add_argument('-M', '--mode', nargs=1, help='mode: example xy or xyz')


args = parser.parse_args()

print(args)

N = int( args.number[0])
name = args.mode[0][1:]

fig = plt.figure()

if name == 'xz': 
    print("ddddddddddddd")
    plt.axis('equal')
    plt.xlim(-5 ,5)
    plt.ylim(-5,10)

    for i in range(0,N):
        db = pd.read_csv( 'files/photon%d.csv'%i, names=['t','x','y','z','pt','px','py','pz'])  
        #print(db)
        plt.plot( db['x'],db['z']  )
if name == 'xy': 
    plt.axis('equal')
    plt.xlim(-5 ,5)
    plt.ylim(-5,10)

    for i in range(0,N):
        db = pd.read_csv( 'files/photon%d.csv'%i, names=['t','x','y','z','pt','px','py','pz'])  
        #print(db)
        plt.plot( db['x'],db['y']  )
elif name == 'xyz':
    ax = fig.add_subplot(111, projection='3d')

    #plt.axis('equal')
    plt.xlim(-5 ,5)
    plt.ylim(-5,5)
    ax.set_zlim(-5,5)

    for i in range(0,N):
        db = pd.read_csv( 'files/photon%d.csv'%i, names=['t','x','y','z','pt','px','py','pz'])  
        #print(db)
        ax.plot( db['x'], db['y'],db['z'])


plt.show()