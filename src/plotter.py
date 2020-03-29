import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import argparse



parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-N', '--number' , type=int, nargs='+',
                    help='number of csv files')


args = parser.parse_args()


N = int( args.number[0])

fig = plt.figure()


circle1 = plt.Circle((0, 0), 1, color='black')
plt.gcf().gca().add_artist(circle1)
plt.axis('equal')
plt.xlim(-5 ,5)
plt.ylim(-5,10)

for i in range(0,N):
    db = pd.read_csv( 'files/photon%d.csv'%i)  
    print(db)
    plt.plot( db['r']*np.sin(db['h'])*np.cos(db['p']),db['r']*np.cos(db['h'])  )


plt.show()