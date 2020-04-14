import os.path
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


fileE = 'energy.txt'
outputE = open(fileE,'r')
#  infile.readline() # skip the first line not first line here (yet)
tijd = []
relEnergy = []
for line in outputE:
    words = line.split()
    # words[0]: tijd, words[1]: relEnergy
    tijd.append(float(words[0]))
    relEnergy.append(float(words[1]))

outputE.close()
plt.figure(101)
plt.plot(tijd,relEnergy,'-')
plt.xlabel(' $t$ ',size=16)
plt.ylabel(' $|H(t)-H(0)|/H(0)$ ',size=16)
# plt.axes([0,10,0.001,0.002])
plt.yticks([0.0010082,0.0010083,0.0010084,0.0010085,0.0010086])
# plt.axes(xlim=(0, 10), ylim=(0.0010082, 0.0010086), autoscale_on=False) # [0 10 0.0010082 0.0010086])

plt.show(block=True)
plt.pause(0.001)
plt.gcf().clear()
plt.show(block=False)

print("Finished program!")

