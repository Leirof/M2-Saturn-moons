import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from datetime import datetime
import matplotlib.animation as animation
#----------------------------------------------------------------------------------
time_1 = datetime.now()

G = 6.6742e-11
MS = 5.6834e26
M_Mimas = 3.7493e19
M_Thetys = 6.17449e20
M_Titan = 1.3452e23
dt = 86400.0

nt = 20
S = np.zeros((nt,9))
S1 = np.zeros((nt,9))
S2 = np.zeros((nt,9))
S3 = np.zeros((nt,9))

Mimas_pos = np.zeros(3)
Mimas_vel = np.zeros(3)
Mimas_acc = np.zeros(3)
Mimas_acc_new = np.zeros(3)

Thetys_pos = np.zeros(3)
Thetys_vel = np.zeros(3)
Thetys_acc = np.zeros(3)
Thetys_acc_new = np.zeros(3)

Titan_pos = np.zeros(3)
Titan_vel = np.zeros(3)
Titan_acc = np.zeros(3)
Titan_acc_new = np.zeros(3)

Mimas_pos[:] = [-5.810375154371507E+07, 1.457217395207454E+08, -7.653232466058011E+07]
Mimas_vel[:] = [-1.412505918711402E+04, 1.326791568773807E+03, 6.383898069507170E+02]
Mimas_acc[:] = [0,0,0]
Mimas_acc_new[:] = [0,0,0]

Thetys_pos[:] = [8.483395081361104E+07, -2.451844175469960E+08, 1.265076985489711E+08]
Thetys_vel[:] = [8.409967324758574E+03, 8.936994802880729E+03, -5.540300346589594E+03]
Thetys_acc[:] = [0,0,0]
Thetys_acc_new[:] = [0,0,0]

Titan_pos[:] = [1.080416731555666E+09, 5.169709280480295E+08, -3.742024710406277E+08]
Titan_vel[:] = [-5.660865246081125E+03, 1.298776350694247E+04, -6.211781636886458E+03]
Titan_acc[:] = [0,0,0]
Titan_acc_new[:] = [0,0,0]

Mimas_data = np.zeros((nt,3))

for t in range(1,nt):
    Mimas_acc_new += G*(MS*M_Mimas)/(np.linalg.norm(Mimas_pos)**2)*Mimas_pos
    Mimas_vel += (Mimas_acc + Mimas_acc_new)*0.5*dt
    Mimas_acc += Mimas_acc_new
    Mimas_pos += Mimas_vel*dt + Mimas_acc*0.5*dt**2
    Mimas_data[t] = Mimas_pos

fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111,projection='3d')
ax.scatter(0,0,0,s=10,c='r')
ax.scatter(Mimas_data[:,0],Mimas_data[:,1],Mimas_data[:,2],s=5)
plt.show()

time_2 = datetime.now()
ex_time = time_2-time_1
print('Running time : ',ex_time)