import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from datetime import datetime
import matplotlib.animation as animation
#----------------------------------------------------------------------------------
time_1 = datetime.now()

# Data extracted from the horizon results
DT = 600 # in seconds
T = 2500

#----------------------------------------------------------------------------------------------------
# Define the class body

class Body():

    all = []

    def __init__(self, name, pos, vel, acc, gm):
        self.name = name
        self.pos = pos
        self.vel = vel
        self.acc = acc
        self.gm = gm
        self.pos_evol = np.zeros((T,3))
        self.vel_evol = np.zeros((T,3))
        self.acc_evol = np.zeros((T,3))
        Body.all.append(self)

    def compute_acceleration(self, other_body=None):
        if other_body is None:
            other_body=Body.all
        elif isinstance(other_body, Body):
            other_body = [other_body]

        acc = np.array([0,0,0])
        for other_body in Body.all:
            if self.name != other_body.name:
                d = other_body.pos - self.pos
                r = np.sqrt(np.sum(d**2))
                acc = acc + d * other_body.gm / (r**3)
        return acc

#----------------------------------------------------------------------------------------------------
# Define the bodies

mimas = Body(
    name = 'Mimas',
    pos  = np.array([-5.810375154371507E+07, 1.457217395207454E+08, -7.653232466058011E+07])/1000,
    vel  = np.array([-1.412505918711402E+04, 1.326791568773807E+03, 6.383898069507170E+02])/1000,
    acc  = np.array([0,0,0]),
    gm   = 2.503489)

thetys = Body(
    name = 'Thetys',
    pos  = np.array([8.483395081361104E+07, -2.451844175469960E+08, 1.265076985489711E+08])/1000,
    vel  = np.array([8.409967324758574E+03, 8.936994802880729E+03, -5.540300346589594E+03])/1000,
    acc  = np.array([0,0,0]),
    gm   = 41.21)

titan = Body(
    name = 'Titan',
    pos  = np.array([1.080416731555666E+09, 5.169709280480295E+08, -3.742024710406277E+08])/1000,
    vel  = np.array([-5.660865246081125E+03, 1.298776350694247E+04, -6.211781636886458E+03])/1000,
    acc  = np.array([0,0,0]),
    gm   = 8978.14)

saturn = Body(
    name = 'Saturn',
    pos  = np.array([0,0,0]),
    vel  = np.array([0,0,0]),
    acc  = np.array([0,0,0]),
    gm   = 37931206.2)

#----------------------------------------------------------------------------------------------------
# Define the function to calculate the acceleration

for body in [mimas, thetys, titan]:
    body.acc = body.compute_acceleration(saturn)

for t in range(0,T):
    
    for body in [mimas, thetys, titan]:
        # print("\n----------\n")
        # print("At t:")
        # print(f"   x  = {p[0]:.3e}, y  = {p[1]:.3e}, z  = {p[2]:.3e}")
        # print(f"   vx = {v[0]:.3e}, vy = {v[1]:.3e}, vz = {v[2]:.3e}")
        # print(f"   ax = {a[0]:.3e}, ay = {a[1]:.3e}, az = {a[2]:.3e}")

        body.pos = body.pos + body.vel * DT + body.acc/2 * DT**2
        body.vel =            body.vel      + body.acc/2 * DT
        body.acc = body.compute_acceleration(saturn)
        body.vel =            body.vel      + body.acc/2 * DT

        body.pos_evol[t] = body.pos

fig = plt.figure(figsize=(7,7))
# ax = plt.subplot(121, projection='3d')
# plt.title('Orbits of the bodies')
# plt.grid()
# plt.scatter(0,0,0, marker='o', color='r')
# ax.scatter(mimas.pos_evol[0,:,0], mimas.pos_evol[0,:,2], c = acc_xz[0,:], cmap = "gnuplot")
# plt.scatter(mimas.pos_evol[:,0], mimas.pos_evol[:,1], mimas.pos_evol[:,2], marker='o', color='b')
# plt.scatter(thetys.pos_evol[:,0], thetys.pos_evol[:,1], thetys.pos_evol[:,2], marker='o', color='g')
# plt.scatter(titan.pos_evol[:,0], titan.pos_evol[:,1], titan.pos_evol[:,2], marker='o', color='y')

plt.subplot(111)
plt.title('Evolution of the x coordinate')
for body in [mimas, thetys, titan]:
    plt.plot(np.arange(T)*DT/86400.0, body.pos_evol[:,0], label=body.name)

plt.legend()
plt.show()

time_2 = datetime.now()
ex_time = time_2-time_1
print('Running time : ',ex_time)