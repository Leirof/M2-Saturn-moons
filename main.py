import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from datetime import datetime
import matplotlib.animation as animation
from copy import deepcopy as copy

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
                r = other_body.pos - self.pos
                d = np.sqrt(np.sum(r**2))
                acc = acc + r * other_body.gm / (d**3)
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


for body in [mimas, thetys, titan]:
    body.acc = body.compute_acceleration(saturn)

#----------------------------------------------------------------------------------------------------
# Leap frog integrator

# for t in range(0,T):
    
#     for body in [mimas, thetys]:

#         body.pos = body.pos + body.vel * DT + body.acc/2 * DT**2
#         body.vel =            body.vel      + body.acc/2 * DT
#         body.acc = body.compute_acceleration(saturn)
#         if body is titan:
#             body.acc *= 5.5
#         body.vel =            body.vel      + body.acc/2 * DT

#         body.pos_evol[t] = body.pos

#----------------------------------------------------------------------------------------------------
# Runge kutta integrator

# Acceleration function : f(t,y) -> (v,a)     with y = (p,v)
def f(t, y) -> np.array:
    pos = y[:3]
    vel = y[3:]

    other_body=saturn

    acc = np.array([0,0,0])
    for other_body in Body.all:
        if body.name != other_body.name:
            r = other_body.pos - pos
            d = np.sqrt(np.sum(r**2))
            acc = acc + r * other_body.gm / (d**3)

    return np.concatenate((vel, acc))

# Runge kutta 4th order : rk4(t,dt,y,evaluate) -> y     with y = (p,v)
def rk4(t,dt,y,evaluate) -> np.array:

    k1 = dt * evaluate(t, y) 
    k2 = dt * evaluate(t + 0.5*dt, y + 0.5*k1)
    k3 = dt * evaluate(t + 0.5*dt, y + 0.5*k2)
    k4 = dt * evaluate(t + dt, y + k3)
    
    y_new = y + (1/6.)*(k1+ 2*k2 + 2*k3 + k4)

    return y_new


for t in range(0,T):
    
    for body in [mimas, thetys]:

        y = np.concatenate((body.pos, body.vel))
        y = rk4(t, DT, y, f)

        body.pos = y[:3]
        body.vel = y[3:]

        body.pos_evol[t] = body.pos


#----------------------------------------------------------------------------------------------------
# Plot the results

fig = plt.figure(figsize=(7,7))

plt.subplot(231)
plt.title('Evolution of the x coordinate')
for body in Body.all:
    plt.plot(np.arange(T)*DT/86400.0, body.pos_evol[:,0], label=body.name)
plt.xlabel('Time [days]')
plt.ylabel('x [km]')
plt.legend()

plt.subplot(232)
plt.title('Evolution of the y coordinate')
for body in Body.all:
    plt.plot(np.arange(T)*DT/86400.0, body.pos_evol[:,1], label=body.name)
plt.xlabel('Time [days]')
plt.ylabel('y [km]')
plt.legend()

plt.subplot(233)
plt.title('Evolution of the z coordinate')
for body in Body.all:
    plt.plot(np.arange(T)*DT/86400.0, body.pos_evol[:,2], label=body.name)
plt.xlabel('Time [days]')
plt.ylabel('z [km]')
plt.legend()

ax = plt.subplot(235, projection='3d')
ax.set_title('Orbits of the bodies')
for body in Body.all:
    ax.plot(body.pos_evol[:,0], body.pos_evol[:,1], body.pos_evol[:,2], marker='o', label=body.name)

plt.show()

time_2 = datetime.now()
ex_time = time_2-time_1
print('Running time : ',ex_time)