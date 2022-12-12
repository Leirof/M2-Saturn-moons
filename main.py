import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from datetime import datetime
import matplotlib.animation as animation
from copy import deepcopy as copy

#----------------------------------------------------------------------------------------------------
# Init

start_time = datetime.now()

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
        self.pos_evol[0] = self.pos
        self.vel_evol = np.zeros((T,3))
        self.vel_evol[0] = self.vel
        self.acc_evol = np.zeros((T,3))
        self.acc_evol[0] = self.acc
        Body.all.append(self)

    # def compute_acceleration(self, other_body=None):
    #     if other_body is None:
    #         other_body=Body.all
    #     elif isinstance(other_body, Body):
    #         other_body = [other_body]

    #     acc = np.array([0,0,0])
    #     for other_body in Body.all:
    #         if self.name != other_body.name:
    #             r = other_body.pos - self.pos
    #             d = np.sqrt(np.sum(r**2))
    #             acc = acc + r * other_body.gm / (d**3)
    #     return acc

#----------------------------------------------------------------------------------------------------
# Define the bodies

mimas = Body(
    name = 'Mimas',
    pos  = np.array([-1.813557294137345E+05, 3.492995866271013E+04, 5.054192147215439E+03]),
    vel  = np.array([-2.445920348718876E+00, 1.416736794928916E+01, 1.873213454946043E-02]),
    acc  = np.array([0,0,0]),
    gm   = 2.503489)

thetys = Body(
    name = 'Thetys',
    pos  = np.array([-1.889929479196472E+05, -2.260662471700913E+05, -7.020174976279632E+02]),
    vel  = np.array([8.706811449016753E+00, -7.280896606319645E+00, -2.143307520131600E-01]),
    acc  = np.array([0,0,0]),
    gm   = 41.21)

titan = Body(
    name = 'Titan',
    pos  = np.array([2.956179324256932E+05, 1.212611254175863E+06, -6.683800852915738E+02]),
    vel  = np.array([-5.322398770209340E+00, 1.190373084298060E+00, -3.847254100240810E-02]),
    acc  = np.array([0,0,0]),
    gm   = 8978.14)

saturn = Body(
    name = 'Saturn',
    pos  = np.array([0,0,0]),
    vel  = np.array([0,0,0]),
    acc  = np.array([0,0,0]),
    gm   = 37931206.2)


# for body in [mimas, thetys, titan]:
#     body.acc = body.compute_acceleration(saturn)


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
#         body.vel_evol[t] = body.vel

#----------------------------------------------------------------------------------------------------
# Runge kutta integrator

# Acceleration function : f(t,y) -> (v,a)     with y = (p,v)
def f(t, y, body = None) -> np.array:
    pos = y[:3]
    vel = y[3:]

    gm = body.gm if body is not None else 0

    r = pos
    d = np.sqrt(np.sum(r**2))
    acc = - r * (saturn.gm + gm) / (d**3)

    return np.concatenate((vel, acc))

# Runge kutta 4th order : rk4(t,dt,y,evaluate) -> y     with y = (p,v)
def rk4(t,dt,y,evaluate, body = None) -> np.array:

    if body is thetys and t==dt:
        print(y)

    k1 = evaluate(t, y, body = body)
    if body is thetys and t==dt:
        print("k1", k1)
    k2 = evaluate(t + 0.5*dt, y + 0.5*k1*dt, body = body)
    if body is thetys and t==dt:
        print("k2", k2)
    k3 = evaluate(t + 0.5*dt, y + 0.5*k2*dt, body = body)
    if body is thetys and t==dt:
        print("k3", k3)
    k4 = evaluate(t + dt, y + dt*k3, body = body)
    if body is thetys and t==dt:
        print("k4", k4)
    
    y_new = y + (1./6.)*(k1+ 2.*k2 + 2.*k3 + k4) * dt

    if body is thetys and t==dt:
        print(y_new)

    return y_new


for t in range(1,T):
    
    if t==1:
        print("----------")
        print(t)
    
    for body in [mimas, thetys, titan]:


        y = np.concatenate((body.pos, body.vel))
        if body is thetys and t==1:
            print(y)
        y = rk4(t*DT, DT, y, f, body = body)
        if body is thetys and t==1:
            print(y)

        body.pos = y[:3]
        body.vel = y[3:]

        body.pos_evol[t] = body.pos
        body.vel_evol[t] = body.vel


# body = thetys
# for i in range(4):
#     print("------------------")
#     print(body.pos_evol[i])
#     print(body.vel_evol[i])




#----------------------------------------------------------------------------------------------------
# Plot the results

fig = plt.figure(figsize=(7,7))

plt.subplot(231)
plt.title('Evolution of the x coordinate')
for body in [mimas, thetys,saturn]:
    plt.plot(np.arange(T)*DT/86400.0, body.pos_evol[:,0], label=body.name)
plt.xlabel('Time [days]')
plt.ylabel('x [km]')
plt.legend()

plt.subplot(232)
plt.title('Evolution of the y coordinate')
for body in [mimas, thetys,saturn]:
    plt.plot(np.arange(T)*DT/86400.0, body.pos_evol[:,1], label=body.name)
plt.xlabel('Time [days]')
plt.ylabel('y [km]')
plt.legend()

plt.subplot(233)
plt.title('Evolution of the z coordinate')
for body in [mimas, thetys,saturn]:
    plt.plot(np.arange(T)*DT/86400.0, body.pos_evol[:,2], label=body.name)
plt.xlabel('Time [days]')
plt.ylabel('z [km]')
plt.legend()

plt.subplot(234)
plt.title('Evolution of the distance saturn - thetys')
plt.plot(np.arange(T), np.sqrt(np.sum((saturn.pos_evol - thetys.pos_evol)**2, axis=1)), label='Distance')
plt.xlabel('Time [days]')
plt.ylabel('Distance [km]')
plt.legend()

ax = plt.subplot(235, projection='3d')
ax.set_title('Orbits of the bodies')
for body in [mimas, thetys,saturn]:
    ax.plot(body.pos_evol[:,0], body.pos_evol[:,1], body.pos_evol[:,2], marker='o', label=body.name)

end_time = datetime.now()
ex_time = end_time - start_time
print('Running time : ',ex_time)

plt.show()
