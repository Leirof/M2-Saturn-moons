import numpy as np
import matplotlib.pyplot as plt

def pos_dyn(bodies, T, DT):
    fig = plt.figure(figsize=(15,5))

    # Evolution of the x coordinate
    plt.subplot(131)
    plt.title('Evolution of the x coordinate')
    for body in bodies:
        plt.plot(np.arange(T)*DT/86400.0, body.pos_dyn[:,0], label=body.name, c=body.color)
    plt.xlabel('Time [days]')
    plt.ylabel('x [km]')
    plt.legend()

    # Evolution of the y coordinate
    plt.subplot(132)
    plt.title('Evolution of the y coordinate')
    for body in bodies:
        plt.plot(np.arange(T)*DT/86400.0, body.pos_dyn[:,1], label=body.name, c=body.color)
    plt.xlabel('Time [days]')
    plt.ylabel('y [km]')
    plt.legend()

    # Evolution of the z coordinate
    plt.subplot(133)
    plt.title('Evolution of the z coordinate')
    for body in bodies:
        plt.plot(np.arange(T)*DT/86400.0, body.pos_dyn[:,2], label=body.name, c=body.color)
    plt.xlabel('Time [days]')
    plt.ylabel('z [km]')
    plt.legend()
    plt.show()

def orbits_3D(bodies):
    fig = plt.figure(figsize=(5,5))

    ax = plt.subplot(111, projection='3d')
    ax.set_title('Orbits of the bodies')
    for body in bodies:
        ax.plot(body.pos_dyn[:,0], body.pos_dyn[:,1], body.pos_dyn[:,2], marker='o', label=body.name, c=body.color)
    
    plt.show()

def orbits(bodies):
    plt.figure(figsize=(15,5))

    # Orbit on plane x,y
    plt.subplot(131)
    plt.title('Orbit on plane x,y')
    for body in bodies:
        plt.plot(body.pos_dyn[:,0], body.pos_dyn[:,1], label=body.name, c=body.color) 
    plt.xlabel('x [km]')
    plt.ylabel('y [km]')
    plt.legend()

    # Orbit on plane x,z
    plt.subplot(132)
    plt.title('Orbit on plane x,z') 
    for body in bodies:
        plt.plot(body.pos_dyn[:,0], body.pos_dyn[:,2], label=body.name, c=body.color)
    plt.xlabel('x [km]')
    plt.ylabel('z [km]')    
    plt.legend()

    # Orbit on plane y,z
    plt.subplot(133)

    plt.title('Orbit on plane y,z')
    for body in bodies:
        plt.plot(body.pos_dyn[:,1], body.pos_dyn[:,2], label=body.name, c=body.color)
    plt.xlabel('y [km]')
    plt.ylabel('z [km]')
    plt.legend()
   
    plt.show()

def vel_dyn(bodies, T, DT):
    fig = plt.figure(figsize=(15,5))

    # Evolution of the x velocity
    plt.subplot(131)
    plt.title('Evolution of the x velocity')
    for body in bodies:
        plt.plot(np.arange(T)*DT/86400.0, body.vel_dyn[:,0], label=body.name, c=body.color)
    plt.xlabel('Time [days]')
    plt.ylabel('Vx [km/s]')
    plt.legend()

    # Evolution of the y velocity
    plt.subplot(132)
    plt.title('Evolution of the y velocity')
    for body in bodies:
        plt.plot(np.arange(T)*DT/86400.0, body.vel_dyn[:,1], label=body.name, c=body.color)
    plt.xlabel('Time [days]')
    plt.ylabel('Vy [km/s]')
    plt.legend()

    # Evolution of the z velocity
    plt.subplot(133)
    plt.title('Evolution of the z velocity')
    for body in bodies:
        plt.plot(np.arange(T)*DT/86400.0, body.vel_dyn[:,2], label=body.name, c=body.color)
    plt.xlabel('Time [days]')
    plt.ylabel('Vz [km/s]')
    plt.legend()
    plt.show()

def orbital_elements(bodies, T, DT):
    plt.figure(figsize=(30,5*len(bodies)))
    
    for i, body in enumerate(bodies):
        plt.subplot(len(bodies),6,i*6+1)
        plt.plot(np.arange(T)*DT/86400.0, body.orb_dyn[:,0], label=body.name, c=body.color)
        plt.title('Evolution of the semimajor axis')
        plt.xlabel('Time [days]')
        plt.ylabel('a [km]')
        plt.legend()

        plt.subplot(len(bodies),6,i*6+2)
        plt.plot(np.arange(T)*DT/86400.0, body.orb_dyn[:,1], label=body.name, c=body.color)
        plt.title('Evolution of the eccentricity')
        plt.xlabel('Time [days]')
        plt.ylabel('e')
        plt.legend()

        plt.subplot(len(bodies),6,i*6+3)
        plt.plot(np.arange(T)*DT/86400.0, body.orb_dyn[:,2], label=body.name, c=body.color)
        plt.title('Evolution of the inclination')
        plt.xlabel('Time [days]')   
        plt.ylabel('i [deg]')
        plt.legend()

        plt.subplot(len(bodies),6,i*6+4)
        plt.plot(np.arange(T)*DT/86400.0, body.orb_dyn[:,3], label=body.name, c=body.color)  
        plt.title('Evolution of the longitude of the ascending node')
        plt.xlabel('Time [days]')
        plt.ylabel(r'$\Omega$ [deg]')   
        plt.legend()    

        plt.subplot(len(bodies),6,i*6+5)
        plt.plot(np.arange(T)*DT/86400.0, body.orb_dyn[:,4], label=body.name, c=body.color)
        plt.title('Evolution of the argument of periapsis')
        plt.xlabel('Time [days]')
        plt.ylabel(r'$\varpi$ [deg]')   
        plt.legend()

        plt.subplot(len(bodies),6,i*6+6)  
        plt.plot(np.arange(T)*DT/86400.0, body.orb_dyn[:,5], label=body.name, c=body.color)
        plt.title('Evolution of the mean anomaly')  
        plt.xlabel('Time [days]')
        plt.ylabel(r'$\lambda$ [deg]')
        plt.legend()

    plt.show()
