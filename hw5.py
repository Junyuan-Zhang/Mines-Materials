import numpy as np
import matplotlib.pyplot as plt
from scipy import signal as sg

# a = [0,0,1,0,3]
# b = [2,1,2]
# c = sg.convolve(a,b, mode = 'same')
# print(c)

def Dirichlet_bc(bcl,bcr,Q,end,H):
    si = 1
    M =200
    dz = si/M
    dt = (dz**2/(2*2.5))
    a = dt/(dz**2)
    #define convolution kernel
    K = [a, 1-2*a, a]
    tt = np.zeros(int(end/dt+1))
    C0 = np.ones(int(1/dz+1))
    #initial condition
    zz = 0
    for i in range(len(C0)):
        C0[i] = 1 - 1*zz
        zz = zz+dz
    Ct = []
    axis = np.arange(0, si+dz, dz)
    ss  = np.zeros(len(axis))
    # plt.plot(axis, C0, label=f'Initial condition with H={H}', linewidth=1)
    for i in range(0, int(end/dt+1)):  
        Cnew = sg.convolve(C0,K, mode = 'same') + dt*Q
        Cnew[0] = bcl
        Cnew[-1] = bcr
        tt[i] = tt[i-1]+dt
        autoplt = round(tt[i],6)
        C0 = Cnew
        # # activate if plot the dimensionless concentration vs dimeensionless distance
        # if autoplt== end: 
        #     plt.plot(axis, C0, label=f'Dimension time t={autoplt}, H={H}', linewidth=1)
        #      for j in range(len(axis)):
        #         ss[j] = -H*(axis[j]**2)/2+((H/2)-1)*axis[j]+1
        #     plt.plot(axis, ss, label='Analytical solution in steady state', linewidth=1)



    #     # # make every point at x=1 in the timescale into a list 
    #     Ct.append(Cnew[axis == 0.5])
    # #transform list to array(To calculate)
    # Ctt = np.array(Ct)
    # plt.plot(tt,Ctt, label=f'Mid depth temperature at H={H}')
    return C0


def Neumann_bc(ic,D,k,Q,Qrad,bcr,si,end,L):
    dx = 1e-2
    dt = 1e-5
    QQ = Q*k/D
    a = dt/(dx**2)
    #define convolution kernel
    K = [a, 1-2*a, a]
    tt = np.zeros(int(end/dt+1))
    C0 = np.zeros(int(si/dx+1))*ic
    C0[-1] = bcr
    Ct = []
    phi = (dx*Qrad*L)/(k)
    axis = np.arange(0, si+dx, dx)
    for i in range(0, int(end/dt+1)):  
        Cnew = sg.convolve(C0,K, mode = 'same')+dt*QQ
        Cnew[0] = ((2*dt)/dx**2)*(phi + C0[1] - C0[0])+C0[0]+dt*QQ
        Cnew[-1] = bcr
        tt[i] = tt[i-1]+dt
        autoplt = round(tt[i],6)
        C0 = Cnew
        # activate if plot the dimensionless concentration vs dimeensionless distance
        if autoplt== end:
            plt.plot(axis, C0, label=f'After {end*1000} years temperature')
        #     DD = 1e-2
        #     plt.plot(axis, C0, label=f'dimension time t={autoplt}')
        #     q = -DD*(Cf/xw)*((C0[1]-C0[0])/dx) 
        #     print(q)
        # make every point at x=1 in the timescale into a list 
    #     Ct.append(Cnew[axis == 0])
    # #transform list to array(To calculate)
    # Ctt = np.array(Ct)
    # #convert dimensionless value to real value
    # act_tt = tt*(xw**2/D)/86400
    # act_c = Ctt*Cf
    # plt.plot(act_tt,act_c,label=f'Concentration(ppm) at the well with D={D}')
    # #print the time when c reach 100ppm
    # act_c = np.around(act_c, 2)            
#1day = 3600*24=86400s
def Dirichlet_bc_special(bcl,bcr,Q,end,H,Cinitial):
    si = 1
    M =200
    dz = si/M
    dt = (dz**2/(2*2.5))
    a = dt/(dz**2)
    #define convolution kernel
    K = [a, 1-2*a, a]
    tt = np.zeros(int(end/dt+1))
    axis = np.arange(0, si+dz, dz)
    for i in range(len(Cinitial)):
        if axis[i] > 19/40 and axis[i] < 21/40:
            Cinitial[i] = 1.2
    C0 = Cinitial
    Ct = []
    Tm = np.ones(len(axis))
    ss  = np.zeros(len(axis))
    if end == 0.01:
        plt.plot(axis, C0, label=f'Initial condition with H={H}', linewidth=1)
        plt.plot(axis, Tm, label='Tm', linewidth=1)
    for i in range(0, int(end/dt+1)):  
        Cnew = sg.convolve(C0,K, mode = 'same') + dt*Q
        Cnew[0] = bcl
        Cnew[-1] = bcr
        tt[i] = tt[i-1]+dt
        autoplt = round(tt[i],6)
        C0 = Cnew
                # # activate if plot the dimensionless concentration vs dimeensionless distance
        if autoplt== end: 
            plt.plot(axis, C0, label=f'Dimension time t={autoplt}, H={H}', linewidth=1)
        #      for j in range(len(axis)):
        #         ss[j] = -H*(axis[j]**2)/2+((H/2)-1)*axis[j]+1
        #     plt.plot(axis, ss, label='Analytical solution in steady state', linewidth=1)



    #     # # make every point at x=1 in the timescale into a list 
    #     Ct.append(Cnew[axis == 0.5])
    # #transform list to array(To calculate)
    # Ctt = np.array(Ct)
    # plt.plot(tt,Ctt, label=f'Mid depth temperature at H={H}')
# #bottom boundary condition value 
# T0=1
# #ceil boundary condition value
# T1=0
# #dimensionless Q
# Q=1
# #H value
# H1=2
# H2=1
# H3=4
# #guessing a time that reach steady state
# end = 1

# Cinitial = Dirichlet_bc(T0,T1,Q*H1,end,H1)
# Dirichlet_bc_special(T0,T1,Q*H1,0.00001,H1,Cinitial)
# Dirichlet_bc_special(T0,T1,Q*H1,0.0001,H1,Cinitial)
# Dirichlet_bc_special(T0,T1,Q*H1,0.01,H1,Cinitial)
#diffusion coefficient
D = 5e-4 #m2/s
k=3  #J/(m-deg-s)
Q=2.25e-7 #J/m3
Qrad = 1.36e-5 #J/(s-m^2)
ic = 0
bcr = 0
si = 1
L = 1
end = 1#actually 1000 years
Neumann_bc(ic,D,k,Q,Qrad,bcr,si,1,L)
#Neumann_bc(ic,D,k,Q,Qrad,bcr,si,0.8,L)
plt.xlabel('Dimensionless depth')
#plt.xlabel('Dimensionless time')
plt.ylabel('Dimensionless temperature')
plt.legend()

#plt.savefig('fig1.png')
plt.show()
