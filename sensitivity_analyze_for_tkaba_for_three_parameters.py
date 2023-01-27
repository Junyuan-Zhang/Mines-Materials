from io import RawIOBase
import numpy as np
from matplotlib import pyplot as plt
from scipy import special as sp
from SALib.sample import saltelli
from SALib.analyze import sobol
# package of Modified Bessel function of the second kind of real order v  is kv(v, z)
# In the calculation, use zero order, so the function is sp.kv(0,z)

plt.rcParams['mathtext.default'] = 'regular'
def dimensionr(B1, raba):
    aa = (1/B1)*(np.sqrt(0.1))*raba
    return aa

def omegaD(Bimp, raba, Kaba, alpha, gamma, Q, B1):
    omega = Bimp/((raba**2)*Kaba)
    bb = alpha*gamma*Q/B1
    res = omega*bb
    return res

# alpha1 and alpha2 will need in sensitivity analysis


# def alpha1(rinf, B1):
#     cc = rinf/B1
#     return cc


# def alpha2(K1h, K2h):
#     dd = K1h/K2h
#     return dd

# def alpha3(K1h, B1, Q):
#     ee = K1h*(B1**2)/Q


def frD(rD):
    f = ((-np.log(rD)+sp.kv(0, np.pi*rD))*(1/(2*np.pi)))
    return f

def dimensionlog(rinf, rcoe):
    rcoeD = np.log(rinf/rcoe)
    return rcoeD
#(Bimp, raba, Kaba, Q, alpha, beta, gamma)
def coe_E(Bimp, raba, Kaba, Q, alpha, beta, gamma):
    rinf = 1000#fixed
    B1 = rinf/alpha
    rD = dimensionr(B1, raba)
    oD = omegaD(Bimp, raba, Kaba, alpha, gamma, Q, B1)
    f = frD(rD)
    rabaD = dimensionlog(rinf, raba)
    E1 = np.pi*rabaD
    E2 = (alpha**2)*beta*gamma
    E3 = (oD+alpha*f)**2
    E = E1*E2/E3
    return E


def coe_F(Bimp, raba, Kaba, Q, alpha, gamma,rinj):
    rinf = 1000
    B1 = rinf/alpha
    rD = dimensionr(B1, raba)
    oD = omegaD(Bimp, raba, Kaba, alpha, gamma, Q, B1)
    f = frD(rD)
    rinjD = dimensionlog(rinf, rinj)
    F1 = rinjD/2
    F2 = alpha/(oD+alpha*f)
    F = -F1*F2
    return F


def final_sol(rinj, Bimp, raba, Q, alpha, beta, gamma, tKaba):
    Kaba = 10**tKaba
    E=coe_E(Bimp, raba, Kaba, Q, alpha, beta, gamma)
    F=coe_F(Bimp, raba, Kaba, Q, alpha, gamma,rinj)
    sol = 0.5*E-F-0.5*np.sqrt((E**2)-4*E*F)
    return sol

problem = {
    'num_vars':3,
    'names':['alpha','beta','gamma'],
    'bounds':[[10,1000],[0.01,100],[1e-4,100]]
}

MS = 18
BS = 18

plt.rc('font', size=BS,family='Times New Roman')          # controls default text sizes
plt.rc('axes', titlesize=BS)     # fontsize of the axes title
plt.rc('axes', labelsize=BS)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MS)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MS)    # fontsize of the tick labels
plt.rc('legend', fontsize=MS)    # legend fontsize
plt.rc('figure', titlesize=BS)

param_values = saltelli.sample(problem, 2**14)

# set input
# thickness of impermeable layer
Bimp = 10
# radius of well
raba = 0.15
# Hydraulic conductivity of abadoned well
tKaba = np.arange(0,3,0.01)
# tKaba = np.arange(2,3,0.01)
#Kimp = 100
#np.arange(10,10000,1)
# approach to infinite distance
# rinf = 1000
# Hrizontal Hydraulic conductivity of upper layer(unconfined aquifer)
#aa2=K1h/K2h
#aa2 = 3
# injection rate
Q = 1000
# thickness of lower layer(confined aquifer)
#aa1=rinf/B1
#aa1 = np.arange(1,10,0.1)
# # Hrizontal Hydraulic conductivity of lower layer(confined aquifer)
# K1h = 1
# # distance from injection well to abandoned well
rinj = 100
# # Vertical Hydraulic conductivity of lower layer(confined aquifer)
# Kz = 0.1*K1h

#final_sol(Bimp, raba, Q, B1, alpha, beta, gamma, tKaba)
q = np.array([final_sol(rinj, Bimp, raba, Q, *params, tKaba) for params in param_values])
#Q=final_sol(Bimp, raba, Kimp, rinf, K2h, Q, B1, K1h, rinj, Kz)

sobol_indices = [sobol.analyze(problem, Q) for Q in q.T]


'''
aa1=100
aa2=0.1

Q = final_sol(Bimp, raba, Kimp, rinf, Q, K1h, rinj, Kz, aa1, aa2)
plt.plot(rinj, Q, label='method f(rD)')

plt.title("Leakage rate change with Confined aquifer thickness")
plt.xlabel("Confined aquifer thickness(m)")
plt.ylabel("leakage rate($m^3/s$)")
plt.legend()
plt.show()
'''


S1s = np.array([s['S1'] for s in sobol_indices])
S2sa = np.array([s['S2'][0,1] for s in sobol_indices])
S2sb = np.array([s['S2'][0,2] for s in sobol_indices])
S2sc = np.array([s['S2'][1,2] for s in sobol_indices])

fig = plt.figure(constrained_layout=True, dpi=150)


for i in range(3):
    if i==0:
        plt.plot(tKaba, S1s[:, i],
            label=rf'$\alpha$',linewidth = 5)
    elif i==1:
        plt.plot(tKaba, S1s[:, i],
            label=rf'$\beta$',linewidth = 3)
    else:
        plt.plot(tKaba, S1s[:, i],
            label=rf'$\gamma$',linewidth = 5)


plt.legend(loc='upper right')
plt.ylabel('First-order indices', fontname = 'Times New Roman') 
plt.xlabel('$log_{10}(K_{aba})$', fontname = 'Times New Roman')  
plt.grid(True)
#plt.title('First-order indices of each variables relative to $B_{imp}$(m)')
fig = plt.figure(constrained_layout=True, dpi=150)
plt.plot(tKaba, S2sa, label=rf'$\alpha$ and $\beta$',linewidth = 5)
plt.plot(tKaba, S2sb, label=rf'$\alpha$ and $\gamma$',linewidth = 5)
plt.plot(tKaba, S2sc, label=rf'$\beta$ and $\gamma$',linewidth = 5)
#gs = fig.add_gridspec(3, 1)
plt.legend(loc='upper right')
plt.ylabel('Second-order indices', fontname = 'Times New Roman') 
plt.xlabel('$log_{10}(K_{aba})$', fontname = 'Times New Roman')  
plt.grid(True)
#ax0 = fig.add_subplot(gs[:, 0])


'''
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax3 = fig.add_subplot(gs[2])

for i, ax in enumerate([ax1, ax2, ax3]):
    ax.plot(raba, S1s[:, i],
            label=r'S1$_\mathregular{{{}}}$'.format(problem["names"][i]),
            color='black')
    ax.set_xlabel("raba")
    ax.set_ylabel("First-order Sobol index")

    ax.set_ylim(0, 1.04)

    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()

    ax.legend(loc='upper right')
'''

#ax0.plot(raba, np.mean(q, axis=0), label="Mean", color='black')

# in percent
#prediction_interval = 95

#ax0.fill_between(raba,
#                 np.percentile(q, 50 - prediction_interval/2., axis=0),
#                 np.percentile(q, 50 + prediction_interval/2., axis=0),
#                 alpha=0.5, color='black',
#                 label=f"{prediction_interval} % prediction interval")

#ax0.set_xlabel("x")
#ax0.set_ylabel("y")
#ax0.legend(title=r"$y=a+b\cdot x^2$",
#           loc='upper center')._legend_box.align = "left"

plt.show()
