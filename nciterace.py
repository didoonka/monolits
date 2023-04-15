# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 18:12:02 2023

@author: DiDoonka
"""
import csv
import numpy as np
import matplotlib.pyplot as plt

VM = 0.000479293
VD = 0.00104267
A = -0.02479
B = 0.00986
C = 0.03214
Kf = 0.0000000000000000867
L = 0.189
ID = 0.075
p = 25
vis = 0.095
fi0 = 0.01
fiF_basic = 0.5

def float_range(start, stop, step):
    curr = start
    while curr < stop:
        yield curr
        curr += step

with open("csvaam.csv", "r") as inputaam:
    reader = csv.reader(inputaam)
    am_all = list(reader)

a, m = [], []

for a_or_m in am_all[0]:
    if am_all[0].index(a_or_m) % 2 == 1:
        m.append(float(a_or_m))
    else:
        a.append(float(a_or_m))

def VG_calc(F, tG):
    return F * tG

def B_calc(fiF, fi0, VG):
    return (fiF - fi0) / VG

def tM_calc(VM, F):
    return (VM / F)

def u_calc(L, tM, ):
    return (L * 1000 /(tM * 60))

def H_calc(A, B, C, u):
    return (A + (B / u) + C * u)

def N_calc(L, H):
    return(L * 1000 / H)

def ke_calc(fi0, VM, B):
    ke_pom = []
    for j in range(len(a)):
        ke_pom.append(1 / ((2.31 * m[j] * B * VM) + (10 ** (m[j] * fi0 - a[j]))))
    return ke_pom

def wg_calc(VM, N, ke_l):
    wg_sum = 0
    for j in range(len(a) - 1):
        wg = ((4 * VM * (1 + ke_l[j + 1])) / N ** 0.5)
        wg_sum += wg
        wg_avr = wg_sum / (len(a) - 1)
    return wg_avr

def VR_calc(fi0, VM, VD, B):
    return ((1 / (m[-1] * B)) * np.log10(2.31 * m[-1] * B * (VM * 10 ** (a[-1] - (m[-1] * fi0)) - VD) + 1) + VM + VD)

def Fd_calc(F):
    return F * 0.0000000000017

def pressure_calc(F, Kf, L, ID, vis):
    return (F * vis * L / (Kf * np.pi * ((ID / 1000) / 2) ** 2)) / 1000000

def nc_tG_calc(wgavr, VG):
    return (VG / wgavr)

def fit_F_dependency(tG):
    F_all = []
    for i in float_range(0.00006, 0.0008, 0.000002):
        F_all.append(i)
    nc = []
    for i in range(len(F_all)):
        VG_need = VG_calc(F_all[i], tG)
        B_need = B_calc(fiF_basic,fi0,VG_need)
        tM_need = tM_calc(VM, F_all[i])
        u_need = u_calc(L, tM_need)
        H_need = H_calc(A, B, C, u_need)
        N_need = N_calc(L, H_need)
        ke_all = ke_calc(fi0, VM, B_need)
        VR_all = []
        wg_sum = 0
        for j in range(len(a) - 1):
            VR_all.append(( 1 / (m[j + 1] * B_need)) * np.log10(2.31 * m[j + 1] * B_need * (VM * 10 ** (a[j + 1] - (m[j + 1] * fi0)) - VD) + 1) + VM + VD)
            wg = ((4 * VM * (1 + ke_all[j + 1])) / (N_need) ** 0.5)
            wg_sum += wg
        wg_avr = wg_sum / len(a)
        nc.append(VG_need / wg_avr)
    for k, item in enumerate(F_all):
        F_all[k] = item * 1000
    plt.plot(F_all, nc)

def iteration(tG):
    nc_max = 0
    F_max = 0
    fiF_max = 0
    for F in float_range(0.00005, 0.00015, 0.000001):
        for fiF in float_range(0.1, 0.5, 0.01):
            nc_next = nc_tG_calc(wg_calc(VM,N_calc(L,H_calc(A,B,C,u_calc(L,tM_calc(VM,F)))),ke_calc(fi0,VM,B_calc(fiF,fi0, VG_calc(F,tG)))), VG_calc(F, tG))
            VG = VG_calc(F, tG)
            VRn = VR_calc(fi0, VM, VD, B_calc(fiF, fi0, VG_calc(F, tG)))
            p = pressure_calc(Fd_calc(F), Kf, L, ID, vis)
            if nc_next > nc_max:
                if VRn <= VG:
                    if p <= 25:
                        nc_max = nc_next
                        F_max = F
                        fiF_max = fiF
                        print(fiF)
                        print(VRn)
                        print(VG)
                        print(p)
                        print(F)
    return nc_max, F_max, fiF_max

def iteration_plotting():
    plt.plot(tG_all, nC_max_all)
    plt.xlabel("tG (min)")
    plt.ylabel("nC")
    plt.show()

def plotting_f_dependency():
    fit_F_dependency(120)
    fit_F_dependency(90)
    fit_F_dependency(60)
    plt.xlabel("F (ul $min^{-1}$)")
    plt.ylabel("nC")
    plt.legend(["120 min", "90 min", "60 min"], loc="upper right")
    plt.show()

plotting_f_dependency()

#iteration conditions
nC_max_all = []
tG_all = []
for tG in range(30, 500, 30):
    nC_max_all.append(iteration(tG)[0])
    tG_all.append(tG)

iteration_plotting()















