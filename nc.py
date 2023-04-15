# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 18:57:36 2023

@author: DiDoonka
"""
import csv
import numpy as np
import matplotlib.pyplot as plt

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

fi0 = 0.01
fiF = 0.5
VM = 0.001061
VD = 0.002785

A = -5.4085
B = 2.2476
C = 3.633402
KF = 0.0000000000000000713
L = 0.192


def fit(tG):
    F_all = []
    for i in float_range(0.0001, 0.0008, 0.000002):
        F_all.append(i)

    VG_all, B_all, tM_all, u_all, H_all, N_all, ke_all, VR_all, wg_all, nc = [],[],[],[],[],[],[],[],[],[]
    for i in range(len(F_all)):
        VG_all.append(F_all[i]*tG)
        B_all.append((fiF-fi0)/VG_all[i])
        tM_all.append(VM/F_all[i])
        u_all.append(L * 1000/(tM_all[i] * 60))
        H_all.append(A+(B/u_all[i])+C*u_all[i])
        N_all.append(L * 1000/H_all[i])
        ke_pom = []
        for j in range(len(a)):
            ke_pom.append(1/((2.31*m[j]*B_all[i]*VM)+(10**(m[j]*fi0-a[j]))))
        ke_all.append(ke_pom)
        VR_pom = []
        for j in range(len(a)):
            VR_pom.append((1/(m[j]*B_all[i]))*np.log10(2.31 * m[j] * B_all[i] * (VM * 10**(a[j] - (m[j] * fi0)) - VD) + 1) + VM + VD)
        VR_all.append(VR_pom)
        wg_sum = 0
        for j in range(len(a)):
            wg = ((4 * VM * (1 + ke_all[i][j]))/(N_all[i])**0.5)
            wg_sum += wg
        wg_avr = wg_sum / len(a)
        wg_all.append(wg_avr)
        nc.append((VR_all[i][-1] - VR_all[i][0])/wg_all[i])

    for i, item in enumerate(F_all):
        F_all[i] = item * 1000

    plt.plot(F_all, nc)

fit(120)
fit(90)
fit(60)
plt.xlabel("F (ul $min^{-1}$)")
plt.ylabel("nC")
plt.legend(["120 min", "90 min", "60 min"], loc="upper right")
plt.show

















