# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 19:09:06 2022

@author: DiDoonka
"""
import numpy as np
import csv
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# - bud vse cesky nebo vse anglicky
# - konkretnejsi vice vypovidajici nazvy
# - bud se drzet xxx_xxx nebo xxxXxx konvence, nepouzivat oboje - v pythonu dulezitejsi konvence xxx_xxx

# print("Write column length in meters: ")
# L = float(input())
# print("Write inner diameter column in meters: ")
# ID = float(input())
# print("Write viscosity of mobile phase in pa.s: ")
# VIS = float(input())

Dm = 0.00000000017246
L = 0.192
ID = 0.0001
VIS = 0.001

def float_range(start, stop, step):
    curr = start
    while curr < stop:
        yield curr
        curr += step

def add_excel(variable, i):
    data[i][0] += ";"
    data[i][0] += variable

def split_line(i):
    line = data[i][0]
    x = line.split(";")
    return x

def counts():
    count_peptides = ((data[1][0]).count(";") - 6)//2
    count_gradients = len(data) - 1
    return count_gradients, count_peptides

def calc_VG():
    for i in range(1, len(data)):
        x = split_line(i)
        VG = (float(x[3])) * (float(x[2]))
        add_excel(str(VG), i)

def calc_B():
    for i in range(1, len(data)):
        x = split_line(i)
        B = (float(x[1]) - float(x[0])) / float(x[count_peptides*3])
        add_excel(str(round(B,3)), i)

def calc_VD():
    VDavr = 0
    for i in range(1, len(data)):
        x = split_line(i)
        VD = (float(x[6]) * float(x[3]))
        VDavr += VD
    VDavr = VDavr / count_gradients
    return VDavr

def calc_VM():
    VMavr = 0
    for i in range(1, len(data)):
        x = split_line(i)
        VM = (float(x[5]) * float(x[3]))
        VMavr += VM
    VMavr = VMavr / count_gradients
    return VMavr

def VD_VM_average(VDavr, VMavr):
    for i in range(1, len(data)):
        add_excel(str(round(VDavr, 6)), i)
        add_excel(str(round(VMavr, 6)), i)

def calc_VRexp():
    VRs_for_curve = []
    for i in range(1, len(data)):
        VR_for_curve = []
        x = split_line(i)
        for j in range(count_peptides, count_peptides * 2, 1):
            VR = float(x[j]) * float(x[3])
            add_excel(str(round(VR, 6)), i)
            VR_for_curve.append(float(round(VR, 6)))
        VRs_for_curve.append(VR_for_curve)
    y_value = []
    for l in range(count_peptides):
        for k in range(count_gradients):
            y_value.append(VRs_for_curve[k][l])
        y_values.append(y_value)
        y_value = []
    return y_values

def calc_F():
    for i in range(1, count_gradients + 1):
        x = split_line(i)
        F = (float(x[3])) * float(Dm)
        add_excel(str(F), i)

def calc_Kf():
    for i in range(1, count_gradients + 1):
        x = split_line(i)
        Kf = ((float(x[11 + 3 * count_peptides]) * VIS * L) / (float(x[4]) * 1000000 * 3.14159 * ((ID / 2) ** 2)))
        add_excel(str(Kf), i)

def calc_N():
    for i in range(1, count_gradients + 1):
        x = split_line(i)
        N = ((5.545 * (float(x[7]) / float(x[count_peptides + 7])) ** 2))
        add_excel(str(round(N, 0)), i)

def calc_u():
    all_U = []
    for i in range(1, count_gradients + 1):
        x = split_line(i)
        u = (L * 1000 / (float(x[5]) * 60))
        add_excel(str(round(u, 3)), i)
        all_U.append(float(u))
    return all_U

def calc_H():
    all_H = []
    for i in range(1, count_gradients + 1):
        x = split_line(i)
        H = (L * 1000 /  float(x[13 + 3 * count_peptides]))
        add_excel(str(round(H, 3)), i)
        all_H.append(float(H))
    return all_H

def fit_am(X, a, m):
    VRc = []
    for i in range(len(X[0])):
        fiA, B, VM, VD = X[0][i], X[1][i], X[2][i], X[3][i]
        null = (float(m) * float(B))
        logarit = (2.31 * float(m) * float(B) *
                  (float(VM) * 10 ** (float(a) - float(m) * float(fiA)) - float(VD)) + 1)
        if null == 0 or logarit <= 0:
            VRc.append(float("nan"))
            continue
        VR = (1 / null) * np.log10(logarit) + float(VM) + float(VD)
        VRc.append(VR)
    return np.array(VRc)

def find_a_m_RSS_RSE(): # vypise a m RSS RSE
    RSEs = "RSE;"
    RSS = 0
    for i in range(count_peptides):
        RSE = 0
        for j in range(count_gradients):
            RSE += (VRc_all[i][j] - y_values[i][j]) ** 2
        print("RSE = " + str(RSE))
        RSEs += str(round(float(RSE), 9))
        if i != count_peptides - 1:
            RSEs += ";"
        RSS += RSE
    print("RSS = " + str(RSS))
    a = "a;"
    m = "m;"
    for i in range(count_peptides):
        a += str(round(float(am[(2 * i)]), 4))
        if i != count_peptides - 1:
            a += ";"
        m += str(round(float(am[(2 * i + 1)]), 4))
        if i != count_peptides - 1:
            m += ";"
    data.append([a])
    data.append([m])
    data.append([RSEs])
    data.append(["RSS;" + str(round(float(RSS), 9))])

def writing():
    for i in range(1, len(data)):
        writer.writerow([data[i][0]])

def fitovani_Hu(Hu_u, A, B, C):
    Hu_H = A + (B / Hu_u) + (C * Hu_u)
    return np.array(Hu_H)

def am_save():
    with open("csvaam.csv", "w") as write_am:
        writer = csv.writer(write_am, lineterminator='\n')
        writer.writerow(am)

with open("csvdata.csv", "r") as csvinput:
    reader = csv.reader(csvinput)
    data = list(reader)

with open("csvdata.csv", "w") as csvoutput:
    writer = csv.writer(csvoutput, lineterminator='\n')
    count_gradients, count_peptides = counts()
    calc_VG()
    calc_B()
    VD_VM_average(calc_VD(), calc_VM())
    y_values = []
    y_values = calc_VRexp()
    writer.writerow([data[0][0]])

    fiA, B, VM, VD = [], [], [], []

    for i in range(1, count_gradients + 1): # výběr dat pro fit
        x = split_line(i)
        fiA.append(float(x[0]))
        B.append(float(x[count_peptides * 2 + 8]))
        VM.append(float(x[count_peptides * 2 + 10]))
        VD.append(float(x[count_peptides * 2 + 9]))

    fiA, B, VM, VD, y_values = np.array(fiA), np.array(B), np.array(VM), np.array(VD), np.array(y_values)

    am = []
    for i in range(count_peptides): # fit
        y_data = y_values[i]
        popt, pcov = curve_fit(fit_am, (fiA, B, VM, VD), y_data, bounds = (0, np.inf))
        print("a = " + str(popt[0]) + "\nm = " + str(popt[1]))
        am.append(popt[0])
        am.append(popt[1])

    VRc_all = []
    for j in range(count_peptides):
        a = am[(2 * j)]
        m = am[(2 * j + 1)]
        VRc_all.append(fit_am((fiA, B, VM, VD), a, m).tolist())

    if len(data) == count_gradients + 1:
        find_a_m_RSS_RSE()
    calc_F()
    calc_Kf()
    calc_N()
    all_U = calc_u()
    all_H = calc_H()
    popt2,pcov = curve_fit(fitovani_Hu, all_U, all_H)

    scope_H,  scope_U = [], []
    for U in float_range(0.5, 1.55, 0.05):
        scope_H.append(popt2[0] + (popt2[1] / U) + (popt2[2] * U))
        scope_U.append(U)

    plt.plot(scope_U, scope_H)
    plt.scatter(all_U, all_H)
    plt.xlabel("u (mm $s^{-1}$)")
    plt.ylabel("H (mm)")
    write_ABC = ("A;" + str(round(popt2[0], 3)) + ";B;" + str(round(popt2[1], 3)) + ";C;" + str(round(popt2[2], 4)))
    data.append([write_ABC])

    am_save()
    writing()