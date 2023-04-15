# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 19:09:06 2022

@author: DiDoonka
"""
import numpy as np
import csv
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


# print("Write column length in meters: ")
# L = float(input())
# print("Write inner diameter column in milimeters: ")
# ID = float(input())
# print("Write viscosity of mobile phase in pa.s: ")
# VIS = float(input())

Dm = 0.00000000017246
L = 189
ID = 0.075
VIS = 0.95

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
    splitted_line = line.split(";")
    return splitted_line

def counts():
    count_peptides = (data[1][0].count(";") - 6) // 2
    count_gradients = len(data) - 1
    return count_gradients, count_peptides

def calc_VG():
    for i in range(1, count_gradients + 1):
        splitted_line = split_line(i)
        fiA.append(float(splitted_line[0]))
        VG_calc = float(splitted_line[3]) * float(splitted_line[2])
        add_excel(str(VG_calc), i)

def calc_B():
    for i in range(1, count_gradients + 1):
        splitted_line = split_line(i)
        B_calc = (float(splitted_line[1]) - float(splitted_line[0])) / float(splitted_line[count_peptides * 3])
        B.append(B_calc)
        add_excel(str(round(B_calc, 3)), i)

def calc_VD():
    VD_avr = 0
    for i in range(1, count_gradients + 1):
        splitted_line = split_line(i)
        VD_calc = float(splitted_line[6]) * float(splitted_line[3])
        VD_avr += VD_calc
    VD_avr = VD_avr / count_gradients
    return VD_avr

def calc_VM():
    VM_avr = 0
    for i in range(1, count_gradients + 1):
        splitted_line = split_line(i)
        VM_calc = float(splitted_line[5]) * float(splitted_line[3])
        VM_avr += VM_calc
    VM_avr = VM_avr / count_gradients
    return VM_avr

def VD_VM_average(VD_avr, VM_avr):
    for i in range(1, count_gradients + 1):
        VD.append(VD_avr)
        VM.append(VM_avr)
        add_excel(str(round(VD_avr, 6)), i)
        add_excel(str(round(VM_avr, 6)), i)

def calc_VRexp_form_ydata():
    VRs_for_curve = []
    for i in range(1, count_gradients + 1):
        VR_for_curve = []
        splitted_line = split_line(i)
        for j in range(count_peptides, count_peptides * 2, 1):
            VR = float(splitted_line[j]) * float(splitted_line[3])
            add_excel(str(round(VR, 6)), i)
            VR_for_curve.append(float(round(VR, 6)))
        VRs_for_curve.append(VR_for_curve)
    y_value = []
    for i in range(count_peptides):
        for j in range(count_gradients):
            y_value.append(VRs_for_curve[j][i])
        y_values.append(y_value)
        y_value = []
    return y_values

def calc_F_Kf():
    Kf_avr = 0
    for i in range(1, count_gradients + 1):
        splitted_line = split_line(i)
        F = (float(splitted_line[3])) * float(Dm)
        add_excel(str(F), i)
        splitted_line = split_line(i)
        Kf = (float(splitted_line[11 + 3 * count_peptides]) * VIS * L) / (float(splitted_line[4]) * 1000000 * 3.14159 * ((ID / 2) ** 2))
        Kf_avr += Kf
    Kf_avr = Kf_avr / count_gradients
    for j in range(1, count_gradients + 1):
        add_excel(str(Kf_avr), j)

def calc_N():
    for i in range(1, count_gradients + 1):
        splitted_line = split_line(i)
        N = 5.545 * ((float(splitted_line[7]) / float(splitted_line[count_peptides + 7])) ** 2)
        print(float(splitted_line[7]))
        print(float(splitted_line[count_peptides + 7]))
        add_excel(str(round(N, 0)), i)

def calc_u():
    all_U = []
    for i in range(1, count_gradients + 1):
        splitted_line = split_line(i)
        u = L / (float(splitted_line[5]) * 60)
        add_excel(str(round(u, 3)), i)
        all_U.append(float(u))
    return all_U

def calc_H():
    all_H = []
    for i in range(1, count_gradients + 1):
        splitted_line = split_line(i)
        H = L /  float(splitted_line[13 + 3 * count_peptides])
        add_excel(str(round(H, 3)), i)
        all_H.append(float(H))
    return all_H

def fit_am(X, a, m):
    VRc = []
    for i in range(len(X[0])):
        fiA, B, VM, VD = X[0][i], X[1][i], X[2][i], X[3][i]
        cond_zero = float(m) * float(B)
        cond_logarit = (2.31 * float(m) * float(B) *
                  (float(VM) * 10 ** (float(a) - float(m) * float(fiA)) - float(VD))
                        + 1)
        if cond_zero == 0 or cond_logarit <= 0 or cond_logarit >= 1e+150:
            VRc.append(float("nan"))
            continue
        VR = (1 / cond_zero) * np.log10(cond_logarit) + float(VM) + float(VD)
        VRc.append(VR)
    return np.array(VRc)

def find_a_m_RSS_RSE():
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
    a_save = "a;"
    m_save = "m;"
    for i in range(count_peptides):
        a_save += str(round(float(am[2 * i]), 4))
        m_save += str(round(float(am[2 * i + 1]), 4))
        if i != count_peptides - 1:
            a_save += ";"
            m_save += ";"
    data.append([a_save])
    data.append([m_save])
    data.append([RSEs])
    data.append(["RSS;" + str(round(float(RSS), 9))])

def write_data():
    for i in range(1, len(data)):
        writer.writerow([data[i][0]])

def fit_Hu(Hu_u, Hu_A, Hu_B, Hu_C):
    Hu_H = Hu_A + (Hu_B / Hu_u) + (Hu_C * Hu_u)
    return np.array(Hu_H)

def am_save():
    with open("csvaam.csv", "w") as write_am:
        writer = csv.writer(write_am, lineterminator='\n')
        writer.writerow(am)

def Hu_from_first(par_1, par_2, par_3):
    calc_N()
    all_U = calc_u()
    all_H = calc_H()
    popt2, pcov = curve_fit(fit_Hu, all_U, all_H)

    scope_H,  scope_U = [], []
    for U in float_range(par_1, par_2, par_3):
        scope_H.append(popt2[0] + (popt2[1] / U) + (popt2[2] * U))
        scope_U.append(U)

    plt.plot(scope_U, scope_H)
    plt.scatter(all_U, all_H)
    plt.xlabel("u (mm $s^{-1}$)")
    plt.ylabel("H (mm)")
    write_ABC = ("A;" + str(round(popt2[0], 5)) + ";B;" + str(round(popt2[1], 5)) + ";C;" + str(round(popt2[2], 5)))
    data.append([write_ABC])


if __name__ == "__main__":
    with open("csvdata.csv", "r") as input_file:
        reader = csv.reader(input_file)
        data = list(reader)

    with open("csvdata.csv", "w") as output_file:
        writer = csv.writer(output_file, lineterminator='\n')
        count_gradients, count_peptides = counts()
        fiA, B, VM, VD = [], [], [], []
        calc_VG()
        calc_B()
        VD_VM_average(calc_VD(), calc_VM())
        y_values = []
        y_values = calc_VRexp_form_ydata()
        writer.writerow([data[0][0]])
        fiA, B, VM, VD, y_values = np.array(fiA), np.array(B), np.array(VM), np.array(VD), np.array(y_values)
        am = []
        VRc_all = []
        for i in range(count_peptides):
            y_data = y_values[i]
            popt, pcov = curve_fit(fit_am, (fiA, B, VM, VD), y_data, bounds=(0, np.inf))
            print("a = " + str(popt[0]) + "\nm = " + str(popt[1]))
            am.append(popt[0])
            am.append(popt[1])
            a = (popt[0])
            m = (popt[1])
            VRc_all.append(fit_am((fiA, B, VM, VD), a, m).tolist())

        find_a_m_RSS_RSE()
        calc_F_Kf()

        Hu_from_first(0.05, 2.55, 0.05)

        am_save()
        write_data()
