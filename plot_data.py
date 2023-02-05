import csv
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from math import ceil

def plot_objective(folder):
    file1 = folder+"Branch_and_cut_sol.csv"    
    file2 = folder+"dualisation_sol.csv"
    file3 = folder+"heuristique_sol.csv"    
    file4 = folder+"plans_coupants_sol.csv"
    file5 = folder+"statique_sol.csv"
    statique = dict()
    relaxation = dict()
    dualisation = dict()
    BB = dict()
    heuristique = dict()
    PC = dict()
    csv_reader = csv.reader(open(file5), delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        elif line_count == 1:
            line_count += 1
        else:
            statique[row[0]]= float(row[2])
            line_count += 1

    csv_reader = csv.reader(open(file1), delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        elif line_count == 1:
            line_count += 1
        else:
            BB[row[0]]= float(row[2])
            relaxation[row[0]] = max(relaxation.get(row[0], 0), float(row[4]))
            line_count += 1
    
    csv_reader = csv.reader(open(file2), delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        elif line_count == 1:
            line_count += 1
        else:
            dualisation[row[0]]= float(row[2])
            relaxation[row[0]] = max(relaxation.get(row[0], 0), float(row[4]))

            line_count += 1

    
    csv_reader = csv.reader(open(file3), delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        elif line_count == 1:
            line_count += 1
        else:
            val = row[2].split(" / ")[1]
            heuristique[row[0]]= float(val)
            line_count += 1

    
    csv_reader = csv.reader(open(file4), delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        elif line_count == 1:
            line_count += 1
        else:
            PC[row[0]]= float(row[2])
            relaxation[row[0]] = max(relaxation.get(row[0], 0), float(row[4]))
            line_count += 1


    fig, ax = plt.subplots(figsize=(15,7))
    plt.yscale('log')
    for i, val in enumerate(statique.values()):
        if i ==0:
            plt.hlines(val, i-0.3, i+0.3, "r",label="statique")
        else:
            plt.hlines(val, i-0.3, i+0.3, "r")
    for i, val in enumerate(relaxation.values()):
        if i ==0:
            plt.hlines(val, i-0.3, i+0.3, "g",label="relaxation")
        else:
            plt.hlines(val, i-0.3, i+0.3, "g") 
    plt.plot(heuristique.keys(), [heuristique.get(key, 0) for key in heuristique.keys()],"o",color=mcolors.TABLEAU_COLORS["tab:orange"], label="heuristique")
    plt.plot(dualisation.keys(), [dualisation.get(key, 0) for key in dualisation.keys()], "o", color=mcolors.TABLEAU_COLORS["tab:purple"], label="dualisation")    
    plt.plot(BB.keys(), [BB.get(key, 0) for key in BB.keys()],"o", color = mcolors.TABLEAU_COLORS["tab:cyan"], label="Branch&Bound")
    plt.plot(PC.keys(), [PC.get(key,0) for key in PC.keys()], "o", color = mcolors.TABLEAU_COLORS["tab:olive"], label="plans coupants")   
    plt.xticks(rotation=90)
    plt.legend()
    fig.tight_layout()
    plt.show()



def plot_nb_instances(folder, gap_max):
    gap_max = gap_max/100
    file1 = folder+"Branch_and_cut_sol.csv"    
    file2 = folder+"dualisation_sol.csv"
    file3 = folder+"heuristique_sol.csv"    
    file4 = folder+"plans_coupants_sol.csv"

    time_max = 0
    relaxation = dict()
    dualisation = dict()
    BB = dict()
    heuristique = dict()
    PC = dict()

    csv_reader = csv.reader(open(file1), delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        elif line_count == 1:
            line_count += 1
        else:
            BB[row[0]]= (float(row[2]), float(row[1]))
            time_max = max(time_max, float(row[1]))
            relaxation[row[0]] = max(relaxation.get(row[0], 0), float(row[4]))
            line_count += 1
    
    csv_reader = csv.reader(open(file2), delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        elif line_count == 1:
            line_count += 1
        else:
            dualisation[row[0]]= (float(row[2]), float(row[1]))
            time_max = max(time_max, float(row[1]))
            relaxation[row[0]] = max(relaxation.get(row[0], 0), float(row[4]))
            line_count += 1

    
    csv_reader = csv.reader(open(file3), delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        elif line_count == 1:
            line_count += 1
        else:
            val = row[2].split(" / ")[1]
            heuristique[row[0]]= (float(val), float(row[1]))
            relaxation[row[0]] = relaxation.get(row[0], 0)
            time_max = max(time_max, float(row[1]))
            line_count += 1

    
    csv_reader = csv.reader(open(file4), delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        elif line_count == 1:
            line_count += 1
        else:
            PC[row[0]]= (float(row[2]), float(row[1]))
            relaxation[row[0]] = max(relaxation.get(row[0], 0), float(row[4]))
            time_max = max(time_max, float(row[1]))
            line_count += 1


    fig, ax = plt.subplots(figsize=(7,7))
    time = np.linspace(0, ceil(time_max)+3, ceil(time_max)+4)
    nb_dual = [0 for t in time]
    nb_BB = [0 for t in time]
    nb_PC = [0 for t in time]
    nb_heur = [0 for t in time]

    for key,elem in dualisation.items():
        if 1-relaxation[key]/elem[0] <= gap_max:
            t = upper_tick_from_time(time, elem[1])
            for i in range(len(time)):
                if i >= t:
                    nb_dual[i] +=1
        else:
            print(key, elem, relaxation[key], 1-relaxation[key]/elem[0])
    for key,elem in heuristique.items():
        if 1-relaxation[key]/elem[0] <= gap_max:
            t = upper_tick_from_time(time, elem[1])
            for i in range(len(time)):
                if i >= t:
                    nb_heur[i] +=1
    for key,elem in BB.items():
        if 1-relaxation[key]/elem[0] <= gap_max:
            t = upper_tick_from_time(time, elem[1])
            print(elem[1], t)
            for i in range(len(time)):
                if i >= t:
                    nb_BB[i] +=1
    for key,elem in PC.items():
        if 1-relaxation[key]/elem[0] <= gap_max:
            t = upper_tick_from_time(time, elem[1])
            for i in range(len(time)):
                if i >= t:
                    nb_PC[i] +=1

    plt.plot(time, nb_dual, color=mcolors.TABLEAU_COLORS["tab:purple"], label="dualisation")    
    plt.plot(time, nb_BB, color = mcolors.TABLEAU_COLORS["tab:cyan"], label="Branch&Bound")
    plt.plot(time, nb_heur,color=mcolors.TABLEAU_COLORS["tab:orange"], label="heuristique")
    plt.plot(time, nb_PC, color = mcolors.TABLEAU_COLORS["tab:olive"], label="plans coupants")   
    plt.legend()
    plt.xlabel("Temps (s)")
    plt.ylabel("Nombre d'instances résolues")
    fig.tight_layout()
    plt.show()

def upper_tick_from_time(time_table, time):
    a = 0
    b = len(time_table)-1
    while b-a>1:
        middle = (b+a)//2
        if time_table[middle] < time:
            a = middle
        else:
            b = middle
    return b

plot_nb_instances("30 s all/", 100)
#plot_nb_instances("5min/", 50)