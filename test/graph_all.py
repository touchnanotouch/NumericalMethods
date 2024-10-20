import re

import numpy as np
import random as rd
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd


# n_s = [10, 25, 50, 75, 100, 250]
# n_s = [10, 25, 50, 75, 100, 250, 500]
n_s = [10, 25, 50, 75, 100, 250, 500, 750, 1000]

infos = {}

for n in n_s:
    times = {}
    accuracies = {}

    with open(f"info_all_n{n}.txt", "r") as f:
        lines = f.readlines()

    found_time = False
    fount_acc = False

    for line in lines:
        if found_time:
            if "gauss" in line:
                times["gauss"] = (float(line.split(":")[1].strip()[:-1]))
            elif "thomas:" in line:
                times["thomas"] = (float(line.split(":")[1].strip()[:-1]))
            elif "lu" in line:
                times["lu"] = (float(line.split(":")[1].strip()[:-1]))
            elif "qroots" in line:
                times["qroots"] = (float(line.split(":")[1].strip()[:-1]))
                found_time = False
            else:
                if "si:" in line:
                    times["si"] = (float(line.split(":")[1].strip()[:-1]))
                elif "seidel" in line:
                    times["seidel"] = (float(line.split(":")[1].strip()[:-1]))
                elif "sor" in line:
                    times["sor"] = (float(line.split(":")[1].strip()[:-1]))
                elif "res" in line:
                    times["res"] = (float(line.split(":")[1].strip()[:-1]))
                elif "grad" in line:
                    times["grad"] = (float(line.split(":")[1].strip()[:-1]))
                    found_time = False
        elif fount_acc:
            if "gauss" in line:
                accuracies["gauss"] = (float(line.split(":")[1].strip()[:-1]))
            elif "thomas:" in line:
                accuracies["thomas"] = (float(line.split(":")[1].strip()[:-1]))
            elif "lu" in line:
                accuracies["lu"] = (float(line.split(":")[1].strip()[:-1]))
            elif "qroots" in line:
                accuracies["qroots"] = (float(line.split(":")[1].strip()[:-1]))
                fount_acc = False
            else:
                if "si:" in line:
                    accuracies["si"] = (float(line.split(":")[1].strip()[:-1]))
                elif "seidel" in line:
                    accuracies["seidel"] = (float(line.split(":")[1].strip()[:-1]))
                elif "sor" in line:
                    accuracies["sor"] = (float(line.split(":")[1].strip()[:-1]))
                elif "res" in line:
                    accuracies["res"] = (float(line.split(":")[1].strip()[:-1]))
                elif "grad" in line:
                    accuracies["grad"] = (float(line.split(":")[1].strip()[:-1]))
                    fount_acc = False

        if line.startswith("Iterative methods (execution time)") or line.startswith("Exact methods (execution time)"):
            found_time = True
            continue
        elif line.startswith("Iterative methods (accuracy)") or line.startswith("Exact methods (accuracy)"):
            fount_acc = True
            continue

    infos[n] = {"times": times, "accuracies": accuracies}

times_data = []
accuracies_data = []

for n in n_s:
    for algo in ["gauss", "thomas", "lu", "qroots", "si", "seidel", "sor", "res", "grad"]:
        times_data.append({"n": n, "algo": algo, "time": infos[n]["times"][algo]})
        accuracies_data.append({"n": n, "algo": algo, "accuracy": infos[n]["accuracies"][algo]})

times_df = pd.DataFrame(times_data)
accuracies_df = pd.DataFrame(accuracies_data)

times_df.to_csv("times.csv")
accuracies_df.to_csv("accuracies.csv")

fig = px.line(times_df, x="n", y="time", color="algo", markers=True)

fig.update_layout(
    title="График времени выполнения всех реализованных методов",
    legend=dict(title="Методы"),
    xaxis_title="Количество переменных в СЛАУ (n)",
    yaxis_title="Время выполнения (s)",
    yaxis=dict(tickformat=".2ef")
)

fig1 = px.line(accuracies_df, x="n", y="accuracy", color="algo", markers=True)

fig1.update_layout(
    title="График точности решений СЛАУ, полученных всеми реализованными методами",
    legend=dict(title="Методы"),
    xaxis_title="Количество переменных в СЛАУ (n)",
    yaxis_title="Точность решения (%)",
    yaxis=dict(tickformat=".2ef")
)

fig.show()
fig1.show()
