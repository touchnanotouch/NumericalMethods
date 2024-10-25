import re

import numpy as np
import random as rd
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd


n_s = [10, 25, 50, 75, 100]
# n_s = [10, 25, 50, 75, 100, 250]
# n_s = [10, 25, 50, 75, 100, 250, 500]
# n_s = [10, 25, 50, 75, 100, 250, 500, 750, 1000]

infos = {}

for n in n_s:
    times = {}
    accuracies = {}

    with open(f"info_exact_n{n}.txt", "r") as f:
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

        if line.startswith("Exact methods (execution time)"):
            found_time = True
            continue
        elif line.startswith("Exact methods (accuracy)"):
            fount_acc = True
            continue

    infos[n] = {"times": times, "accuracies": accuracies}

times_data = []
accuracies_data = []

for n in n_s:
    for algo in ["thomas", "lu", "qroots"]:
        times_data.append({"n": n, "algo": algo, "time": infos[n]["times"][algo]})
        accuracies_data.append({"n": n, "algo": algo, "accuracy": infos[n]["accuracies"][algo]})

times_df = pd.DataFrame(times_data)
accuracies_df = pd.DataFrame(accuracies_data)

fig = px.line(times_df, x="n", y="time", color="algo", markers=True)

fig.update_layout(
    title="График времени выполнения итерационных методов",
    legend=dict(title="Методы"),
    xaxis_title="Количество переменных в СЛАУ (n)",
    yaxis_title="Время выполнения (s)",
    yaxis=dict(tickformat=".2ef")
)

fig1 = px.line(accuracies_df, x="n", y="accuracy", color="algo", markers=True)

fig1.update_layout(
    title="График точности решений СЛАУ, полученных итерационными методами",
    legend=dict(title="Методы"),
    xaxis_title="Количество переменных в СЛАУ (n)",
    yaxis_title="Точность решения (%)",
    yaxis=dict(tickformat=".2ef")
)

fig.write_html("exact_times.html")
fig1.write_html("exact_accuracies.html")

fig.show()
fig1.show()
