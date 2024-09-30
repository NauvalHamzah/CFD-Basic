import copy
import numpy as np
import os
import csv

import cores.CoreMethod as CoreMethod
import cores.PlotMethod as PlotMethod

#NMax 201 non isentropic
X201=PlotMethod.CSVdat('CSV/X201.csv')
A201=PlotMethod.CSVdat('CSV/A201.csv')
rho201_2=PlotMethod.CSVdat('CSV/2rho201.csv')
T201_2=PlotMethod.CSVdat('CSV/2T201.csv')
V201_2=PlotMethod.CSVdat('CSV/2V201.csv')
P201_2=PlotMethod.CSVdat('CSV/2P201.csv')
M201_2=PlotMethod.CSVdat('CSV/2M201.csv')
mdot201_2=PlotMethod.CSVdat('CSV/2mdot201.csv')

#NMax 201 Isentropic
rho201_1=PlotMethod.CSVdat('CSV/1rho201.csv')
T201_1=PlotMethod.CSVdat('CSV/1T201.csv')
V201_1=PlotMethod.CSVdat('CSV/1V201.csv')
P201_1=PlotMethod.CSVdat('CSV/1P201.csv')
M201_1=PlotMethod.CSVdat('CSV/1M201.csv')
mdot201_1=PlotMethod.CSVdat('CSV/1mdot201.csv')

#NMax 21 non isentropic
X21=PlotMethod.CSVdat('CSV/X21.csv')
A21=PlotMethod.CSVdat('CSV/A21.csv')
rho21_2=PlotMethod.CSVdat('CSV/2rho21.csv')
T21_2=PlotMethod.CSVdat('CSV/2T21.csv')
V21_2=PlotMethod.CSVdat('CSV/2V21.csv')
P21_2=PlotMethod.CSVdat('CSV/2P21.csv')
M21_2=PlotMethod.CSVdat('CSV/2M21.csv')
mdot21_2=PlotMethod.CSVdat('CSV/2mdot21.csv')

#NMax 21 Isentropic
rho21_1=PlotMethod.CSVdat('CSV/1rho21.csv')
T21_1=PlotMethod.CSVdat('CSV/1T21.csv')
V21_1=PlotMethod.CSVdat('CSV/1V21.csv')
P21_1=PlotMethod.CSVdat('CSV/1P21.csv')
M21_1=PlotMethod.CSVdat('CSV/1M21.csv')
mdot21_1=PlotMethod.CSVdat('CSV/1mdot21.csv')

#NMax 21 Isentropic Analytical
rho21_0=PlotMethod.CSVdat('CSV/0rho21.csv')
T21_0=PlotMethod.CSVdat('CSV/0T21.csv')
V21_0=PlotMethod.CSVdat('CSV/0V21.csv')
P21_0=PlotMethod.CSVdat('CSV/0P21.csv')
M21_0=PlotMethod.CSVdat('CSV/0M21.csv')
mdot21_0=PlotMethod.CSVdat('CSV/0mdot21.csv')

PlotMethod.Plotdata3(X21, X21, X201, rho21_0, rho21_1, rho201_1, 2, 'Density Isentropic', r'$\rho/\rho_{0}$')
PlotMethod.Plotdata3(X21, X21, X201, T21_0, T21_1, T201_1, 2, 'Temperature Isentropic', r'$T/T_{0}$')
PlotMethod.Plotdata3(X21, X21, X201, P21_0, P21_1, P201_1, 2, 'Pressure Isentropic', r'$P/P_{0}$')
PlotMethod.Plotdata3(X21, X21, X201, V21_0, V21_1, V201_1, 2, 'Non-Dimensional Velocity Isentropic', r'$V/a_{0}$')
PlotMethod.Plotdata3(X21, X21, X201, M21_0, M21_1, M201_1, 2, 'Mach Number Isentropic', r'$M$')
PlotMethod.Plotdata3(X21, X21, X201, mdot21_0, mdot21_1, mdot201_1, 2, 'Non-Dimensional Mass Flow Rate Isentropic', r'$(\rho A V)/\rho_{0} A^{*} a_{0}$')

PlotMethod.Plotdata2(X21, X201, rho21_2, rho201_2, 2, 'Density Non-Isentropic', r'$\rho/\rho_{0}$')
PlotMethod.Plotdata2(X21, X201, T21_2, T201_2, 2, 'Temperature Non-Isentropic', r'$T/T_{0}$')
PlotMethod.Plotdata2(X21, X201, P21_2, P201_2, 2, 'Pressure Non-Isentropic', r'$P/P_{0}$')
PlotMethod.Plotdata2(X21, X201, V21_2, V201_2, 2, 'Non-Dimensional Velocity Non-Isentropic', r'$V/a_{0}$')
PlotMethod.Plotdata2(X21, X201, M21_2, M201_2, 2, 'Mach Number Non-Isentropic', r'$M$')
PlotMethod.Plotdata3(X21, X21, X201, mdot21_0, mdot21_2, mdot201_2, 2, 'Non-Dimensional Mass Flow Rate Non-Isentropic', r'$(\rho A V)/\rho_{0} A^{*} a_{0}$')




