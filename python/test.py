import pandas as pd
import time
import rpy2.robjects as robjects
from rpy2.robjects import IntVector
from rpy2.robjects.packages import importr
import os


os.environ['R_HOME'] = r'C:\Program Files\R\R-4.5.1' # Set local R path

robjects.r("library(simARG)")

n = IntVector([20])
robjects.globalenv['n'] = n
rho_site = 10 / 1e5
delta = 30
L = IntVector([1e5])
robjects.globalenv['L'] = L
R_time_vec = []
Rcpp_time_vec = []
clean_time_vec = []
for i in range(1, 101):
    robjects.r(f'set.seed({i})')

    start_time_R = time.time()
    robjects.r(f'result <- simbac_ARG(n, {rho_site}, L, {delta}, node_max = 1000)')
    end_time_R = time.time()
    R_time_vec.append(end_time_R - start_time_R)

    start_time_Rcpp = time.time()
    robjects.r(f'result <- simbac_ARG.decimal(n, {rho_site}, L, {delta}, node_max = 1000)')
    end_time_Rcpp = time.time()
    Rcpp_time_vec.append(end_time_Rcpp - start_time_Rcpp)

    start_time_clean = time.time()
    robjects.r('gc()')
    end_time_clean = time.time()
    clean_time_vec.append(end_time_clean - start_time_clean)

    print(f'Complete {i} iterations')

print(R_time_vec)
print(Rcpp_time_vec)
