import pandas as pd
import time
import rpy2.robjects as robjects
from rpy2.robjects import IntVector
from rpy2.robjects.packages import importr
import os


os.environ['R_HOME'] = r'C:\Program Files\R\R-4.5.1' # Set local R path

simARG = importr("simARG")
n = 20
rho_site = 10 / 1e5
delta = 30
R_time_vec = []
Rcpp_time_vec = []
clean_time_vec = []
for i in range(1, 101):
    robjects.r(f'set.seed({i})')

    start_time_R = time.time()
    result = simARG.simbac_ARG(n, rho_site, int(1e5), delta, node_max = 1000)
    end_time_R = time.time()
    R_time_vec.append(end_time_R - start_time_R)

    start_time_Rcpp = time.time()
    result = simARG.simbac_ARG_decimal(n, rho_site, int(1e5), delta, node_max = 1000)
    end_time_Rcpp = time.time()
    Rcpp_time_vec.append(end_time_Rcpp - start_time_Rcpp)

    start_time_clean = time.time()
    robjects.r('gc()')
    end_time_clean = time.time()
    clean_time_vec.append(end_time_clean - start_time_clean)

    print(f'Complete {i} iterations')


print(R_time_vec)
print(Rcpp_time_vec)
