import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np 

#path="C:/Users/Miki/uni/CPPrakt/ISING_V4/ISING_EXP4/IsingVS/"
path="output_files/testing_therm_params/"
#determine optimal multihit parameter and thermal steps

def autocorrelation(x, max_lag=None):
    x = np.asarray(x)
    x = x - np.mean(x)
    n = len(x)

    if max_lag is None:
        max_lag = n // 10

    var = np.var(x)
    acf = np.empty(max_lag + 1)

    for lag in range(max_lag + 1):
        acf[lag] = np.mean(x[:n-lag] * x[lag:]) / var

    return acf

def tau_int_from_acf(acf):
    tau = 0.5
    for c in acf[1:]:
        if c <= 0:
            break
        tau += c
    return tau



def exercise3a():
    '''
    

    '''

    df = pd.read_table(path+"testing_metro_L=128_beta=0.440_MH=2_start=false.txt") 
    #df = pd.read_table(path+"testing_metro_L=128_beta=0.500_MH=2.txt") 
    e = np.array(df.iloc[:, 1])
    m = np.array(df.iloc[:, 3])
    it = np.arange(0,len(m), step = 1)

    fig, ax = plt.subplots()
    ax.scatter(it, m, marker="s", facecolors = "none", edgecolors="b")
    ax.set_xlabel("no. sweeps")
    ax.set_ylabel("magnetisation per particle")
    ax.set_title(r'$\beta = 0.6, N = 5$')
    ax.grid(True)
    plt.show()

    fig, ax = plt.subplots()
    ax.scatter(it, e, marker="s", facecolors = "none", edgecolors="b")
    ax.set_xlabel("no. sweeps")
    ax.set_ylabel("energy density")
    ax.set_title(r'$\beta = 0.6, N = 5$')
    ax.grid(True)
    plt.show()
    
def corr_plot():
    df = pd.read_table(path+"testing_metro_L=128_beta=0.440_MH=1_start=true.txt") 
    
    e = np.array(df.iloc[13000:14000, 1])
    abs_m = np.array(df.iloc[13000:14000, 3])
    it = np.arange(0,len(abs_m), step = 1)
    print(abs_m)
    fig, ax = plt.subplots()
    acf_abs_m = autocorrelation(abs_m, max_lag=300)
    #print(acf_e)
    #print(np.isnan(acf_e).any())
    #print(np.min(e), np.max(e))
    #print(np.var(abs_m))
    plt.plot(acf_abs_m, label="ACF |m|", color = "r")
    plt.grid()
    plt.show()
    
    df = pd.read_table(path+"testing_metro_L=128_beta=0.440_MH=1_start=false.txt") 
    
    e = np.array(df.iloc[13000:14000, 1])
    abs_m = np.array(df.iloc[13000:14000, 3])
    it = np.arange(0,len(abs_m), step = 1)
    print(abs_m)
    fig, ax = plt.subplots()
    acf_e = autocorrelation(e, max_lag=300)
    #print(acf_e)
    #print(np.isnan(acf_e).any())
    #print(np.min(e), np.max(e))
    #print(np.var(abs_m))
    plt.plot(acf_e, label="ACF |m|", color = "r")
    plt.grid()
    plt.show()
    tau_e = tau_int_from_acf(acf_e)
    tau_abs_m = tau_int_from_acf(acf_abs_m)

    print("tau_int(e) =", tau_e)
    print("tau_int(|m|) =", tau_abs_m)


exercise3a()
#corr_plot()



