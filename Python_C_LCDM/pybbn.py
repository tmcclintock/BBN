import ctypes as ctypes
import numpy as np
import matplotlib.pyplot as plt
import time

#run the program once
bbn_lib = ctypes.cdll.LoadLibrary("src/pybbn.so")
eta0_guess = 6.19e-10 #Baryon to photon ratio
driver = bbn_lib.driver #function pointer
driver.argtyps = [ctypes.c_int,ctypes.c_int,ctypes.c_double]
driver.restype = ctypes.c_int
tp,pi=1,1 #to_print and print_increment
print "Printing output (1 is true): ",tp
if tp:
    print "Printing increment: ",pi
print "Eta0 guess: ",eta0_guess

totalnnuc = ctypes.c_int.in_dll(bbn_lib,'totalnnuc_alias').value
Y_final = np.zeros(totalnnuc,dtype=np.double) #Final abundances
Y_final_ct = Y_final.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
start = time.time()
driver(ctypes.c_int(tp),ctypes.c_int(1),ctypes.c_double(eta0_guess),Y_final_ct)
end = time.time()
print "BBN calculation time = %f"%(end-start)
Y_final = np.array(Y_final_ct[0:totalnnuc])

print Y_final

#get the data
if(tp==1):
    dynamics_data = np.genfromtxt("output_files/dynamics.txt",names=True)
    abundances_data = np.genfromtxt("output_files/abundances.txt",names=True)
    
    T9 = dynamics_data['T9']
    p = abundances_data['p']
    n = abundances_data['n']
    D = abundances_data['D']
    T = abundances_data['T']
    He3 = abundances_data['He3']
    He4 = abundances_data['He4']
    Be7 = abundances_data['Be7']
    Li7 = abundances_data['Li7']
    Li6 = abundances_data['Li6']
    
    #The input data are raw number per volume
    #abundances are measured w.r.t. the ratio with H
    #Y is an abundance, X is a mass fraction
    Yn = n/p
    YD = D/p
    XHe4 = He4*4
    YHe3 = (He3+T)/p
    YLi7 = (Be7+Li7)/p
    YLi6 = Li6/p
    
    plt.loglog(T9,Yn,label=r"$Y_n$")
    plt.loglog(T9,YD,label=r"$D$")
    plt.loglog(T9,XHe4,label=r"$X_\alpha$")
    plt.loglog(T9,YHe3,label=r"$Y_{\,^3\!\,He}$")
    plt.loglog(T9,YLi7,label=r"$Y_{\,^7\!\,Li}$")
    plt.loglog(T9,YLi6,label=r"$Y_{\,^6\!\,Li}$")
    
    plt.xlim(max(T9),min(T9))
    plt.ylim(1e-25,10)
    plt.legend(loc='lower right')
    
    plt.show()
    
