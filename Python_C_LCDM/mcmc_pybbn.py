import ctypes as ctypes
import numpy as np
import matplotlib.pyplot as plt
import emcee

#run the program once
bbn_lib = ctypes.cdll.LoadLibrary("src/pybbn.so")
driver = bbn_lib.driver #function pointer
driver.argtyps = [ctypes.c_int,ctypes.c_int,ctypes.c_double]
driver.restype = ctypes.c_int

to_print, print_increment = 0,1 #1 is true, for to_print

#Declare the array that holds the final number densities (mislabeled as Y)
totalnnuc = ctypes.c_int.in_dll(bbn_lib,'totalnnuc_alias').value
Y_final = np.zeros(totalnnuc,dtype=np.double)
Y_final_ct = Y_final.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

#driver(ctypes.c_int(tp),ctypes.c_int(1),ctypes.c_double(eta0_guess),Y_final_ct)
#Y_final = np.array(Y_final_ct[0:totalnnuc])

obs = np.array([2.53e-5,0.2465,1.6e-10]) #D,He4,Li7
sigma = np.array([0.04e-5,0.0097,.3e-10])#D,He4,Li7
norm = 0.5*np.log(2*np.pi*sigma**2)


def lnprob(params,obs,sigma,norm):
    #The guess 
    leta0_guess = params[0]
    #The priors
    if (leta0_guess > -8 or leta0_guess < -12):
        return -1e99 #a bad likelihood
    eta0_guess = 10**leta0_guess

    #Run our model to get the final abundances
    Y_final = np.zeros(totalnnuc,dtype=np.double)
    Y_final_ct = Y_final.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    driver(ctypes.c_int(to_print),ctypes.c_int(print_increment),\
           ctypes.c_double(eta0_guess),Y_final_ct)
    Y_final = np.array(Y_final_ct[0:totalnnuc])
    Y = Y_final

    #Calculate abundances
    YD = Y[2]/Y[1] #D
    YHe3 = (Y[4]+Y[3])/Y[1] #He3
    XHe4 = Y[5]*4 #He4 mass fraction
    YLi7 = (Y[6]+Y[7])/Y[1] #Li7
    YLi6 = Y[8]/Y[1] #Li6

    #Compare to observations

    #Make an array of our model
    model = np.array([YD,XHe4,YLi7])

    like = -(model-obs)**2/sigma**2/2.0 - norm

    return np.sum(like)

#do the emcee stuff
ndim, nwalkers = 1, 8
eta0_best = 6.19e-10
leta0_best = np.log10(eta0_best)
p0 =[leta0_best+1e-4*np.random.rand(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,args=(obs,sigma,norm))

#Run the burnin
nburnin = 100
pow,prob,state = sampler.run_mcmc(p0,nburnin)

#Reset the sampler
sampler.reset()

#Now run the actual sampling
ntrials = 200
sampler.run_mcmc(pow,ntrials,rstate0=state)

#Print out some useful information
print "Mean acceptance fraction: ",np.mean(sampler.acceptance_fraction)
print "Autocorrelation time: ", sampler.get_autocorr_time()

#Save and output the chain
flatchain = sampler.flatchain
flatlnprobs = sampler.flatlnprobability
output_file = open("output_files/chain.txt","w")
output_file.write("eta\tloglike\n")
for i in range(len(flatlnprobs)):
    output_file.write(str(flatchain[i,0])+"\t"+str(flatlnprobs[i])+"\n")
