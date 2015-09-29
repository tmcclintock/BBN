import ctypes
import numpy as np
import matplotlib.pyplot as plt

forward_titles = [r"$n\rightarrow p$",#0
                  r"$^3H\rightarrow e^-+\nu+^3He$",#1
                  r"$^8Li\rightarrow e^- + \nu + 2^4He$",#2
                  r"$^{12}B\rightarrow e^- + \nu + ^{12}C$",#3
                  r"$^{14}C \rightarrow e^- +\nu + ^{14}N$",#4
                  r"$^8B\rightarrow e^+ \nu+2^{4}He$",#5
                  r"$^{11}C\rightarrow e^++\nu+^{11}B$",#6
                  r"$^{12}N\rightarrow e^++\nu+^{12}C$",#7
                  r"$^{13}N\rightarrow e^+|\nu+^{13}C$",#8
                  r"$^{14}O\rightarrow e^+|\nu+^{14}N$",#9
                  r"$^{15}O\rightarrow e^++\nu+^{15}N$",#10, LAST WEAK ONE
                  r"$p+n\rightarrow \gamma+D$",#11
                  r"$D+p\rightarrow\gamma+\!\,^3\!He$",#12
                  r"$D+D\rightarrow n+\!\,^3\!He$",#13
                  r"$D+D\rightarrow p+T$",#14
                  r"$^3\!He+n\rightarrow p+T$",#15
                  r"$T+D\rightarrow n+\!\,^4\!He$",#16
                  r"$^3\!He+D\rightarrow p+\!\,^4\!He$",#17
                  r"$^3\!He+\!\,^4\!He\rightarrow\gamma+\!\,^7\!Be$",#18
                  r"$T+\!\,^4\!He\rightarrow\gamma+\!\,^7\!Li$",#19
                  r"$^7\!Be+n\rightarrow p+\!\,^7\!Li$",#20
                  r"$^7\!Li+p\rightarrow\!\,^4\!He+\!\,^4\!He$",#21
                  r"$D+n\rightarrow\gamma+T$",#22
                  r"$^3\!He+n\rightarrow\gamma+\!\,^4\!He$",#23
                  r"$^7\!Be+n\rightarrow\!\,^4\!He+\!\,^4\!He$"#24
                  r"$T+p\rightarrow\gamma+\!\,^4\!He$",#25
                  r"$^3\!He+\!\,^3\!He\rightarrow p+p+\!\,^4\!He$",#26
                  r"$^7\!Li+D\rightarrow n+\!\,^4\!He+\!\,^4\!^4He$",#27
                  r"$^7\!Be+D\rightarrow p+\!\,^4\!He+\!\,^4\!^4He$",#28
                  r"$D+\!\,^4\!He\rightarrow p+\!\,^6\!Li$",#29
                  r"$^6\!Li+n\rightarrow\gamma+\!\,^7\!Li$",#30
                  r"$^6\!Li+n\rightarrow\!\,^4\!He+T$",#31
                  r"$^6\!Li+p\rightarrow\gamma+\!\,^7\!Be$",#32
                  r"$^6\!Li+p\rightarrow\!\,^4\!He+\!\,^3\!He$"#33
                  ]

def test_rates():
    print "Testing nuclear reaction rates"
    #Load the library
    rr_lib = ctypes.cdll.LoadLibrary('../Reaction_rates/reaction_rates.so')
    rr = rr_lib.reaction_rates #function pointer
    #Set the argument types (double T9, double*f, double*r)
    rr.argtyps = [ctypes.c_double,\
                  ctypes.POINTER(ctypes.c_double),\
                  ctypes.POINTER(ctypes.c_double)]
    #Get the number of reactions
    nreac = ctypes.c_int.in_dll(rr_lib, 'totalnreac_alias').value
    print "\tNumber of reactions = ",nreac

    #Create the temperature array
    ntemps = 100
    T = 10**(np.arange(np.log10(100),np.log10(0.1),\
                         (np.log10(0.1)-np.log10(100))/ntemps))
    
    #Create the forward and backwards arrays
    forward = np.zeros((ntemps,nreac))
    reverse = np.zeros((ntemps,nreac))

    #Loop over temperatures to get the rates
    for i,T9 in zip(range(0,len(T)),T):
        f = np.zeros(nreac,dtype=np.double)
        fin = f.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        r = np.zeros((nreac),dtype=np.double)
        rin = r.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        #Call the reaction rates function
        rr(ctypes.c_double(T9),fin,rin)
        f = np.array(fin[0:nreac])
        r = np.array(rin[0:nreac])
        forward[i,:]=f
        reverse[i,:]=r
        continue

    for i in range(0,len(forward_titles)):#nreac):
        try:
            plt.clf()
            plt.loglog(T,forward[:,i])
            plt.xlabel("$T_9 [GK]$")
            plt.ylabel("rate")
            plt.title(forward_titles[i])
            plt.xlim(max(T),min(T))
            plt.gcf().savefig("reaction_pngs/forwards/forward_reaction_"\
                              +str(i)+".png")
            plt.clf()
            plt.loglog(T,reverse[:,i])
            plt.xlabel("$T_9 [GK]$")
            plt.ylabel("rate")
            plt.title("Reverse: "+forward_titles[i])
            plt.xlim(max(T),min(T))
            plt.gcf().savefig("reaction_pngs/reverses/reverse_reaction_"\
                              +str(i)+".png")
            plt.clf()
            plt.loglog(T,forward[:,i])
            plt.loglog(T,reverse[:,i],ls=":")
            plt.xlabel("$T_9 [GK]$")
            plt.ylabel("rate")
            plt.title("Both: "+forward_titles[i])
            plt.xlim(max(T),min(T))
            plt.gcf().savefig("reaction_pngs/boths/both_reaction_"\
                              +str(i)+".png")
            plt.clf()
            #plt.show()
        except ValueError:
            continue
        continue

if __name__ == "__main__":
    test_rates()
