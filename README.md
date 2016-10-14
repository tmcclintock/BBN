# BBN

This is my implementation of a Big Bang Nucleosythesis calculation. The
advantage over other codes is that it uses a higher order integrator (RK4)
than some others which are lower order. On top of that the
integration variable is time rather than redshift, which to me
is a bit easier to understand.

The directory Old_LCDM_Version contains my original implementation
and isn't commented well. I don't recommend using this.

The directory Python_C_LCDM contains a C implementation that
is compiled into a shared library to be used in Python.

# Compilation
Just cd into Python_C_LCDM/src/ and compile with $make.

# Usage
Look at the file Python_C_LCDM/pybbn.py for how to interface with
the library. One day in the hopefully not-too distant future
an actual wrapper routine will be written that can do this
more cleanly.

Essentially, you pass in the baryon to photon ratio (eta0)
as well as an array for the final abundances (Y_final). The two
other arguments are to_print (tp) and the print_increment (pi)
which control how often output is printed to file, if at all.

The real interesting data gets dumped into Python_C_LCDM/output_files/
where you can find dynamics.txt and abundances.txt. 

The dynamics.txt file has as a function of 
time the temperature, hubble constant,
and electron chemical potential (I think I remember that's what it's called).

The abundances.txt file has abundances of all species for each step writen
out to the dynamics file (so as a function of time).