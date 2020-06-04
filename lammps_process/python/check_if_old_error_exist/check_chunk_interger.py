import sys
import numpy as np

import read_setting as rr
import plotfigure.plotchunk1 as pp

def checkNchunkint():
    if int(rr.logfile["freq_ave_chunk_momentum_mass_field"]) != 0:
        ob1 = pp.chunkstress(1, rr.lammps_directory)
        for step_inloop in ob1.allsteps:
            n_line_0 = int(int(step_inloop - ob1.step_first_in_file)/ob1.d_step)*(ob1.n_line_in_a_step+1) + 4
            n_line_1 = int(n_line_0 + ob1.n_line_in_a_step)
            ## select data
            Nchunklist = [ob1.lines[t].split()[ob1.header.index('Ncount')] for t in range(n_line_0, n_line_1)]
            for x in Nchunklist:
                try: 
                    int(x)
                except ValueError:
                    print("one of chunk in step {step_inloop} contain {number} Ncount".format(step_inloop=step_inloop, number=x))
                    if float(x) < 0.1:
                        breakpoint()
                        Nchunkarray = np.asarray(Nchunklist, dtype=np.float64, order='F')
                        totalNchunk = np.sum(Nchunkarray)
                        print(totalNchunk)
                        sys.exit("Ncount (number of atom in chunk) not int")
        print("all chunk has integer number of atoms")

# main exclusive
if __name__ == "__main__":
    checkNchunkint()