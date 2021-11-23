from matplotlib import pyplot as plt

def Plotting_results(x,y,n,types,ref):


    for i in range(0,n):

        plt.plot(x,(y[i,:]-ref)/ref*100)
        plt.xlabel("q(a.u.)")
        plt.ylabel("%\Delta I")
        plt.legend(types[i])
