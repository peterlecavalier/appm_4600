from mypkg.my2DPlot import my2DPlot as myplt
import numpy as np
plt = myplt(lambda x : np.sin(x) + 1,0.,10.)
plt.labels('x','y')
plt.addPlot(lambda x : np.cos(x) + 1)
plt.dotted()
plt.color('black')
plt.logy()
plt.logx()
plt.save('figure.pdf')
plt.show()