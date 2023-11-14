#!/usr/bin/python3

import numpy as np
import pyqtgraph as pg
import os

DATA_DIR = os.path.expanduser('~/data/surya2/serial')

data = np.loadtxt(os.path.join(DATA_DIR, 'run.dat'))

t  = data[:,0]
c1 = data[:,1]
c2 = data[:,2]
c3 = data[:,3]
c4 = data[:,4]
c5 = data[:,5]

w = pg.GraphicsLayoutWidget()

# Add subplots
w.addPlot(x=t, y=c1, row=0, col=0)
w.addPlot(x=t, y=c2, row=1, col=0)
w.addPlot(x=t, y=c3, row=2, col=0)
w.addPlot(x=t, y=c4, row=3, col=0)
w.addPlot(x=t, y=c5, row=4, col=0)

w.show()

pg.Qt.QtGui.QGuiApplication.instance().exec_()
