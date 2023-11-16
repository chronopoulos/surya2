#!/usr/bin/python3

from PyQt5 import QtWidgets
import pyqtgraph as pg

import numpy as np
from scipy.interpolate import griddata

import os

DATA_DIR = os.path.expanduser('~/data/surya2/serial')
DATA_FILE = os.path.expanduser('diffrot.dat')

N = 128 + 2
RMIN = 0.382794500E+11
RMAX = 0.695990000E+11

class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, data):
        super().__init__()

        t =     data[:,0]
        r =     data[:,1]
        omega =   data[:,2]

        x = r * np.sin(t)
        y = RMAX - r*np.cos(t)

        xf = np.linspace(0, RMAX, 2*N)
        yf = np.linspace(0, 2*RMAX, 4*N)
        xxf, yyf = np.meshgrid(xf, yf, indexing='ij')

        omega_rect = griddata((x, y), omega, (xxf, yyf), method='linear')

        self.imageView = pg.ImageView()
        self.imageView.setImage(omega_rect)
        self.setCentralWidget(self.imageView)

        self.resize(600, 800)


if __name__ == '__main__':

    data = np.loadtxt(os.path.join(DATA_DIR, DATA_FILE))

    app = QtWidgets.QApplication([])
    main = MainWindow(data)
    main.show()
    app.exec()
