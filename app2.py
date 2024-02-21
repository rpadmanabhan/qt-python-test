import sys

from PySide6.QtWidgets import QApplication, QMainWindow
import pyqtgraph as pg
import numpy as np

class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        N = 500
        self.data = np.empty(N, dtype=[('hour', int), ('temp', float)])
        self.data['hour'] = np.arange(0, N)
        self.data['temp'] = np.linspace(35, 100, N)
        print(f"Done generating {N} datapoints")
        self.scatterPlotWidget = pg.ScatterPlotWidget()
        self.scatterPlotWidget.setData(self.data)
        self.scatterPlotWidget.setFields([('hour', {'units': 'h'}), ('temp', {'units': 'C'})])
        self.setCentralWidget(self.scatterPlotWidget)
        self.scatterPlotWidget.show()


app = QApplication(sys.argv)
w = MainWindow()
w.show()
app.exec()