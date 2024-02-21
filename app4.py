import sys

from PySide6.QtWidgets import QApplication, QMainWindow
from PySide6.QtCore import QObject, QRunnable, QThread, QThreadPool, Signal, Slot
import pyqtgraph as pg
import numpy as np

class ROIWorkerSignals(QObject):
    finished = Signal()
    result = Signal(list)
    error = Signal(str)

class ROIWorker(QRunnable):

    def __init__(self, data, roiBounds):
        super(ROIWorker, self).__init__()
        self.data = data
        self.roiBounds = roiBounds
        self.signals = ROIWorkerSignals()

    @Slot()
    def run(self):
        selected = []

        for i, (x, y) in enumerate(zip(self.data['pressure'], self.data['temp'])):
            if self.roiBounds.contains(pg.Point(x, y)):
                selected.append(i)
        self.signals.result.emit(selected)

def generate_points(N, K):
    points_per_cluster = N // K

    # Generate random means and covariance matrices for each cluster
    cluster_means = np.random.uniform(low=5, high=50, size=(K, 2))
    cluster_covs = np.random.uniform(low=0.5, high=2, size=(K, 2, 2))

    # Generate data points for each cluster
    data = np.empty(N, dtype=[('pressure', float), ('temp', float)])
    for i in range(K):
        start_index = i * points_per_cluster
        end_index = (i + 1) * points_per_cluster
        mvn_samples = np.random.multivariate_normal(mean=cluster_means[i], cov=cluster_covs[i], size=points_per_cluster)
        data['pressure'][start_index:end_index] = mvn_samples[:, 0]
        data['temp'][start_index:end_index] = mvn_samples[:, 1]

    print(f"Done generating {N} datapoints with {K} clusters.")
    return data

class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.threadpool = QThreadPool()

        N = 100_000  # Total number of points
        K = 5       # Number of clusters
        self.data = generate_points(N, K)

        sp = pg.ScatterPlotItem(size = 10, pen = pg.mkPen(None), brush = pg.mkBrush(255, 255, 255, 20),
                                hoverable=True, hoverSymbol='o', hoverSize=8,
                                hoverBrush=pg.mkBrush('g'), tip = self.scatter_point_tip)
        sp.addPoints(x = self.data['pressure'], y = self.data['temp'], data = np.arange(N))

        view = pg.GraphicsLayoutWidget()
        self.setCentralWidget(view)

        plt = view.addPlot()
        plt.addItem(sp)

        self.roi = pg.RectROI([8, 8], size = 5)
        self.roi.sigRegionChangeFinished.connect(self.onROIChanged)
        plt.addItem(self.roi)

    def scatter_point_tip(self, x, y, data):
        return f'Point {data} = ({x:0.2f}, {y:0.2f})'
    
    def onROIChanged(self):
        roiBounds = self.roi.parentBounds()
        worker = ROIWorker(self.data, roiBounds)
        worker.signals.result.connect(self.updateGUIWithSelectedPoints)
        self.threadpool.start(worker)

    def updateGUIWithSelectedPoints(self, selected):
        # This function runs in the main thread
        # Update your GUI based on the selected points here
        print(f"Selected points: {len(selected)}")
        # Remember to update GUI elements like plots here
                
app = QApplication(sys.argv)
w = MainWindow()
w.show()
app.exec()