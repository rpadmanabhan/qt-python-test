import sys

from PySide6.QtWidgets import QApplication, QMainWindow
import pyqtgraph as pg
import numpy as np

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
        N = 300_000  # Total number of points
        K = 10       # Number of clusters
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
        plt.addItem(self.roi)

    def scatter_point_tip(self, x, y, data):
        return f'Point {data} = ({x:0.2f}, {y:0.2f})'

app = QApplication(sys.argv)
w = MainWindow()
w.show()
app.exec()