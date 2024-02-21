import sys

from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QMenu, QDockWidget, QListWidget, QInputDialog, QPushButton,
    QLayout, QVBoxLayout, QWidget, QAbstractItemView
)
from PySide6.QtCore import (
    QObject, QRunnable, QThread, QThreadPool, Qt, Signal, Slot,
)


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

class MainAppSignals(QObject):
    selectedPointsUpdated = Signal(list)


class CustomPlotWidget(pg.GraphicsLayoutWidget):
    def __init__(self, labels_list_widget, *args, **kwargs):
        super(CustomPlotWidget, self).__init__(*args, **kwargs)
        self.labels_list_widget = labels_list_widget
        self.selected_points = []  # Initialize with an empty list
        self.labels_dict = {}
        self.setMouseTracking(True)  # Optional: For more interactive features
    
    def mousePressEvent(self, event):
        if event.button() == Qt.RightButton:
            self.showContextMenu(event.pos())
        else:
            super(CustomPlotWidget, self).mousePressEvent(event)
    
    def showContextMenu(self, position):
        menu = QMenu()
        assign_label_action = QAction('Assign Label', self)
        assign_label_action.triggered.connect(self.assignLabelToSelectedPoints)
        menu.addAction(assign_label_action)
        menu.exec_(self.mapToGlobal(position))
    
    def assignLabelToSelectedPoints(self):
        # This slot will be triggered when the QAction is activated
        label, ok = QInputDialog.getText(self, "Label Points", "Enter label name:")
        if ok and label:
            self.labels_dict[label] = self.selected_points
            self.updateLabelList()

    def updateLabelList(self):
        self.labels_list_widget.clear()
        for label, points in self.labels_dict.items():
            self.labels_list_widget.addItem(f"{label}: {len(points)}")

    def updateSelectedPoints(self, selected_points):
        self.selected_points = selected_points
        # Handle the update, e.g., refresh the view if necessary


class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.threadpool = QThreadPool()
        self.selected_points = []

        # Initialization code...
        self.signals = MainAppSignals()
        # Setup the dock widget for displaying labels, run DE, etc.
        self.dock_widget = QDockWidget("Analysis Options", self)        
        self.dock_layout        = QVBoxLayout()
        self.dock_layout.setSpacing(100)
        self.multi_widget       = QWidget()
        self.labels_list_widget = QListWidget()
        self.labels_list_widget.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.labels_list_widget.setFixedSize(200, 100)

        self.runDE_widget       = QPushButton("Run Differential Expression", disabled = True)

        self.dock_layout.addWidget(self.labels_list_widget)
        self.dock_layout.addWidget(self.runDE_widget)

        self.multi_widget.setLayout(self.dock_layout)

        self.dock_widget.setWidget(self.multi_widget)

        self.addDockWidget(Qt.RightDockWidgetArea, self.dock_widget)


        N = 300_000  # Total number of points
        K = 5       # Number of clusters
        self.data = generate_points(N, K)

        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')
        sp = pg.ScatterPlotItem(size = 6, pen = pg.mkPen(None), brush = pg.mkBrush(0, 0, 0, 20),
                                hoverable=True, hoverSymbol='o', hoverSize=8,
                                hoverBrush=pg.mkBrush('g'), tip = self.scatter_point_tip)
        sp.addPoints(x = self.data['pressure'], y = self.data['temp'], data = np.arange(N))

        view = CustomPlotWidget(self.labels_list_widget)
        self.setCentralWidget(view)
        self.signals.selectedPointsUpdated.connect(view.updateSelectedPoints)


        plt = view.addPlot()
        plt.addItem(sp)

        self.roi = pg.RectROI(
            [8, 8], size = 5, centered=True, sideScalers=True,
            pen = pg.mkPen('r', width = 3, style = Qt.DotLine),
            hoverPen = pg.mkPen('r', width = 3, style = Qt.DashLine),
        )
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
        self.selected_points = selected
        self.signals.selectedPointsUpdated.emit(selected)
        # Remember to update GUI elements like plots here
                
app = QApplication(sys.argv)
w = MainWindow()
w.show()
app.exec()