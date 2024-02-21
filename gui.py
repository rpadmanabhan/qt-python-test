import os
import sys

from PySide6.QtGui import QAction, QColor, QFont
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QMenu, QDockWidget, QListWidget, QInputDialog, QPushButton,
    QLayout, QVBoxLayout, QHBoxLayout, QWidget, QAbstractItemView, QFileDialog,
    QRadioButton, QListWidgetItem, QDialog, QLabel
)
from PySide6.QtCore import (
    QObject, QRunnable, QThread, QThreadPool, Qt, Signal, Slot
)
import pyqtgraph as pg
import numpy as np

from backend import load_h5mu, compute_tsne, compute_umap

# Define colors as constants
COLOR_WHITE = (0, 0, 0)   # Black
COLOR_RED = (255, 0, 0)         # Red
COLOR_GREEN = (0, 255, 0)       # Green
COLOR_BLUE = (0, 0, 255)        # Blue
COLOR_YELLOW = (255, 255, 0)    # Yellow
COLOR_CYAN = (0, 255, 255)      # Cyan
COLOR_MAGENTA = (255, 0, 255)   # Magenta
COLOR_ORANGE = (255, 165, 0)    # Orange
COLOR_PURPLE = (128, 0, 128)    # Purple
COLOR_TEAL = (0, 128, 128)      # Teal
PREDEFINED_COLORS = [
    COLOR_WHITE, COLOR_RED, COLOR_GREEN, COLOR_BLUE, COLOR_YELLOW,
    COLOR_CYAN, COLOR_MAGENTA, COLOR_ORANGE, COLOR_PURPLE, COLOR_TEAL
]

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')


class ROIWorkerSignals(QObject):
    finished = Signal()
    result = Signal(list)
    error = Signal(str)

class ROIWorker(QRunnable):

    def __init__(self, data, roiBounds, dimred_plot):
        super(ROIWorker, self).__init__()
        self.data = data
        self.dimred_plot = dimred_plot
        self.roiBounds = roiBounds
        self.signals = ROIWorkerSignals()

    def run(self):
        selected_cells = []
        for i, point_data in enumerate(self.dimred_plot.data):
            x, y = point_data[0], point_data[1]
            if self.roiBounds.contains(pg.Point(x, y)):
                selected_cells.append(point_data[7])
        self.signals.result.emit(selected_cells)

class MainAppSignals(QObject):
    selectedPointsUpdated = Signal(list)

class UmapRunnable(QRunnable):
    def __init__(self, adata, color_map):
        super().__init__()
        self.adata = adata
        self.color_map = color_map

        self.dialog = QDialog()
        self.layout = QVBoxLayout()
        self.graph_widget = pg.GraphicsLayoutWidget()
        self.label = QLabel("")

    def run(self):
        # Perform UMAP computation here
        self.update_progress("UMAP Computation in Progress...")

        # Perform UMAP computation
        compute_umap(self.adata)

        # Draw the UMAP scatter plot in the same window
        self.draw_dimred_scatter(dimred='UMAP')

        # Display "Done computing UMAP" message
        self.update_progress("Done computing UMAP")



    def draw_dimred_scatter(self, dimred='t-SNE'):
        self.layout.addWidget(self.graph_widget)
        self.cur_dimred_type = dimred

        scatter_plot = pg.ScatterPlotItem(
            size=6, pen=pg.mkPen(None), brush=pg.mkBrush(0, 0, 0, 20),
            hoverable=True, hoverSymbol='o', hoverSize=4,
            hoverBrush=pg.mkBrush('w'), tip=self.scatter_point_tip)
        if dimred == 't-SNE':
            scatter_plot.addPoints(
                x=self.adata.obsm['X_tsne'][:, 0],
                y=self.adata.obsm['X_tsne'][:, 1],
                data=self.adata.obs.index,
                brush=[self.color_map[cell_type] for cell_type in self.adata.obs['Cell_Type_Experimental']]
            )
        else:
            scatter_plot.addPoints(
                x=self.adata.obsm['X_umap'][:, 0],
                y=self.adata.obsm['X_umap'][:, 1],
                data=self.adata.obs.index,
                brush=[self.color_map[cell_type] for cell_type in self.adata.obs['Cell_Type_Experimental']]
            )

        # Add the scatter plot to the layout
        self.current_plot = self.graph_widget.addPlot()
        self.current_plot.addItem(scatter_plot)
            

    def scatter_point_tip(self, x, y, data):
        return f'Cell {data}'

    def update_progress(self, message):
        self.label.setText(message)
        QApplication.processEvents()

        if not self.dialog.isVisible():
            self.layout.addWidget(self.label)
            self.dialog.setLayout(self.layout)
            self.dialog.setWindowTitle("UMAP Progress")
            self.dialog.show()


class TsneRunnable(QRunnable):
    def run(self):
        # Perform t-SNE computation here
        # For now, let's just open a new window
        dialog = QDialog()
        layout = QVBoxLayout()
        label = QLabel("t-SNE Computation in Progress...")
        layout.addWidget(label)
        dialog.setLayout(layout)
        dialog.exec()


class CustomPlotWidget(pg.GraphicsLayoutWidget):
    def __init__(self, labels_list_widget, labels_dict, *args, **kwargs):
        super(CustomPlotWidget, self).__init__(*args, **kwargs)
        self.labels_list_widget = labels_list_widget
        self.labels_dict = labels_dict
        self.selected_points = []
        self.setMouseTracking(True)
    
    def mousePressEvent(self, event):
        if event.button() == Qt.RightButton:
            self.showContextMenu(event.pos())
        else:
            super(CustomPlotWidget, self).mousePressEvent(event)
    
    def showContextMenu(self, position):
        menu = QMenu()
        assign_label_action = QAction('Assign Label', self)
        assign_label_action.triggered.connect(
            self.assignLabelToSelectedPoints)
        menu.addAction(assign_label_action)
        menu.exec_(self.mapToGlobal(position))
    
    def assignLabelToSelectedPoints(self):
        # This slot will be triggered when the QAction is activated
        label, ok = QInputDialog.getText(
            self, "Label Points", "Enter label name:")
        if ok and label:
            label_disp = f'{label}: {len(self.selected_points)} cells'
            self.labels_dict[label_disp] = self.selected_points
            self.labels_list_widget.addItem(label_disp)

    def updateSelectedPoints(self, selected_points, dimred_plot):
        self.selected_points = selected_points

class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()

        
        # Common app related
        #----------------------------------------------------------------------------------------        
        self.threadpool = QThreadPool()
        self.selected_file = None
        self.dataset = {}
        self.signals = MainAppSignals()
        self.labels_dict = {}
        self.dimred_buttons_added = False
        self.setWindowTitle('tSNE/UMAP visualizer')

        # Side Dock
        #----------------------------------------------------------------------------------------
        self.dock_widget = QDockWidget('', self)
        self.dock_container = QWidget()
        self.dock_layout = QVBoxLayout()
        self.addDockWidget(Qt.RightDockWidgetArea, self.dock_widget)

        self.dock_container.setLayout(self.dock_layout)
        self.dock_widget.setWidget(self.dock_container)

        self.file_button = QPushButton("Load h5mu file")
        self.file_dialog = QFileDialog(self)
        self.file_dialog.setNameFilters(['h5mu (*.h5mu)', 'h5ad (*.h5ad)'])
        self.file_button.clicked.connect(self.load_file)
        self.dock_layout.addWidget(self.file_button)

        self.dataset_list_widget = QListWidget()
        self.dataset_list_widget.setSelectionMode(
            QAbstractItemView.ExtendedSelection)
        self.dataset_list_widget.setFixedSize(200, 50)
        self.dock_layout.addWidget(self.dataset_list_widget)

        self.info_section_layout = QVBoxLayout()
        self.info_section_widget = QWidget()
        self.info_section_widget.setLayout(self.info_section_layout)
        self.dock_layout.addWidget(self.info_section_widget)

        self.info_list_widget = QListWidget()
        self.info_list_widget.setSelectionMode(
            QAbstractItemView.ExtendedSelection)
        self.info_list_widget.setFixedSize(200, 200)
        #self.info_list_widget.setStyleSheet("background-color: darkgray;")
        self.info_list_widget.itemSelectionChanged.connect(self.enable_dimred_buttons)
        self.info_list_widget.itemSelectionChanged.connect(self.handleInfoListSelection)
        self.info_section_layout.addWidget(self.info_list_widget)

        self.dock_button_layout  = QHBoxLayout()
        self.dock_button_widgets = QWidget()
        self.dock_button_widgets.setLayout(self.dock_button_layout)
        self.dock_layout.addWidget(self.dock_button_widgets)

        self.select_cells       = QPushButton(
            'Select Cells', disabled = True)        
        self.dock_button_layout.addWidget(self.select_cells)
        self.compute_UMAP       = QPushButton(
            'Compute UMAP', disabled = True)
        self.compute_TSNE       = QPushButton(
            'Compute t-SNE', disabled = True)
        self.select_cells.clicked.connect(self.add_roi_selector)
        self.compute_UMAP.clicked.connect(self.compute_umap)
        self.compute_TSNE.clicked.connect(self.compute_tsne)       
        self.dock_button_layout.addWidget(self.compute_UMAP)
        self.dock_button_layout.addWidget(self.compute_TSNE)

        # Graphing related stuff       
        #----------------------------------------------------------------------------------------
        self.graph_widget = CustomPlotWidget(self.info_list_widget, self.labels_dict)
        self.current_plot = None
        self.scatter_plot = None
        self.roi          = None
        self.setCentralWidget(self.graph_widget)
        self.signals.selectedPointsUpdated.connect(
            lambda selected_points: \
                self.graph_widget.updateSelectedPoints(
                    selected_points, self.scatter_plot)
        )

    def onROIChanged(self):
        roiBounds = self.roi.parentBounds()
        worker = ROIWorker(self.dataset, roiBounds, self.scatter_plot)
        worker.signals.result.connect(self.updateGUIWithSelectedPoints)
        self.threadpool.start(worker)

    def updateGUIWithSelectedPoints(self, selected):
        # This function runs in the main thread
        # Update your GUI based on the selected points here
        print(f"Selected points: {len(selected)}")
        self.selected_points = selected
        self.signals.selectedPointsUpdated.emit(selected)

    def add_roi_selector(self):
        self.roi = pg.RectROI(
            [8, 8], size = 5, centered=True, sideScalers=True,
            pen = pg.mkPen('r', width = 3, style = Qt.DotLine),
            hoverPen = pg.mkPen('r', width = 3, style = Qt.DashLine),
        )
        self.roi.sigRegionChangeFinished.connect(self.onROIChanged)
        self.current_plot.addItem(self.roi)

    def draw_dimred_scatter(self, dimred = 't-SNE'):
        self.cur_dimred_type = dimred
        self.graph_widget.clear()
        unique_cell_types = self.dataset['rna'].obs['Cell_Type_Experimental'].unique()
        self.color_map = {cell_type: (*color, 60) for color, cell_type in zip(PREDEFINED_COLORS, unique_cell_types)}

        self.scatter_plot = pg.ScatterPlotItem(
            size = 6, pen = pg.mkPen(None), brush = pg.mkBrush(0, 0, 0, 20),
            hoverable=True, hoverSymbol='o', hoverSize=4,
            hoverBrush=pg.mkBrush('w'), tip = self.scatter_point_tip)
        if dimred == 't-SNE':
            self.scatter_plot.addPoints(
                x = self.dataset['rna'].obsm['X_tsne'][:,0],
                y = self.dataset['rna'].obsm['X_tsne'][:,1],
                data = self.dataset['rna'].obs.index,
                brush=[self.color_map[cell_type] for cell_type in self.dataset['rna'].obs['Cell_Type_Experimental']]

            )
        else:
            self.scatter_plot.addPoints(
                x = self.dataset['rna'].obsm['X_umap'][:,0],
                y = self.dataset['rna'].obsm['X_umap'][:,1],
                data = self.dataset['rna'].obs.index,
                brush=[self.color_map[cell_type] for cell_type in self.dataset['rna'].obs['Cell_Type_Experimental']]
            )
        self.current_plot = self.graph_widget.addPlot()
        self.current_plot.addItem(self.scatter_plot)


    def handleInfoListSelection(self):
        # Get the currently selected items in the info list widget
        selected_items = self.info_list_widget.selectedItems()

        if len(selected_items) == 0:
            # No item is selected, repaint the graph
            self.repaintGraph()
        else:
            # Extract the text (cell types/labels) of the selected items
            selected_cell_types = [item.text() for item in selected_items]

            # Get the indices of points corresponding to the selected cell types
            selected_indices = []
            for cell_type in selected_cell_types:
                if cell_type in self.labels_dict:
                    selected_indices.extend(self.labels_dict[cell_type])

            # Highlight the selected points in the scatter plot
            self.highlightSelectedPoints(selected_indices)

    def repaintGraph(self):
        # Repaint the graph with the original colors from the color map
        #self.draw_dimred_scatter(dimred = self.cur_dimred_type)
        self.scatter_plot.setPointsVisible(True)

    def highlightSelectedPoints(self, selected_indices):
        # Highlight the selected points
        self.scatter_plot.setPointsVisible(
            self.dataset['rna'].obs.index.isin(selected_indices))

    def scatter_point_tip(self, x, y, data):
        return f'Cell {data}'

    def load_file(self):
        self.file_dialog.exec()
        self.selected_file = self.file_dialog.selectedFiles()[0]
        self.dataset_list_widget.clear()
        fname_base = os.path.basename(self.selected_file)
        self.dataset_list_widget.addItem(fname_base)
        self.dataset = load_h5mu(self.selected_file)
        self.add_cell_types_to_info_list()

        if not self.dimred_buttons_added:
            self.tsne_button = QRadioButton('t-SNE')
            self.umap_button = QRadioButton('UMAP')
            self.tsne_button.clicked.connect(
                lambda: self.draw_dimred_scatter(dimred='t-SNE'))
            self.umap_button.clicked.connect(
                lambda: self.draw_dimred_scatter(dimred='UMAP'))
            self.info_section_layout.addWidget(self.tsne_button)
            self.info_section_layout.addWidget(self.umap_button)
            self.dimred_buttons_added = True

        self.draw_dimred_scatter('t-SNE')
        self.select_cells.setEnabled(True)

    def add_cell_types_to_info_list(self):
        # Clear existing items from the info list widget
        self.info_list_widget.clear()

        # Add items to the info list widget
        for i, cell_type in enumerate(
            self.dataset['rna'].obs['Cell_Type_Experimental'].unique()):
            item = QListWidgetItem(cell_type)
            color = PREDEFINED_COLORS[i]
            item.setForeground(QColor(*color))
            font = QFont()
            font.setBold(True)
            item.setFont(font)
            self.info_list_widget.addItem(item)
            self.labels_dict[cell_type] = self.dataset['rna'].obs[self.dataset['rna'].obs['Cell_Type_Experimental'] == cell_type].index

    def enable_dimred_buttons(self):
        self.compute_UMAP.setEnabled(True)
        self.compute_TSNE.setEnabled(True)

    def compute_umap(self):
        # This function will be called when the Compute UMAP button is clicked
        # Perform UMAP computation in a background thread
        selected_items = self.info_list_widget.selectedItems()
        # Extract the text (cell types/labels) of the selected items
        selected_cell_types = [item.text() for item in selected_items]
        # Get the indices of points corresponding to the selected cell types
        selected_indices = []
        for cell_type in selected_cell_types:
            if cell_type in self.labels_dict:
                selected_indices.extend(self.labels_dict[cell_type])

        runnable = UmapRunnable(
            self.dataset['rna'][selected_indices].copy(),
            self.color_map
        )
        self.threadpool.start(runnable)

    def compute_tsne(self):
        # This function will be called when the Compute t-SNE button is clicked
        # Perform t-SNE computation in a background thread
        runnable = TsneRunnable()
        self.threadpool.start(runnable)


app = QApplication(sys.argv)
w = MainWindow()
w.show()
app.exec()