from PyQt5 import QtCore, QtWidgets, uic
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as Canvas
from matplotlib.figure import Figure

# Matplotlib canvas class to create figure 
class MplCanvas(Canvas):
    def __init__(self):
        self.fig = Figure()
        #self.ax = self.fig.add_subplot(1,1,1)
        Canvas.__init__(self, self.fig)
        Canvas.setSizePolicy(self, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        #Canvas.updateGeometry(self)
            
# Matplotlib widget
class MplWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)   # Inherit from QWidget
        #Container Widget   
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.qlayout = QtWidgets.QHBoxLayout()
        self.qlayout.setSpacing(10)
        self.setLayout(self.qlayout)
        self.setContentsMargins(0,0,0,0)
        #Call Canvas Class
        self.canvas = MplCanvas()          
        self.canvas.setMinimumSize(self.canvas.sizeHint())
        #Scroll Area Properties     
        self.scroll=QtWidgets.QScrollArea(self)    
        self.scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.scroll.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.scroll.setWidgetResizable(False)
        self.scroll.setWidget(self.canvas)
        self.qlayout.addWidget(self.scroll)
        

       

        
  
