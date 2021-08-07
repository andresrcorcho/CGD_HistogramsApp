import matplotlib as mpl
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as Canvas
from matplotlib.offsetbox import AnchoredText
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter
from PyQt5 import QtCore, QtWidgets, uic, QtGui
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog, QMessageBox
from PyQt5.QtGui import QKeyEvent,QMouseEvent
from PyQt5.QtCore import QEvent, Qt
from PyQt5.QtTest import QTest
import numpy as np
import sys
import histograms
from functions import PDP, KDEp, KDE_PDP, loadData, kde_scipy, peakdet, kde_difussion,KDEadap_PDP
from adjustText import adjust_text
from pylab import figure, show, close
from io import StringIO
import ntpath
import pickle
import os
import pandas as pd
from pathlib import Path
from Cfunctions import GetisZeroOrString

#From FBS code
from fbs_runtime.application_context.PyQt5 import ApplicationContext
from PyQt5.QtWidgets import QMainWindow

#Application
class GeochronologyPlots(QtWidgets.QMainWindow, histograms.Ui_Geochronology):
    resized = QtCore.pyqtSignal()
    pressed=QtCore.pyqtSignal()
    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        #Application Icon
        self.setWindowIcon(QtGui.QIcon('icon.ico'))
        #Resized window event
        self.resized.connect(self.timerUpdate)
        #Clicked tab event
        self.InitProperties()
        self.Methods.currentChanged.connect(self.timerUpdate)
        #Setting Conections
        #Bw Sliders Setup
        #Initial Value
        self.BW.setText(str(self.Bw.value()))
        self.TR.setText(str(self.tratio.value()))
        #KDE Adaptative inactive by defect
        self.Methods.setTabEnabled(2,False)
        self.Methods.setTabEnabled(1,False)
        #Text Size
        self.tsize.setText(str(self.TSize.value()))
        #Events-Sliders
        #BW Slider
        self.Bw.sliderMoved.connect(self.updatelabel)
        self.Bw.sliderPressed.connect(self.updatelabel)
        self.Bw.valueChanged.connect(self.updatelabel)
        self.Bw.sliderReleased.connect(self.updatelabel)
        #Ticks Slider
        self.tratio.sliderMoved.connect(self.updateTicks)
        self.tratio.sliderPressed.connect(self.updateTicks)
        self.tratio.valueChanged.connect(self.updateTicks)
        self.tratio.sliderReleased.connect(self.updateTicks)
        #Peaks delta slider
        self.delta.sliderMoved.connect(self.updateDelta)
        self.delta.sliderPressed.connect(self.updateDelta)
        self.delta.valueChanged.connect(self.updateDelta)
        self.delta.sliderReleased.connect(self.updateDelta)
        #TextSize Slider
        self.TSize.sliderMoved.connect(self.updateTsize)
        self.TSize.sliderPressed.connect(self.updateTsize)
        self.TSize.valueChanged.connect(self.updateTsize)
        self.TSize.sliderReleased.connect(self.updateTsize)
        #Integrated plot Fuction
        self.plotData.clicked.connect(self.plotDensity)
        #Histogram option state
        self.Hist.clicked.connect(self.hist)
        #Detector slider state
        self.peakdetect.clicked.connect(self.peakDec)
        #Select File
        self.LoadData.clicked.connect(self.selectFiles)
        #Peak Numbers are shown?
        self.peakdetect.clicked.connect(self.peakLabels)
        #Load Status
        self.loadSession.clicked.connect(self.loadStatus)
        #Reset Status
        self.resetFields.clicked.connect(self.resetStatus)
        #Reset only files
        self.clearSlots.clicked.connect(self.clearStatus)
        #Save Status
        self.saveSession.clicked.connect(self.saveStatus)
        #PDP/KDE are Enabled?
        self.PDPstatus.clicked.connect(self.disableKDE)
        self.KDEstatus.clicked.connect(self.disableKDE)
        #Adjust of peak labels?
        self.peakLabel.clicked.connect(self.labelDetect)
        #shared mode protection
        #self.sharedXY.clicked.connect(self.shareEvent)
        #self.YAxisTicks.clicked.connect(self.shareEvent)
        self.Hist.clicked.connect(self.shareEvent)
        #self.geoScale.clicked.connect(self.shareEvent)
        self.flipPosition.clicked.connect(self.flipPositions)
        #self.DecimalX.clicked.connect(self.shareEvent)
        #self.DecimalY.clicked.connect(self.shareEvent)
        self.customBw.clicked.connect(self.customBwEv)
        #Expand Figure Mode
        #If Option is Checked
        self.exoandStatus.clicked.connect(self.expandMode)
        #Labels from peaks
        self.peakLabel.setEnabled(False)
        #Data Storage
        self.DataAges=[]
        self.DataErrors=[]
        self.Names=[]
        self.NSamples=[]
        self.BWs=[]
        self.indicesCounter=0
        self.initiated=False
        self.Nplots=0
        #Plot Status
        self.plotStatus=False
        #BW calculation
        self.B_ratio=float(20/50)
        self.Rbw=round(self.B_ratio*self.Bw.value(),1)
        self.BW.setText(str(self.Rbw))
        #Status of Loaded Files
        self.FileStatus=[]
        for i in range(1,18):
            self.FileStatus.append(False)
        #Files Counter
        self.NFiles=0
        #Initialize Tabs
        self.started()
        self.InitProperties()
    
    #After the app is initialized
    def started(self):
        self.timer = QtCore.QTimer()
        self.timer.setInterval(1)
        self.timer.start(0)
        self.timer.timeout.connect(self.InitProperties)
        
    def customBwEv(self):
        if self.customBw.isChecked()==True:
            self.Bw.setEnabled(False)
            self.AdaptativeBw.setEnabled(False)
            self.AdaptativeBw.setChecked(False)
        else:
            if self.AdaptativeBw.isChecked()==False:
                self.Bw.setEnabled(True)
                self.AdaptativeBw.setEnabled(True)
                self.customBw.setEnabled(True)
            else:
                self.customBw.setChecked(False)
                self.customBw.setEnabled(False)
                self.Bw.setEnabled(False)
    
    #If windows is resized a signal is emited
    def resizeEvent(self, event):
        self.resized.emit()
        return super(GeochronologyPlots, self).resizeEvent(event)
    
    def timerUpdate(self):
        self.timer = QtCore.QTimer()
        self.timer.setInterval(60)
        self.timer.start(0)
        self.timer.timeout.connect(self.UpdateProperties)
        
    def InitProperties(self):
        #Initialize
        maxW=self.getCurrentTabSize(self.KDEf)[0]-17
        maxH=self.getCurrentTabSize(self.KDEf)[1]-17
        self.KDEf.scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.PDPf.scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.KDEa.scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        #Get Maximum Widget Dimensions
        self.KDEf.scroll.setGeometry(QtCore.QRect(0, 0, maxW+17, maxH+17))
        self.PDPf.scroll.setGeometry(QtCore.QRect(0, 0, maxW+17, maxH+17))
        self.KDEa.scroll.setGeometry(QtCore.QRect(0, 0, maxW+17, maxH+17))
        self.KDEf.canvas.resize(maxW,maxH)
        self.PDPf.canvas.resize(maxW,maxH)
        self.KDEa.canvas.resize(maxW,maxH)
        
    def UpdateProperties(self):
        #Update Dimensions
        tabs=['KDEf','PDPf','KDEa']
        #Current Tab Dimensions
        counter=0
        #Get Current Widget Index
        index=self.Methods.currentIndex()
        
        #Check the Bw choice
        self.customBwEv()
        #Disables shared axes for plots of only one dataset
        indices=self.getIndices()
        if len(indices)<2:
            self.sharedXY.setEnabled(False)
            self.sharedXY.setChecked(False)
        else:
            self.sharedXY.setEnabled(True)
            
        #Disables the flip button when two datasets are selected
        if len(indices)==2:
            self.flipPosition.setEnabled(True)
        else:
            self.flipPosition.setEnabled(False)
            
        for i in tabs:
            base=eval('self.'+i)
            #Get Maximum scroll dimensions
            maxW=self.getCurrentTabSize(base)[0]-17
            maxH=self.getCurrentTabSize(base)[1]-17
            if self.KDEstatus.isChecked==False and counter==2:
                break
            elif self.PDPstatus.isChecked==False and counter==1:
                break
            if self.exoandStatus.isChecked()==False:
                #Update App events
                #Modifyng of Scroll bar dimensions
                base.scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
                #Get Maximum Widget Dimensions
                base.scroll.setGeometry(QtCore.QRect(0, 0, maxW+17, maxH+17))
                base.canvas.resize(maxW,maxH)
                #Update Size of Plots
                self.UpdateFigSize(base)
            else:
                #Update Scroll Size
                self.UpdateFigSize(base)
            #Adjust plot to the maximum plot area
            self.adjustToSize(base)
            #Counter
            counter+=1
    
    def getCurrentTabSize(self,Widg):
        width = float(Widg.frameGeometry().width())
        height =float(Widg.frameGeometry().height())
        return width,height
    
    def disableKDE(self):
        if self.KDEstatus.isChecked()==True:
            self.Methods.setTabEnabled(2,True)
            self.PDPstatus.setChecked(True)
        else:
            self.Methods.setTabEnabled(2,False)

        if self.PDPstatus.isChecked()==False:
            self.Methods.setTabEnabled(1,False)
            self.Methods.setTabEnabled(2,False)
        else:
            self.Methods.setTabEnabled(1,True)
    
    def UpdateFigSize(self,Widg):
        #Get Frame current Dimensions ,vertical,horizontal,Nsubplots
        width=self.getCurrentTabSize(Widg)[0]-17
        height =self.getCurrentTabSize(Widg)[1]-17
        #Size Factors from sliders
        if float(self.sizeFactorV.value())>1:
            factorV=(height/15.)*float(self.sizeFactorV.value())
        else:
            factorV=0
        if float(self.sizeFactorH.value())>1:
            factorH=(width/15.)*float(self.sizeFactorH.value())
        else:
            factorH=0
        #Changing FigSize Settings-If the Slider value is 1, the window size remains as the default screen area
        #Changing height
        height=height+factorV
        #Changing width
        width=width+factorH
        Widg.canvas.resize(width,height)
        #Widg.canvas.draw()
        Widg.canvas.updateGeometry()
         
    def expandMode(self):
        #Activation/Desactivation of Slider
        if self.exoandStatus.isChecked()==True:
            self.sizeFactorV.setEnabled(True)
            self.sizeFactorH.setEnabled(True)
            self.timerUpdate()
        else:
            self.sizeFactorV.setValue(1)
            self.sizeFactorH.setValue(1)
            self.sizeFactorV.setEnabled(False)
            self.sizeFactorH.setEnabled(False)
            self.timerUpdate()
    
    def adjustToSize(self,base):
        if self.sharedXY.isChecked()==True:
            base.canvas.fig.tight_layout(pad=0.1,h_pad=0.0,w_pad=0.0)
            base.canvas.fig.subplots_adjust(hspace=0.0)
            base.canvas.updateGeometry()
        else:
            base.canvas.fig.tight_layout(pad=0.5,h_pad=0.4,w_pad=0.0)
            base.canvas.updateGeometry()
            
    def shareEvent(self):
        app.app.processEvents()
        self.plotDensity()
        self.timerUpdate()
        
    def KDEdetect(self):
        if self.KDEstatus.isChecked()==True:
            self.Nplots+=1
        if self.PDPstatus.isChecked()==True:
            self.Nplots+=1
        self.UpdateProperties()
        
    def labelDetect(self):
        if self.peakLabel.isChecked()==True:
            self.adjustLabel.setEnabled(True)
        else:
            self.adjustLabel.setEnabled(False)
    
    def peakLabels(self):
        if self.peakdetect.isChecked()==True:
            self.peakLabel.setDisabled(False)
        else:
            self.peakLabel.setDisabled(True)
    
    def combineData(self,data):
        string=""
        for i in range(0,len(data)):
            aux=str(data[i])
            string=string+aux
        return string
    
    def updateTsize(self):
        self.tsize.setText(str(self.TSize.value()))
    
    def getIndices(self):
        indexes=[]
        for i in range (1,18):
            path='self.f'+str(i)+'.isChecked()==True'
            if eval(path):
                indexes.append(i-1)
        return indexes
    
    def updateDelta(self):
        self.peakvalue.setText(str(self.delta.value()))
    
    def selectFiles(self):
        #Get working directory
        userDir=os.path.expanduser('~/')
        currentDir=userDir+'CGD_DataStructure/'
        FilePaths=[]
        # open the dialog and get the selected files
        File, _filter = QtWidgets.QFileDialog.getOpenFileNames(self,'Load Data',currentDir+"/Data/","Data File (*.txt)")
            
        for j in File:
            try:
                self.selectFile(j)
            except StopIteration as why:
                pass
            
    def ClearCanvas(self,canva):
        base=eval("self."+canva+".canvas")
        base.fig.clear()
        base.draw()
    
    def CheckForErrors(self,ages,errors):
        Errors=np.array([])
        error=False
        err1="Zero Value Found in Ages Col!!!"
        err2="Missing or wrong character (text, symbol) Found in Ages Col!!!"
        err3="Zero Value Found in Errors Col!!!"
        err4="Missing or wrong character (text, symbol) Found in Errors Col!!!"
        
        #Look for errors in Ages Column
        error,ErrorsAges=GetisZeroOrString(ages,err1,err2)
#        for i in ages:
#            if i==0:
#                Errors.append(err1)
#                error=True
#            else:
#                try:
#                    aux=int(i)
#                except ValueError:
#                    Errors.append(err2)
#                    error=True
        #Look for errors in Uncertainties Column
        error,ErrorsErrors=GetisZeroOrString(errors,err3,err4)
#        for j in errors:
#            if j==0:
#                Errors.append(err3)
#                error=True
#            else:
#                try:
#                    aux=int(j)
#                except ValueError:
#                    Errors.append(err4)
#                    error=True
        #print(ErrorsAges)
        #print(ErrorsErrors)
        Errors=np.append(Errors,[ErrorsAges])
        Errors=np.append(Errors,[ErrorsErrors])
        #print (Errors)
        #return error, combineArray(Errors)
        return error, Errors
    
    def selectFile(self, File):
        
        if len(File)>0:
            self.initiated=True
        # if a file is selected
        if File:
            #The App have been initialized
#            ages=np.array([])
#            errors=np.array([])
#            self.initiated=True
#            data=loadData(File)
#            for i in range(0,len(data)):
#                age=data[i][0]
#                error=data[i][1]
#
#                #Add Data
#                ages=np.append(ages,[age])
#                errors=np.append(errors,[error])

            #Using Pandas approach
            data=loadData(File)
            ages = data[['Age']].to_numpy().flatten()
            errors =data[['Error']].to_numpy().flatten()
            #Verification of errors in file
            Verif=self.CheckForErrors(ages,errors)
            #print(Verif)
            
            #MsgBox With Found errors
            if Verif[0]==True:
                #String with all found errors
                errorsList=""
                for y in Verif[1]:
                    errorsList=errorsList+y+"<br><br>"
                    
                error_dialog = QtWidgets.QErrorMessage()
                error_dialog.setWindowIcon(QtGui.QIcon('icon.ico'))
                error_dialog.setWindowTitle('Error in loading File!!')
                error_dialog.showMessage("File "+(ntpath.basename(File))+" was not loaded because the following error/errors were found:<br><br>"+
                                        errorsList)
                
                error_dialog.exec_()
                return
            
            else:
            
                #Tittle by each file in the subplot
                
                with open(File) as fp:
                    lines=fp.readlines()
                    #print (lines[0])
                    #print (lines[1])
                    counter=0
                    for line in lines:
                        #print(fp)
                        if counter == 0:
                            lineName = line
                        elif counter == 1:
                            # 30th line
                            bwt=line
                        elif counter>1:
                            break
                        counter=counter+1
                fp.close()
                
                self.Names.append(lineName)
                #Storage of Ages
                self.DataAges.append(ages)
                #Storage of Number of samples in file
                self.NSamples.append(int(np.size(ages)))
                #Storage of uncertainties
                self.DataErrors.append(errors)
                #Storage of number of loaded files
                self.indicesCounter=self.indicesCounter+1
                #Add File Number to slots
                File=ntpath.basename(File)
                #Update Number of loaded Files
                self.NFiles=self.NFiles+1
                #Add to bandwidths
                self.BWs.append(float(bwt))
            # update the field with a new set of data
                for i in range(1,18):
                    root=eval('self.f'+str(i))
                    box=eval('self.f'+str(i)+'.isChecked()==False and self.FileStatus[i-1]==False')
                    if box:
                        root.setEnabled(True)
                        root.setCheckState(True)
                        root.setText(File)
                        self.FileStatus[i-1]=True
                        break
    
    #Get the informacion of the current session
    def getStatus(self):
        #Array wich Stores all Status information
        Status=[]
        #Interface Seetting Status Variables
        FilePathNames=np.array([]) #Directory of Each File
        SlotCheckState=np.array([]) #State of loaded files, checked or not checked?
        #Array which Stores Settings Variables
        VariablesStatus=['self.dataName.text()',
                         'self.Bw.value()',
                         'self.tratio.value()',
                         'float(self.Minin.text())',
                         'float(self.Maxi.text())',
                         'self.savepdf.isChecked()',
                         'self.savepng.isChecked()',
                         'self.Hist.isChecked()',
                         'float(self.bins.text())',
                         'self.peakdetect.isChecked()',
                         'self.delta.value()',
                         'self.KDEstatus.isChecked()',
                         'self.filled.isChecked()',
                         'float(self.TSize.value())',
                         'self.sharedXY.isChecked()',
                         'self.exoandStatus.isChecked()',
                         'self.PDPstatus.isChecked()',
                         'float(self.sizeFactorV.value())',
                         'float(self.sizeFactorH.value())']
       
        #Get Interface Settings Status
        for i in VariablesStatus:
            Status.append(eval(i))
        
        #Get Slots Status and filenames
        for k in range(1,18):
            root=eval('self.f'+str(k))
            FilePathNames=np.append(FilePathNames,[root.text()])
            SlotCheckState=np.append(SlotCheckState,[root.isChecked()])
            
        #Include Slots values and Filenames in status
        Status.append(FilePathNames)
        Status.append(SlotCheckState)
        return Status
        
    #Load and displays a custom session-Data is loaded from a Collection of Data- (List of arrays)
    def setStatus(self,Status):
        self.resetStatus()
        #Local Variables to Store data
        Aux=[]
        FileDir=np.array([])
        CheckSlots=np.array([])
        #Get data from Status dataset
        FileDir=Status[19][:]
        CheckSlots=Status[20][:]
        #Get Status Variables from Status dataset
        counter=0
        for j in Status:
            Aux.append(j)
            counter=counter+1
            if counter>18:
                break
        
        #Set status as loaded status File
        StatusVariables=['self.dataName.setText(Aux[0])',
                         'self.Bw.setValue(Aux[1])',
                         'self.tratio.setValue(Aux[2])',
                         'self.Minin.setText(str(Aux[3]))',
                         'self.Maxi.setText(str(Aux[4]))',
                         'self.savepdf.setChecked(Aux[5])',
                         'self.savepng.setChecked(Aux[6])',
                         'self.Hist.setChecked(Aux[7])',
                         'self.bins.setText(str(Aux[8]))',
                         'self.peakdetect.setChecked(Aux[9])',
                         'self.delta.setValue(Aux[10])',
                         'self.KDEstatus.setChecked(Aux[11])',
                         'self.filled.setChecked(Aux[12])',
                         'self.TSize.setValue(Aux[13])',
                         'self.sharedXY.setChecked(Aux[14])',
                         'self.exoandStatus.setChecked(Aux[15])',
                         'self.PDPstatus.setChecked(Aux[16])',
                         'self.sizeFactorV.setValue(Aux[17])',
                         'self.sizeFactorH.setValue(Aux[18])']
        
        #Change to a new status given as a parameter
        for i in StatusVariables:
            eval(i)
            
        #Disable/Enable optional settings
        #Histogram
        if self.Hist.isChecked()==True:
            self.bins.setEnabled(True)
        else:
            self.bins.setEnabled(False)
        
        #Peak Detector
        if self.peakdetect.isChecked()==True:
            self.delta.setEnabled(True)
            self.peakLabel.setEnabled(True)
        else:
            self.delta.setEnabled(False)
        
        #Look for changes in tabs state
        self.disableKDE()

        #Load Files
        for j in range(1,len(FileDir)):
            if FileDir[j-1]=="":
                break
            else:
                root=eval('self.f'+str(j))
                currentDir=userDir+'CGD_DataStructure/'
                file=currentDir+'/Data/'+FileDir[j-1]
                try:
                    self.selectFile(file)
                    root.setChecked(CheckSlots[j-1])
                except TypeError as why:
                    pass
                
        #Look for changes in expand mode
        self.expandMode()
            
    def resetStatus(self):
        #Reset Canvas
        self.ClearCanvas("KDEf")
        self.ClearCanvas("PDPf")
        self.ClearCanvas("KDEa")
        #Reset Data Variables
        self.DataAges=[]
        self.DataErrors=[]
        self.Names=[]
        self.NSamples=[]
        self.BWs=[]
        self.indicesCounter=0
        self.initiated=False
        self.Nplots=0
        #Status of Loaded Files
        self.FileStatus=[]
        for i in range(1,18):
            root=eval('self.f'+str(i))
            self.FileStatus.append(False)
            root.setEnabled(False)
            root.setChecked(False)
            root.setText("")
        #Files Counter
        self.NFiles=0
        #Reset Status Variables
        self.dataName.setText("Unnamed")
        self.Bw.setValue(7)
        self.tratio.setValue(2)
        self.Minin.setText("0")
        self.Maxi.setText("2000")
        self.savepdf.setChecked(False)
        self.savepng.setChecked(False)
        self.Hist.setChecked(False)
        self.bins.setText("100")
        self.peakdetect.setChecked(False)
        self.delta.setValue(2)
        self.KDEstatus.setChecked(False)
        self.Methods.setTabEnabled(2,False)
        self.Methods.setTabEnabled(1,False)
        self.exoandStatus.setChecked(False)
        self.sizeFactorV.setEnabled(False)
        self.filled.setChecked(True)
        self.TSize.setValue(8)
        self.sharedXY.setChecked(False)
        self.peakLabel.setDisabled(False)
        self.bins.setEnabled(False)
        self.delta.setEnabled(False)
        self.peakLabel.setEnabled(False)
        #Look for changes in tabs status
        self.disableKDE()
        #Update Figure Settings
        self.UpdateProperties()
        
    def clearStatus(self):
        #Reset Canvas
        self.ClearCanvas("KDEf")
        self.ClearCanvas("PDPf")
        self.ClearCanvas("KDEa")
        #Reset Data Variables
        self.DataAges=[]
        self.DataErrors=[]
        self.Names=[]
        self.NSamples=[]
        self.indicesCounter=0
        self.initiated=False
        self.Nplots=0
        #Status of Loaded Files
        self.FileStatus=[]
        for i in range(1,18):
            root=eval('self.f'+str(i))
            self.FileStatus.append(False)
            root.setEnabled(False)
            root.setChecked(False)
            root.setText("")
        #Files Counter
        self.NFiles=0
    
    def saveStatus(self):
        #Get working directory
        userDir=os.path.expanduser('~/')
        currentDir=userDir+'CGD_DataStructure/'
        #Get current Status and ask for File Name
        Status=self.getStatus()
        try:
            save_dialog=QtWidgets.QFileDialog()
            fileName, _filter = save_dialog.getSaveFileName(self,'Save Dataset',currentDir+"/Datasets/","Status File (*.p)")
            #Save Data in binary file
            pickle.dump(Status, open(fileName, "wb" ) )
            
        except FileNotFoundError as why:
            
            pass

    def loadStatus(self):
        userDir=os.path.expanduser('~/')
        currentDir=userDir+'CGD_DataStructure/'
        #Load Data from binary format
        try:
            fileName, _filter = QtWidgets.QFileDialog.getOpenFileName(self,'Load Dataset',currentDir+"/Datasets/","Status File (*.p)")
            Status= pickle.load(open(fileName, "rb" ) )
            self.setStatus(Status)
            
        except FileNotFoundError as why:
            
            pass
        
    def updateTicks(self):
        self.TR.setText(str(self.tratio.value()))
    
    def updatelabel(self):
        self.Rbw=round(self.B_ratio*self.Bw.value(),1)
        self.BW.setText(str(self.Rbw))
    
    def hist(self):
        if self.Hist.isChecked()==True:
            self.bins.setDisabled(False)
        else:
            self.bins.setDisabled(True)
            
    def peakDec(self):
        if self.peakdetect.isChecked()==True:
            self.delta.setDisabled(False)
        else:
            self.delta.setDisabled(True)
    
    def clearFigures(self):
        Plots=['KDEf','PDPf','KDEa']
        for i in Plots:
            base=eval('self.'+i+'.canvas')
            base.fig.clear()
            
    def arraytoScalar(self,data):
        x=np.array([])
        y=np.array([])
    
        for i in range (0,len(data)):
            auxX=round(data[:,0][i],1)
            auxY=round(data[:,1][i],1)
            x=np.append(x,[auxX])
            y=np.append(y,[auxY])
        
        return x,y
            
    def savePlot(self, base,Tabcounter):
        userDir=os.path.expanduser('~/')
        currentDir=userDir+'CGD_DataStructure/'
        #Saving Data
        PlotP=["KDE-Fixed_"+self.dataName.text()+".png","PDP_"+self.dataName.text()+".png","KDE-PDP_"+self.dataName.text()+".png"]
        PlotF=["KDE-Fixed_"+self.dataName.text()+".pdf","PDP_"+self.dataName.text()+".pdf","KDE-PDP_"+self.dataName.text()+".pdf"]
        macWindow=['\\Results\\','/Results/']
        #Path to Folder with Data
        data_folder=(currentDir+"/Results/")
        
        file=Path(data_folder+PlotP[Tabcounter])
        file1=Path(data_folder+PlotF[Tabcounter])
        if self.savepng.isChecked()==True:
            try:
                base.fig.savefig(file)
            except PermissionError as why:
                pass

        if self.savepdf.isChecked()==True:
            try:
                base.fig.savefig(file1)
            except PermissionError as why:
                pass
                
        
    def flipPositions(self):
        self.flipPosition.setEnabled(False)
        cIndxs=self.getIndices()
        #Upper
        PosU=cIndxs[0]
        Path1=eval('self.f'+str(PosU+1))
        #Lower
        PosL=cIndxs[1]
        Path2=eval('self.f'+str(PosL+1))
        #Store data temporaly
        dataAux1=[self.Names[PosU],self.NSamples[PosU],self.DataAges[PosU],self.DataErrors[PosU],self.BWs[PosU]]
        dataAux2=[self.Names[PosL],self.NSamples[PosL],self.DataAges[PosL],self.DataErrors[PosL],self.BWs[PosL]]
        
        #Postions
        cpositions=[PosU,PosL]
        #Fliped positions
        cFpositions=[PosL,PosU]
        #Fliped Data
        flipedData=[dataAux2,dataAux1]
        
        for i,fliped in zip (cpositions, flipedData):
        
            self.Names[i]=fliped[0]
            #Storage of Number of samples in file
            self.NSamples[i]=fliped[1]
            self.DataAges[i]=fliped[2]
            #Storage of uncertainties
            self.DataErrors[i]=fliped[3]
            self.BWs[i]=fliped[4]
            
        Text1=Path1.text()
        Text2=Path2.text()
        
        Path1.setText(Text2)
        Path2.setText(Text1)
        
        #self.shareEvent()
        self.flipPosition.setEnabled(True)
        
    def getLocalNsamples(self,Ages):
        Min=float(self.Minin.text())
        Max=float(self.Maxi.text())
        NumAges=0
        for age in Ages:
            if age>=Min and age<=Max:
                NumAges=NumAges+1
        return NumAges
    
    
    def plotDensity(self):
        #Lock Plot button for prevent duplication of Axex
        self.plotData.setEnabled(False)
        #Plot Status
        #self.plotStatus=True
        #Reset Canvas
        self.Nplots=0
        self.clearFigures()
        if self.initiated==True:
            #Reset widgets and state variables
            self.clearFigures()
            #Vector with indices with data to evaluate
            indices=self.getIndices()
            #Look for Data to plot
            self.KDEdetect()
            #Evaluate Data
            Min=float(self.Minin.text())
            Max=float(self.Maxi.text())
            #Vector to Evaluate data
            x1_grid=np.linspace(Min,Max,2048)
            Plots=['KDEf','PDPf','KDEa']
            
            if self.AdaptativeBw.isChecked()==True:
                calculationKDE='kde_difussion(x1_grid,self.DataAges[i],Min,Max)[0]'
                calculationKDEa_PDP='KDEadap_PDP(self.DataAges[i],self.DataErrors[i],x1_grid,Min,Max)'
            else:
                calculationKDE='KDEp(x1_grid,self.DataAges[i],bandwidth=auxBw)'
                calculationKDEa_PDP='KDE_PDP(self.DataAges[i],self.DataErrors[i],x1_grid,auxBw)'
            
            Calculations=[calculationKDE,'PDP(x1_grid,self.DataAges[i],self.DataErrors[i])',calculationKDEa_PDP]
            Tabcounter=0
            #Plot Bucle
            for i in Plots:
                if Tabcounter==1 and self.KDEstatus.isChecked()==False and self.PDPstatus.isChecked()==False:
                    break
                if Tabcounter==2 and int(self.Nplots)<2:
                    break
                base=eval('self.'+i+'.canvas')
                app.app.processEvents()
                self.Methods.setCurrentWidget(base)
                self.Methods.setCurrentIndex(Tabcounter)
                #Plot Parameters
                self.plotCanvas(Min,Max,indices,base,x1_grid,Calculations[Tabcounter],Tabcounter)
                #Initializing Figure
                app.app.processEvents()
                self.timerUpdate()
                app.app.processEvents() #El UNICO QUE SE NECESITA
                base.draw()
                self.savePlot(base,Tabcounter)
                #Counter of tabs to plot
                Tabcounter=Tabcounter+1
        #Enable Plotting after completing the plot process
        self.plotData.setEnabled(True)
                        
    def plotCanvas(self,Min,Max,indices,base,x1_grid,Calculation,counter):
    #Editable text
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
        #Plot Counter
        plotCounter=0
        #Initializing Tabs for plots
        #Plot Limits
        if self.sharedXY.isChecked()==True:
            #Common X,Y Axes
            bx=base.fig.add_subplot(1,1,1)
            # Turn off axis lines and ticks of the big subplot
            bx.spines['top'].set_color('none')
            bx.spines['bottom'].set_color('none')
            bx.spines['left'].set_color('none')
            bx.spines['right'].set_color('none')
            #bx.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
            bx.tick_params(labelsize=float(self.TSize.value()))
            if self.DecimalY.isChecked()==False:
                bx.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            else:
                bx.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            bx.get_xaxis().set_ticks([])
                
        for i in indices:
            #Use local or global Bandwidth
            if self.customBw.isChecked()==True:
                auxBw=self.BWs[i]
            elif self.customBw.isChecked()==False:
                if self.AdaptativeBw.isChecked()==False:
                    auxBw=self.Rbw
                else:
                    try:
                        auxBw=np.round(eval('kde_difussion(x1_grid,self.DataAges[i],Min,Max)[1]'),2)
                    except ValueError:
                        #If the KDE adaptative calculation fails it will display the message and use the default Bw
                        KDEerror_dialog = QtWidgets.QErrorMessage()
                        KDEerror_dialog.setWindowIcon(QtGui.QIcon('icon.ico'))
                        KDEerror_dialog.setWindowTitle('KDE Error')
                        KDEerror_dialog.showMessage("Bandwidth optimization did not converge.The bandwidth for the: <br><br>"+
                                            self.Names[i]+ "dataset"+" will be set to default")
                                            
                        KDEerror_dialog.exec_()
                        #Set bandwidth to Default
                        auxBw=self.Rbw
                        if Calculation=='kde_difussion(x1_grid,self.DataAges[i],Min,Max)[0]':
                            Calculation='KDEp(x1_grid,self.DataAges[i],bandwidth=auxBw)'
                        elif Calculation=='KDEadap_PDP(self.DataAges[i],self.DataErrors[i],x1_grid,Min,Max)':
                            Calculation='KDE_PDP(self.DataAges[i],self.DataErrors[i],x1_grid,auxBw)'
                        #pass
            #Update Bandwidth label
            #self.updatelabel()
            #Update Number of samples within the time interval specified by the user
            self.NSamples[i]=self.getLocalNsamples(self.DataAges[i])
            #KDE Calculation --
            KDE=eval(Calculation)
            #Subplots Inicialization
            ax=base.fig.add_subplot(len(indices),1,plotCounter+1)
            #Bins (Histogram) with Data
            if self.Hist.isChecked()==True:
                ax1=ax.twinx()
                
                #Filter Data Ages for Histogram
                DataAgesHist=np.array([])
                for age in self.DataAges[i]:
                    if age>=Min and age<=Max:
                        DataAgesHist=np.append(DataAgesHist,[age])
                
                y, x, _ =ax.hist(DataAgesHist,bins=int(float(self.bins.text())),color='cornflowerblue',edgecolor='darkgreen', linewidth=0.8, density=False)
                ax.yaxis.set_major_locator(mpl.ticker.LinearLocator(8))
                #Verification of Shared Axes
                if self.sharedXY.isChecked()==True:
                    ax.set_xlabel("Age (Ma)", fontsize=float(self.TSize.value()))
                    bx.get_xaxis().set_ticks([])
                    anchored_title=AnchoredText(self.Names[i],loc='upper center',pad=0.1,borderpad=0.1,frameon=False,prop=dict(size=float(self.TSize.value())*1.2))
                    ax1.add_artist(anchored_title)

                else:
                    #X and Y axis-labels per histogram
                    ax.set_ylabel("Frequency", fontsize=float(self.TSize.value()))
                    ax.set_xlabel("Age (Ma)", fontsize=float(self.TSize.value()))
                    #ax1.yaxis.set_label_position("left")
                    #ax.yaxis.set_label_position("right")
                    #ax.yaxis.tick_right()
                    #ax.yaxis.set_label_position("right")
                    #Title from loadaed file
                    anchored_title=AnchoredText(self.Names[i],loc='upper center',pad=0.1,borderpad=0.1,frameon=False,prop=dict(size=float(self.TSize.value())*1.2))
                    ax1.add_artist(anchored_title)
                            
                ax1.set_xlim([Min, Max])
                ax.set_xlim([Min, Max])
                
                if self.Methods.currentIndex()<2:
                    ax1.set_ylim([0,max(KDE)+(max(KDE)/4)])
                    ax.set_ylim([0,y.max()+(y.max()/4)])
                else:
                    auMax=(KDE[0])
                    if max(KDE[1])>max(auMax):
                        auMax=KDE[1]
                    ax1.set_ylim([0,max(auMax)+(max(auMax)/4)])
                    ax.set_ylim([0,y.max()+(y.max()/4)])
                #Plotting density Function
                if self.Methods.currentIndex()<1:
                    ax1.plot(x1_grid,KDE,label='n={0} - Bw={1}'.format((eval(str(self.NSamples[i]))),eval(str(auxBw))), color='mediumblue',linewidth=0.8)
                    if self.filled.isChecked()==True:
                        ax1.fill_between(x1_grid,KDE,color='powderblue', alpha=.5)
                elif self.Methods.currentIndex()==1:
                    ax1.plot(x1_grid,KDE,label='n={0}'.format((eval(str(self.NSamples[i]))),eval(str(auxBw))), color='mediumblue',linewidth=0.8)
                else:
                    ax1.plot(x1_grid,KDE[0],label='KDE n={0} - Bw={1}'.format((eval(str(self.NSamples[i]))),eval(str(auxBw))), color='mediumblue',linewidth=0.8)
                    ax1.plot(x1_grid,KDE[1],label='PDP n={0}'.format((eval(str(self.NSamples[i])))), color='saddlebrown',linewidth=0.8)
                
                #Axis settings
                ax1.yaxis.set_major_locator(mpl.ticker.LinearLocator(8))
                ax1.legend(loc='upper right',fontsize=int(self.TSize.value()))
                ax1.get_yaxis().set_ticks([])
                
                if self.peakdetect.isChecked()==True and self.Methods.currentIndex()<2:
                    ax2=ax1.twinx()
                    #Delta Settings
                    maxvalue=((max(KDE)-min(KDE))/2)
                    minvalue=max(KDE)-(maxvalue*2)
                    if self.Methods.currentIndex()==0:
                        ratio=(maxvalue)/200
                    elif self.Methods.currentIndex()==1:
                        ratio=ratio=((maxvalue)/200)*5
                    #Peaks Detection
                    maxtab, mintab = peakdet(KDE,(float(self.delta.value())*(ratio)),x1_grid)
                    #Plotting Peaks
                    texts=[]
                    Peaks=ax2.scatter(maxtab[:,0], maxtab[:,1], color="red")
                    if self.peakLabel.isChecked()==True:
                        if self.adjustLabel.isChecked()==False:
                            for f in range(0,len(maxtab[:,0])):
                                if self.DecimalPeaks.isChecked()==True:
                                    anotateKey=(np.round(maxtab[:,0][f],1))
                                else:
                                    anotateKey=int((maxtab[:,0][f]))
                                ax2.annotate(str(anotateKey),(maxtab[:,0][f], maxtab[:,1][f]), size=int(self.TSize.value()))
                        else:
                            texts=[]
                            for txt,x,y in zip(maxtab[:,0], maxtab[:,0],  maxtab[:,1]):
                                if self.DecimalPeaks.isChecked()==True:
                                    txt=str(np.round(txt,1))
                                else:
                                    txt=str(int(txt))
                                texts.append(ax2.text(x, y,txt,size=int(self.TSize.value())))
                            if self.sharedXY.isChecked()==True:
                                adjust_text(texts,maxtab[:,0],maxtab[:,1],ax=ax2,expand_text=(1.05, 1),autoalign='xy',expand_points=(1.01, 1.05),
            force_text=(0.01, 0.25), force_points=(0.01, 0.25))
                            else:
                                adjust_text(texts,maxtab[:,0],maxtab[:,1],ax=ax2,expand_text=(1.05, 1),autoalign='xy',expand_points=(1.01, 1.05),
            force_text=(0.01, 0.25), force_points=(0.01, 0.25))
                    ax2.get_yaxis().set_ticks([])
                    ax2.set_xlim([Min, Max])
                    if self.Methods.currentIndex()<2:
                        ax2.set_ylim([0,max(KDE)+(max(KDE)/4)])
                    else:
                        auMax=(KDE[0])
                        if max(KDE[1])>max(auMax):
                            auMax=KDE[1]
                        
                        ax2.set_ylim([0,max(auMax)+(max(auMax)/4)])
            else:
                #information about axis -from user
                ax3=ax.twinx()
                #Filter Data Ages for Histogram
                DataAgesHist=np.array([])
                for age in self.DataAges[i]:
                    if age>=Min and age<=Max:
                        DataAgesHist=np.append(DataAgesHist,[age])
                ax3.yaxis.set_major_locator(mpl.ticker.LinearLocator(8))
#               #Histogram parameters for getting number of samples
                hst,binss = np.histogram(DataAgesHist, bins=int(len(x1_grid))) #Fixed and based in x1_grid length
                #hst,binss = np.histogram(DataAgesHist, bins=int(float(self.bins.text()))) - Version March 2021
                tickss=np.linspace(0,hst.max()+(hst.max()/4),8)
                #Get Histogram ticks
            
                #print(tickss)
                ax3.yaxis.tick_left()
                ax3.yaxis.set_label_position("left")
                NSamples_ticks = []
                if self.Hist.isChecked()==False:

                    #label_format = '%.1f'
                    counter=0
                    for tick in ax3.get_yticks():
                        tick = tickss[counter]
                        if self.DecimalY.isChecked()==True:
                            NSamples_ticks.append('%.1f' % (tick,))
                        else:
                            NSamples_ticks.append('%.0f' % (tick,))
                        
                        counter=counter+1
                ax3.set_yticklabels(NSamples_ticks)
                
            
                #Verification of Shared Axes
                if self.sharedXY.isChecked()==True:
                    ax.set_xlabel("Age (Ma)", fontsize=float(self.TSize.value()))
                    bx.get_xaxis().set_ticks([])
                    anchored_title=AnchoredText(self.Names[i],loc='upper center',pad=0.1,borderpad=0.1,frameon=False,prop=dict(size=float(self.TSize.value())*1.2))
                    ax.add_artist(anchored_title)
                else:
                    #X and Y axis-labels per histogram
                    ax3.set_ylabel("Frequency",fontsize=int(self.TSize.value()))
                    ax.set_xlabel("Age (Ma)",fontsize=int(self.TSize.value()))
                    #Title from loadaed file
                    anchored_title=AnchoredText(self.Names[i],loc='upper center',pad=0.1,borderpad=0.1,frameon=False,prop=dict(size=float(self.TSize.value())*1.2))
                    ax.add_artist(anchored_title)
                            
                ax.set_xlim([Min, Max])
                if self.Methods.currentIndex()<2:
                    ax.set_ylim([0,max(KDE)+(max(KDE)/4)])
                else:
                    auMax=(KDE[0])
                    if max(KDE[1])>max(auMax):
                        auMax=KDE[1]
                    ax.set_ylim([0,max(auMax)+(max(auMax)/4)])
                
                #Plotting density Function
                if self.Methods.currentIndex()<1:
                    ax.plot(x1_grid,KDE,label='n={0} - Bw={1}'.format((eval(str(self.NSamples[i]))),eval(str(auxBw))), color='mediumblue',linewidth=0.8)
                    if self.filled.isChecked()==True:
                        ax.fill_between(x1_grid,KDE,color='powderblue')
                elif self.Methods.currentIndex()==1:
                    ax.plot(x1_grid,KDE,label='n={0}'.format((eval(str(self.NSamples[i]))),eval(str(auxBw))), color='mediumblue',linewidth=0.8)
                else:
                    ax.plot(x1_grid,KDE[0],label='KDE n={0} - Bw={1}'.format((eval(str(self.NSamples[i]))),eval(str(auxBw))), color='mediumblue',linewidth=0.8)
                    ax.plot(x1_grid,KDE[1],label='PDP n={0}'.format((eval(str(self.NSamples[i])))), color='saddlebrown',linewidth=0.8)
                
                #Axis settings
                ax.yaxis.set_major_locator(mpl.ticker.LinearLocator(8))
                #Change ticks from histogram ones
                ax.legend(loc='upper right',fontsize=int(self.TSize.value()))
                ax.get_yaxis().set_ticks([])
                
                #Peaks Detection
                if self.peakdetect.isChecked()==True and self.Methods.currentIndex()<2:
                    ax2=ax.twinx()
                    #Delta Settings
                    maxvalue=((max(KDE)-min(KDE))/2)
                    minvalue=max(KDE)-(maxvalue*2)
                    if self.Methods.currentIndex()==0:
                        ratio=(maxvalue)/200
                    elif self.Methods.currentIndex()==1:
                        ratio=((maxvalue)/200)*5
                    #Peaks Detection
                    maxtab, mintab = peakdet(KDE,(float(self.delta.value())*(ratio)),x1_grid)
                    #Plotting Peaks
                    texts=[]
                    Peaks=ax2.scatter(maxtab[:,0], maxtab[:,1], color="red")
                    if self.peakLabel.isChecked()==True:
                        if self.adjustLabel.isChecked()==False:
                            for f in range(0,len(maxtab[:,0])):
                                if self.DecimalPeaks.isChecked()==True:
                                    anotateKey=(np.round(maxtab[:,0][f],1))
                                else:
                                    anotateKey=int((maxtab[:,0][f]))
                                ax2.annotate(str(anotateKey),(maxtab[:,0][f], maxtab[:,1][f]), size=int(self.TSize.value()))
                        else:
                            texts=[]
                            for txt,x,y in zip(maxtab[:,0], maxtab[:,0],  maxtab[:,1]):
                                if self.DecimalPeaks.isChecked()==True:
                                    txt=str(np.round(txt,1))
                                    
                                else:
                                    txt=str(int(txt))
                                #txt=str(int(txt))
                                texts.append(ax2.text(x, y,txt,size=int(self.TSize.value())))
                            if self.sharedXY.isChecked()==True:
                                adjust_text(texts,maxtab[:,0],maxtab[:,1],ax=ax2,expand_text=(1.05, 1),autoalign='xy',expand_points=(1.01, 1.05),
            force_text=(0.01, 0.25), force_points=(0.01, 0.25))
                            else:
                                adjust_text(texts,maxtab[:,0],maxtab[:,1],ax=ax2,expand_text=(1.05, 1),autoalign='xy',expand_points=(1.01, 1.05),
            force_text=(0.01, 0.25), force_points=(0.01, 0.25))
                    ax2.get_yaxis().set_ticks([])
                    ax2.set_xlim([Min, Max])
                    if self.Methods.currentIndex()<2:
                        ax2.set_ylim([0,max(KDE)+(max(KDE)/4)])
                    else:
                        auMax=(KDE[0])
                        if max(KDE[1])>max(auMax):
                            auMax=KDE[1]
                        ax2.set_ylim([0,max(auMax)+(max(auMax)/4)])
                   
            #Ticks Adjust
            adjust=(Max-Min)/(float(self.tratio.value()))
            #Arrange ticks Array
            arrangeTicks=np.arange(Min,Max,adjust)
            #Ticks Change formatting
            if self.Hist.isChecked()==False:
                ax.xaxis.set_ticks(arrangeTicks)
                #ax.get_yaxis().set_ticks([])
                #Ticks Float formatting
                if self.DecimalX.isChecked()==False:
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
                else:
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    
                if self.DecimalY.isChecked()==False:
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
                else:
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            else:
                ax1.xaxis.set_ticks(arrangeTicks)
                #Ticks Float formatting
                if self.DecimalX.isChecked()==False:
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
                else:
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    
                if self.DecimalY.isChecked()==False:
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
                else:
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            #Integers???
            #from matplotlib.ticker import MaxNLocator
            #ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            #X-Axis Ticks Modification
            counter=0
            reset=0
            if self.geoScale.isChecked()==True:
                #Colors Geological timescale
                #Pliocene-Present
                ax.axvspan(0, 5.333, facecolor='#FFF9AE', alpha=0.5)
                #Miocene
                ax.axvspan(5.333, 23.03, facecolor='#FFF200', alpha=0.5)
                #Oligocene
                ax.axvspan(23.03, 33.9, facecolor='#FFDAAB', alpha=0.5)
                #Eocene
                ax.axvspan(33.9, 56, facecolor='#FCBC86', alpha=0.5)
                #Paleocene
                ax.axvspan(56, 66, facecolor='#FAB27B', alpha=0.5)
                #Neogene
                #ax.axvspan(0, 23.03, facecolor='#FFDD00', alpha=0.5)
                #Paleogene
                #ax.axvspan(23.03, 66, facecolor='#F9A870', alpha=0.5)
                #Cretaceous divisions
                ax.axvspan(66, 100.5, facecolor='#C7E09B', alpha=0.5)
                ax.axvspan(100.5, 145, facecolor='#94CC78', alpha=0.5)
                #Jurassic
                ax.axvspan(145, 163.5, facecolor='#ABE1FA', alpha=0.5)
                ax.axvspan(163.5, 174.1, facecolor='#71CFEB', alpha=0.5)
                ax.axvspan(174.1, 201.3, facecolor='#00B4EB', alpha=0.5)
                #ax.axvspan(145, 201.3, facecolor='#00B9E7', alpha=0.5)
                #Triassic
                ax.axvspan(201.3, 251.9    , facecolor='#8F53A1', alpha=0.5)
                #Permian
                ax.axvspan(251.9, 298.9, facecolor='#E8644B', alpha=0.5)
                #Carboniferous
                ax.axvspan(298.9, 358.9, facecolor='#68AEB2', alpha=0.5)
                #Devonian
                ax.axvspan(358.9,419.2, facecolor='#CF9C5A', alpha=0.5)
                #Silurian
                ax.axvspan(419.2, 443.8, facecolor='#B3DDCA', alpha=0.5)
                #Ordovician
                ax.axvspan(443.8, 485.4, facecolor='#00A88E', alpha=0.5)
                #Cambrian
                ax.axvspan(485.4, 541, facecolor='#8CAB79', alpha=0.5)
                #Middle_early Paleozoic
                #ax.axvspan(298.9, 541, facecolor='#9EC1A6', alpha=0.5)
                #Proterozoic
                ax.axvspan(541,1000, facecolor='#FBBA63', alpha=0.5)
                ax.axvspan(1000, 1600, facecolor='#FBBB7E', alpha=0.5)
                ax.axvspan(1600, 2500, facecolor='#F06682', alpha=0.5)
                #Archean
                ax.axvspan(2500, 4000, facecolor='#ED2891', alpha=0.5)
            
            #Arrange ticks in the x-axis
            for i in arrangeTicks:
                if reset!=0:
                    ax.xticks = ax.xaxis.get_major_ticks()
                    ax.xticks[counter].label1.set_visible(False)
                counter=counter+1
                reset=reset+1
                if reset==2:
                    reset=0
                    
            #Remove label of last tick in the y-axis
#
            if self.sharedXY.isChecked()==True and self.Hist.isChecked()==False:
                ax3.yticks = ax3.yaxis.get_major_ticks()
                ax3.yticks[-1].label1.set_visible(False)
            elif self.sharedXY.isChecked()==True and self.Hist.isChecked()==True:
                ax.yticks = ax.yaxis.get_major_ticks()
                ax.yticks[-1].label1.set_visible(False)
                    
            #Ticks Font-Size
            
            if self.Hist.isChecked()==False:
                ax3.tick_params(axis = 'both', which = 'major', labelsize = float(self.TSize.value()))
                ax.tick_params(axis = 'both', which = 'major', labelsize = float(self.TSize.value()))
            else:
                ax1.tick_params(axis = 'both', which = 'major', labelsize = float(self.TSize.value()))
                ax.tick_params(axis = 'both', which = 'major', labelsize = float(self.TSize.value()))
                    
            #Check for Shared XY option- Much more plot area
            #print(len(indices),plotCounter)
            if self.sharedXY.isChecked()==True:
                if plotCounter != len(indices)-1 and self.Hist.isChecked()==True:
                    ax1.get_xaxis().set_visible(False)
                    ax1.set_xticklabels([])
                elif plotCounter != len(indices)-1 and self.Hist.isChecked()==False:
                    #x-axis
                    ax.get_xaxis().set_visible(False)
                    ax.set_xticklabels([])
                
                #y-axis
                    
            if self.sharedXY.isChecked()==True:
                
                bx.set_ylabel("Frequency", labelpad=0, fontsize=float(self.TSize.value()))
                bx.yaxis.set_major_locator(mpl.ticker.LinearLocator(8))
                
                #Proportion text size and pad between axis and Y label
                tickLenght=15
                prop=15/8.
                if self.YAxisTicks.isChecked()==True:
                    bx.tick_params(labelcolor="none",length=prop*float(self.TSize.value()),color='white')
                    
                else:
                    bx.get_yaxis().set_ticks([])
                    bx.set_ylabel("Frequency", labelpad=7, fontsize=float(self.TSize.value()))
                    if self.Hist.isChecked()==False:
                        ax3.get_yaxis().set_ticks([])
                    else:
                        ax.get_yaxis().set_ticks([])
            else:
                if self.YAxisTicks.isChecked()==False:
                    if self.Hist.isChecked()==False:
                        ax.get_yaxis().set_ticks([])
                        ax.set_ylabel("Frequency", fontsize=float(self.TSize.value()))
                        ax3.get_yaxis().set_ticks([])
                    else:
                        ax1.get_yaxis().set_ticks([])
                        #ax1.set_ylabel("Frequency", fontsize=float(self.TSize.value()))
                        ax.get_yaxis().set_ticks([])
            
           
            #Increment Plot Counter
            plotCounter+=1
             

#New
import shutil,os
import sys

def copytree2(source,dest):
    #os.mkdir(dest)
    dest_dir = os.path.join(dest,os.path.basename(source))
    shutil.copytree(source,dest_dir)
    

if __name__ == '__main__':
    app= ApplicationContext()       # 1. Instantiate ApplicationContext
    form = GeochronologyPlots()
    
    #Get Resource Folders
    #Experimental recognision of user folder
    userDir=os.path.expanduser('~/')
    #Create Folders for Application deployment if not already created
    newpath = userDir+"CGD_DataStructure/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    #Initalize Directories for Structured Data
    Data_Pa=app.get_resource('./Data/')
    DataSet_Pa=app.get_resource('./Datasets/')
    Results_Pa=app.get_resource('./Results/')

    #Copy Sample data from Application resources folder to User local machine
    if not os.path.exists(newpath+'Data'):
        copytree2(Data_Pa,newpath)
    if not os.path.exists(newpath+'Datasets'):
        copytree2(DataSet_Pa,newpath)
    if not os.path.exists(newpath+'Results'):
        copytree2(Results_Pa,newpath)
    #Launching App
    form.show()
    exit_code = app.app.exec_()      # 2. Invoke appctxt.app.exec_()
    #When closing the app
    sys.exit(exit_code)
