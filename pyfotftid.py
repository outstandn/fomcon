import sys
import numpy as np
import pandas as pd
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.uic import *
from PyQt5.QtWidgets import QMessageBox

from fotf import *
from fomconoptimize import  *
from matplotlib import pyplot as plt
import fotftidgui
import loaddatagui
import trimdatagui

class loadDataClass(QDialog, loaddatagui.Ui_LoadDataForm):
    def __init__(self):
        QDialog.__init__(self)
        loaddatagui.Ui_LoadDataForm.__init__(self)
        self.setupUi(self)
        self.setWindowIcon(QIcon('index.png'))
        self.pushButtonBrowse.clicked.connect(self._browse)
        self.pushButtonOK.clicked.connect(self.close)
        self.pushButtonCancel.clicked.connect(self.close)
        self.lineEditDataPath.textChanged.connect(self._checkText)
        self.lineEditSysName.textChanged.connect(self._checkText)
        self.show()

    def _browse(self):
        fname = QFileDialog.getOpenFileName(self, 'Open Data', QDir.currentPath(), "Excel Files (*.xls *.xlsx)")
        self.lineEditDataPath.setText(fname[0])

    def _checkText(self):
        if (self.lineEditSysName.text() is not (None or "")) and self.lineEditDataPath.text() is not (None or ""):
            self.pushButtonOK.setEnabled(True)
        else:
            self.pushButtonOK.setEnabled(False)

    def closeEvent(self, event):
        close = QMessageBox.question(self, "QUIT", "Are you sure you are done with this form?", QMessageBox.Yes | QMessageBox.No)
        if close == QMessageBox.Yes:
            sender = self.sender().text()
            if sender == "OK":
                pass
            else:
                self.lineEditSysName.clear()
                self.lineEditDataPath.clear()
            event.accept()
        else:
            event.ignore()

class trimDataClass(QDialog,trimdatagui.Ui_TrimDataForm):
    def __init__(self):
        QDialog.__init__(self)
        trimdatagui.Ui_TrimDataForm.__init__(self)
        self.setupUi(self)
        self.setWindowIcon(QIcon('index.png'))
        # self.lineEditT2.editingFinished.connect(self.editedT2)
        self.lineEditT1.editingFinished.connect(self.editedT1)
        self.lineEditNewName.editingFinished.connect(self.editedName)
        self.pushButtonOK.clicked.connect(self.close)
        self.pushButtonCancel.clicked.connect(self.close)
        self.t1edited = False
        self.t2edited = False
        self.nameEdited = False

        self.show()

    def editedT1(self):
        if 0 < float(self.lineEditT1.text()):
            self.t1edited = True
        else:
            self.t1edited = False
        self.edited()

    # def editedT2(self):
    #     if float(self.lineEditT2.text()) > float(self.lineEditT1.text()):
    #         self.t2edited = True
    #     else:
    #         self.t2edited = False
    #     self.edited()
    def editedName(self):
        self.nameEdited = True
        self.edited()

    def edited(self):
        if self.nameEdited or self.t1edited:
            self.pushButtonOK.setEnabled(True)
        else:
            self.pushButtonOK.setEnabled(False)

    def closeEvent(self, event):
        close = QMessageBox.question(self, "QUIT", "Are you sure you are done with this form?", QMessageBox.Yes | QMessageBox.No)
        if close == QMessageBox.Yes:
            sender = self.sender().text()
            if sender == "OK":
                pass
            event.accept()
        else:
            event.ignore()

class iddata():
    def __init(self):
        self.u
        self.v
        self.y

class fotftidguiclass(QMainWindow, fotftidgui.Ui_MainWindow_fotftid):
    def __init__(self):
        QMainWindow.__init__(self)
        fotftidgui.Ui_MainWindow_fotftid.__init__(self)
        self.setupUi(self)
        self.setWindowIcon(QIcon('index.png'))
        self.reloadAllFOTransFunc()
        self.comboBoxData.currentIndexChanged.connect(self.comboBoxDataEmpty)
        self.pushButtonGeneratGuess.clicked.connect(self._GeneratePolynomials)
        self.pushButtonAddData.clicked.connect(self._addData)
        self.pushButtonDeleteData.clicked.connect(self._deleteData)
        self.pushButtonPlotData.clicked.connect(self._plot)
        self.pushButtonTrimData.clicked.connect(self._trim)

        self.show()


    def reloadAllFOTransFunc(self):
        #Startup Config
        for i in optMethod:
            self.comboBoxOptTypeMethod.addItem(str(i), i)

        for i in optAlgo:
            self.comboBoxAlgorithm.addItem(str(i), i)

        for i in optFix:
            self.comboBoxOptFix.addItem(str(i), i)

        self.comboBoxPolesOrZeros.addItems(['Pole Polynomial','Zero Polynomial',])
        self.lineEditCoefLimitLower.setEnabled(False)
        self.lineEditCoefLimitUpper.setEnabled(False)
        self.lineEditExpLimitLower.setEnabled(False)
        self.lineEditExpLimitUpper.setEnabled(False)
        self.plainTextEdit_Zeros.setPlainText('1')
        self._GeneratePolynomials()

    def _GeneratePolynomials(self):
        q = float(self.lineEditCommesurateOder.text())
        order = int(self.lineEditFOTFOrder.text())
        poleOrZero = self.comboBoxPolesOrZeros.currentText()

        firstgen = np.arange(order+1, dtype=float)*q
        nnum = firstgen[::-1]
        num = np.ones_like(nnum)

        # Get fractional polynomial
        poly = poly2str(num, nnum)

        if poleOrZero == "Zero Polynomial":
            self.plainTextEdit_Zeros.setPlainText(poly)
        elif poleOrZero == 'Pole Polynomial':
            self.plainTextEdit_Poles.setPlainText(poly)

    def _addData(self):
        _loadData = loadDataClass()
        _loadData.exec_()

        try:
            _sysname = _loadData.lineEditSysName.text()
            _datapath = _loadData.lineEditDataPath.text()
            _pandasD = pd.read_excel(_datapath)
            _pandasData = iddata()
            _pandasData.y = np.array(_pandasD.y.values)
            _pandasData.u = np.array(_pandasD.u.values)
            _pandasData.t = np.array(_pandasD.t.values)

            del _pandasD
            self.comboBoxData.addItem(_sysname, _pandasData)
        except:
            self.statusbar.showMessage('Data Addition Failed', 7000)
            print('\nData Addition Failed\n')

    def _deleteData(self):
        self.comboBoxData.removeItem(self.comboBoxData.currentIndex())

    def comboBoxDataEmpty(self):
        if self.comboBoxData.count() == 0:
            self.pushButtonDeleteData.setEnabled(False)
            self.pushButtonPlotData.setEnabled(False)
            self.pushButtonTrimData.setEnabled(False)
            self.pushButtonModelValidate.setEnabled(False)
            self.pushButtonIdentify.setEnabled(False)
        else:
            self.pushButtonDeleteData.setEnabled(True)
            self.pushButtonPlotData.setEnabled(True)
            self.pushButtonTrimData.setEnabled(True)
            self.pushButtonModelValidate.setEnabled(True)
            self.pushButtonIdentify.setEnabled(True)

    def _plot(self):
        try:
            currentText = self.comboBoxData.currentText()
            currentData = self.comboBoxData.currentData()
            y = currentData.y
            u = currentData.u
            t = currentData.t

            if y.size != u.size or u.size != t.size or t.size != y.size:
                raise IOError("fotftidguiclass._plot: size of data idd are not the same. Kindly Fix your data")

            plt.figure(dpi=128)
            plt.subplot(2, 1, 1)
            plt.plot(t, y, 'b-')
            plt.title(currentText)
            plt.ylabel('output')
            plt.grid(True, axis='both', which='both')

            plt.subplot(2, 1, 2)
            plt.plot(t, u, 'r-')
            plt.xlabel('time (sec)')
            plt.ylabel('input')
            plt.grid(True, axis='both', which='both')
            plt.show()

            if self.pushButtonTrimData.isEnabled() == False:
                self.pushButtonTrimData.setEnabled(True)
        except:
            print("fotftidguiclass._plot: size of data idd are not the same. Kindly Fix your data")
            self.statusbar.showMessage("fotftidguiclass._plot: size of data idd are not the same. Kindly Fix your data", 7000)

    def _trim(self):
        currentIndex = self.comboBoxData.currentIndex()
        currentText = self.comboBoxData.currentText()
        currentData = self.comboBoxData.currentData()
        y = currentData.y
        u = currentData.u
        t = currentData.t

        _trista = trimDataClass()
        _trista.lineEditNewName.setText(currentText)
        t0=str(t[0])
        # t1=str(t[t.size-1])
        _trista.lineEditT1.setReadOnly(False)
        _trista.lineEditT1.setText(t0)
        # _trista.lineEditT2.setReadOnly(False)
        # _trista.lineEditT2.setText(t1)
        _trista.setFocus()
        _trista.exec_()

        newt1 = float(_trista.lineEditT1.text())
        # newt2 = float(_trista.lineEditT2.text())
        newText = _trista.lineEditNewName.text()

        truncy = np.nonzero(t >= newt1)
        y = y[truncy]
        u = u[truncy]
        t = t[truncy]

        newdata = iddata()
        newdata.y = y
        newdata.u= u
        newdata.t = t

        self.comboBoxData.addItem(newText, newdata)
    def closeEvent(self, event):
        reply = QMessageBox.question(self, "Exit?", "Would you like to exit?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            sys.exit()



if __name__ == "__main__":
    app = QApplication(sys.argv)
    fomcon = fotftidguiclass()
    app.exec_()