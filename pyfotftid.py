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
import fotftidgui
import loaddatagui
import  fotfviewergui

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

class fotftidguiclass(QMainWindow, fotftidgui.Ui_MainWindow_fotftid):
    def __init__(self):
        QMainWindow.__init__(self)
        fotftidgui.Ui_MainWindow_fotftid.__init__(self)
        self.setupUi(self)
        self.setWindowIcon(QIcon('index.png'))
        self.reloadAllFOTransFunc()
        self.comboBoxData.currentIndexChanged.connect(self.comboBoxDataEmpty)
        self.pushButtonDeleteData.clicked(self._deleteData)
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
        self.pushButtonGeneratGuess.clicked.connect(self._GeneratePolynomials)
        self.pushButtonAddData.clicked.connect(self._addData)

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
            _pandasData = pd.read_excel(_datapath)
            self.comboBoxData.addItem(_sysname, _pandasData)
        except:
            self.statusbar.showMessage('Data Addition Fialed', 7000)
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

    def closeEvent(self, event):
        reply = QMessageBox.question(self, "Exit?", "Would you like to exit?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            sys.exit()



if __name__ == "__main__":
    app = QApplication(sys.argv)
    fomcon = fotftidguiclass()
    app.exec_()