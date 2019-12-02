import sys
import numpy as np
import pandas as pd
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.uic import *
from PyQt5.QtWidgets import QMessageBox

from fotf import *
from matplotlib import pyplot as plt
import fotftidgui
import loaddatagui
import trimdatagui
from fomconoptimizegui import *

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
        #self._checkText()

    def _checkText(self):
        if (self.lineEditSysName.text() is not (None or "")) and self.lineEditDataPath.text() is not (None or ""):
            self.pushButtonOK.setEnabled(True)
        else:
            self.pushButtonOK.setEnabled(False)

    def closeEvent(self, event):
        sender = self.sender().text()
        if sender == "Add":
            sender = "Close"

        close = QMessageBox.question(self, "{0}?".format(sender), "Are you sure you would like to '{0}' this form?".format(sender), QMessageBox.Yes | QMessageBox.No)
        if close == QMessageBox.Yes:
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
        self.lineEditT2.textChanged.connect(self.editedT2)
        self.lineEditT1.textChanged.connect(self.editedT1)
        self.lineEditNewName.textChanged.connect(self.editedName)
        self.pushButtonOK.clicked.connect(self.close)
        self.pushButtonCancel.clicked.connect(self.close)
        self.t1edited = False
        self.t2edited = False
        self.nameEdited = False
        self.initialt1 = -1
        self.initialt2 = -1
        self.initialName = ""

        self.show()

    def editedT1(self):
        try:
            if self.initialt1 <= float(self.lineEditT1.text()) < self.initialt2:
                if float(self.lineEditT1.text()) < float(self.lineEditT2.text()):
                    self.t1edited = True
                else:
                    self.t1edited = False
            else:
                self.t1edited = False
            self.edited()
        except:
            self.t1edited = False

    def editedT2(self):
        try:
            if self.initialt1 < float(self.lineEditT2.text()) <= self.initialt2:
                if float(self.lineEditT1.text()) < float(self.lineEditT2.text()):
                    self.t2edited = True
                else:
                    self.t2edited = False
            else:
                self.t2edited = False
            self.edited()
        except:
            self.t2edited = False



    def editedName(self):
        try:
            if self.initialName == self.lineEditNewName.text():
                self.nameEdited = False
            else:
                self.nameEdited = True
            self.edited()
        except:
            self.nameEdited = False


    def edited(self):
        if self.nameEdited and self.t1edited and self.t2edited:
            self.pushButtonOK.setEnabled(True)
        else:
            self.pushButtonOK.setEnabled(False)

    def closeEvent(self, event):
        sender = self.sender().text()
        if sender == "Trim":
            sender = "Close"
        close = QMessageBox.question(self, "{}?".format(sender), "Confirm {0}?".format(sender), QMessageBox.Yes | QMessageBox.No)
        if close == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

class fotftidguiclass(QMainWindow, fotftidgui.Ui_MainWindow_fotftid):
    def __init__(self):
        QMainWindow.__init__(self)
        fotftidgui.Ui_MainWindow_fotftid.__init__(self)
        self.setupUi(self)
        self.setWindowIcon(QIcon('index.png'))

        self.comboBoxData.currentIndexChanged.connect(self._comboBoxDataEmpty)
        self.pushButtonGeneratGuess.clicked.connect(self._GeneratePolynomials)
        self.pushButtonAddData.clicked.connect(self._addData)
        self.pushButtonDeleteData.clicked.connect(self._deleteData)
        self.pushButtonPlotData.clicked.connect(self._plot)
        self.pushButtonTrimData.clicked.connect(self._trim)
        self.checkBoxFixZeros.clicked.connect(self._fixedZero)
        self.checkBoxFixPoles.clicked.connect(self._fixedPole)
        self.checkBoxUseDelay.clicked.connect(self._useDelay)
        self.checkBoxUseCoefLimits.clicked.connect(self._useCoeffLim)
        self.checkBoxUseExpoLimits.clicked.connect(self._useExpoLim)
        self.pushButtonModelStability.clicked.connect(self._stabilityCheck)
        self.pushButtonRoundOff.clicked.connect(self._roundOff)
        self.pushButtonModelValidate.clicked.connect(self._validate)
        self.pushButtonIdentify.clicked.connect(self._identify)
        self.lineEditCommesurateOder.textChanged.connect(self._qChanged)
        self.lineEditFOTFOrder.textChanged.connect(self._fotfOrderChanged)
        self.lineEditCoefLimitLower.textChanged.connect(self._coeffLowerChanger)
        self.lineEditCoefLimitUpper.textChanged.connect(self._coeffUpperChanged)
        self.lineEditExpLimitLower.textChanged.connect(self._expLowerChanged)
        self.lineEditExpLimitUpper.textChanged.connect(self._expUpperChanged)
        self.lineEditOrder.textChanged.connect(self._OustaOrderChanged)
        self.lineEdit_StartFreq.textChanged.connect(self._OustaStartFreqChanged)
        self.lineEdit_StopFreq.textChanged.connect(self._OustaStopFreqChanged)
        self.lineEditLamda.textChanged.connect(self._LamdaChanged)
        self.textEdit_Zeros.textChanged.connect(self._LamdaChanged)
        self.textEdit_Poles.textChanged.connect(self._LamdaChanged)
        self.comboBoxOptTypeMethod.currentIndexChanged.connect(self._oustaloopSelected)
        self._isexpUpperok = self._isexpLowerok = self._iscoeffUpperok = self._iscoeffLowerok =True
        self._isOustaStartFreq = self._isOustaStopFreqOK = self._isOustaOrderOk = self._isLamdaOk = True
        self._isqok = self._isnOk = True
        self._reloadAllFOTransFunc()
        self.show()

    def _oustaloopSelected(self):
        if self.comboBoxOptTypeMethod.currentData() is optMethod.grunwaldLetnikov:
            self.lineEdit_StartFreq.setEnabled(False)
            self.lineEdit_StopFreq.setEnabled(False)
            self.lineEditOrder.setEnabled(False)
            self._isOustaStartFreq = self._isOustaStopFreqOK = self._isOustaOrderOk = True
            try:
                self._ok2Identify()
            except:
                pass
        else:
            self.lineEdit_StartFreq.setEnabled(True)
            self.lineEdit_StopFreq.setEnabled(True)
            self.lineEditOrder.setEnabled(True)
            try:
                self._ok2Identify()
            except:
                pass

    #LineEdit Checks

    def _LamdaChanged(self):
        try:
            x = int(self.lineEditLamda.text())
            z = self.textEdit_Zeros.toPlainText()
            p = self.textEdit_Poles.toPlainText()

            if 0<x<=5 and p is not "" and z is not "":
                self._isLamdaOk= True
                self.pushButtonRoundOff.setEnabled(True)
                self.pushButtonModelStability.setEnabled(True)
            else:
                self._isLamdaOk = False
                self.pushButtonRoundOff.setEnabled(False)
                self.pushButtonModelStability.setEnabled(False)
            self._ok2Identify()
        except:
            self._isLamdaOk = False
            self.pushButtonRoundOff.setEnabled(False)
            self.pushButtonModelStability.setEnabled(False)

    def _OustaStartFreqChanged(self):
        try:
            x = int(self.lineEdit_StartFreq.text())
            if x < 0 :
                self._isOustaStartFreq= True
            else:
                self._isOustaStartFreq = False
            self._ok2Identify()
        except:
            self._isOustaStartFreq = False

    def _OustaStopFreqChanged(self):
        try:
            x = int(self.lineEdit_StopFreq.text())
            if x > 0:
                self._isOustaStopFreqOK = True
            else:
                self._isOustaStopFreqOK = False
            self._ok2Identify()
        except:
            self._isOustaStopFreqOK = False

    def _OustaOrderChanged(self):
        try:
            x = int(self.lineEditOrder.text())
            if x > 0:
                self._isOustaOrderOk = True
            else:
                self._isOustaOrderOk = False
            self._ok2Identify()
        except:
            self._isOustaOrderOk = False

    def _qChanged(self):
        try:
            x = float(self.lineEditCommesurateOder.text())
            if 0.01 <= x <= 2:
                self._isqok = True
            else:
                self._isqok = False
            self._generateOk()
        except:
            self._isqok = False

    def _fotfOrderChanged(self):
        try:
            x = int(self.lineEditFOTFOrder.text())
            if 1 <= x <=10:
                self._isnOk = True
            else:
                self._isnOk = False
            self._generateOk()
        except:
            self._isnOk = False

    def _coeffLowerChanger(self):
        try:
            lower = int(self.lineEditCoefLimitLower.text())
            higher = int(self.lineEditCoefLimitUpper.text())
            if -20 <= lower <=20 and lower < higher:
                self._iscoeffLowerok = True
            else:
                self._iscoeffLowerok = False
            self._ok2Identify()
        except:
            self._iscoeffLowerok = False

    def _coeffUpperChanged(self):
        try:
            lower = int(self.lineEditCoefLimitLower.text())
            higher = int(self.lineEditCoefLimitUpper.text())
            if -20 <= higher <=20 and lower < higher:
                self._iscoeffUpperok = True
            else:
                self._iscoeffUpperok = False
            self._ok2Identify()
        except:
            self._iscoeffUpperok = False

    def _expLowerChanged(self):
        try:
            lower = int(self.lineEditExpLimitLower.text())
            higher = int(self.lineEditExpLimitUpper.text())
            if 0<= lower <= 10 and lower < higher:
                self._isexpLowerok = True
            else:
                self._isexpLowerok = False
            self._ok2Identify()
        except:
            self._isexpLowerok = False

    def _expUpperChanged(self):
        try:
            lower = int(self.lineEditExpLimitLower.text())
            higher = int(self.lineEditExpLimitUpper.text())
            if 0 <= higher <= 10 and lower < higher:
                self._isexpUpperok = True
            else:
                self._isexpUpperok = False
            self._ok2Identify()
        except:
            self._isexpUpperok = False

    def _generateOk(self):
        try:
            if self._isqok and self._isnOk:
                self.pushButtonGeneratGuess.setEnabled(True)
            else:
                self.pushButtonGeneratGuess.setEnabled(False)
        except:
            self.pushButtonGeneratGuess.setEnabled(False)

    def _ok2Identify(self):
        x = self.comboBoxData.count() > 0
        try:
            if x and self._iscoeffLowerok and self._iscoeffUpperok and self._isexpLowerok and self._isexpUpperok \
                    and self._isOustaStartFreq and self._isOustaStopFreqOK and self._isOustaOrderOk and self._isLamdaOk:
                self.pushButtonIdentify.setEnabled(True)
            else:
                self.pushButtonIdentify.setEnabled(False)
        except:
            self.pushButtonIdentify.setEnabled(False)

    #Checkboxes event handlers
    def _fixedZero(self):
        if self.checkBoxFixZeros.isChecked():
            self.textEdit_Zeros.setEnabled(False)
        else:
            self.textEdit_Zeros.setEnabled(True)

    def _fixedPole(self):
        if self.checkBoxFixPoles.isChecked():
            self.textEdit_Poles.setEnabled(False)
        else:
            self.textEdit_Poles.setEnabled(True)

    def _useDelay(self):
        if self.checkBoxUseDelay.isChecked():
            self.lineEdit_Delay.setEnabled(True)
        else:
            self.lineEdit_Delay.setEnabled(False)

    def _useCoeffLim(self):
        if self.checkBoxUseCoefLimits.isChecked():
            self.lineEditCoefLimitLower.setEnabled(True)
            self.lineEditCoefLimitUpper.setEnabled(True)
        else:
            self.lineEditCoefLimitLower.setEnabled(False)
            self.lineEditCoefLimitUpper.setEnabled(False)

    def _useExpoLim(self):
        if self.checkBoxUseExpoLimits.isChecked():
            self.lineEditExpLimitLower.setEnabled(True)
            self.lineEditExpLimitUpper.setEnabled(True)
        else:
            self.lineEditExpLimitLower.setEnabled(False)
            self.lineEditExpLimitUpper.setEnabled(False)

    def _reloadAllFOTransFunc(self):
        #Startup Config
        fix = {"Free Identification": optFix.Free, "Fix Coefficient":optFix.Coeff, "Fix Exponents": optFix.Exp}
        algo = {"Levenberg Marquardt": optAlgo.LevenbergMarquardt,
                "Trust Region Reflective": optAlgo.TrustRegionReflective,
                "Cauchy Robust Loss": optAlgo.RobustLoss,
                "softl1 Robust Loss": optAlgo.Softl1}
        method = {"Grunwald Letnikov": optMethod.grunwaldLetnikov, "Oustaloop": optMethod.oustaloop}
        #method = {"Grunwald Letnikov": optMethod.grunwaldLetnikov}
        for i in fix:
            self.comboBoxOptFix.addItem(i,fix[i])

        for i in algo:
            self.comboBoxAlgorithm.addItem(i,algo[i])

        for i in method:
            self.comboBoxOptTypeMethod.addItem(i,method[i])

        self.comboBoxPolesOrZeros.addItems(['Zero Polynomial','Pole Polynomial'])
        self.comboBoxPolesOrZeros.setCurrentIndex(1)
        # self.lineEditCoefLimitLower.setEnabled(False)
        # self.lineEditCoefLimitUpper.setEnabled(False)
        # self.lineEditExpLimitLower.setEnabled(False)
        # self.lineEditExpLimitUpper.setEnabled(False)
        # self.textEdit_Zeros.setPlainText('1')
        self._GeneratePolynomials()

    def _GeneratePolynomials(self):
        accu = int(self.lineEditLamda.text())
        q = float(self.lineEditCommesurateOder.text())
        order = int(self.lineEditFOTFOrder.text())
        poleOrZero = self.comboBoxPolesOrZeros.currentText()

        firstgen = np.arange(order+1, dtype=float)*q
        nnum = firstgen[::-1]
        num = np.ones_like(nnum)

        # Get fractional polynomial
        poly = poly2str(num, nnum, eps=accu)

        if poleOrZero == "Zero Polynomial":
            self.textEdit_Zeros.setPlainText(poly)
        elif poleOrZero == 'Pole Polynomial':
            self.textEdit_Poles.setPlainText(poly)

    def _addData(self):
        _loadData = loadDataClass()
        _loadData.exec_()

        try:
            _sysname = _loadData.lineEditSysName.text()
            _datapath = _loadData.lineEditDataPath.text()
            _pandasD = pd.read_excel(_datapath)
            _pandasData = idData()
            _pandasData.y = np.array(_pandasD.y.values)
            _pandasData.u = np.array(_pandasD.u.values)
            _pandasData.t = np.array(_pandasD.t.values)

            del _pandasD
            self.comboBoxData.addItem(_sysname, _pandasData)
            self.comboBoxData.setCurrentIndex(int(self.comboBoxData.count())-1)
        except:
            self.statusbar.showMessage('Data Addition Failed', 7000)
            print('\nData Addition Failed\n')

    def _deleteData(self):
        self.comboBoxData.removeItem(self.comboBoxData.currentIndex())

    def _comboBoxDataEmpty(self):
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
        #read current name and Data
        currentText = self.comboBoxData.currentText()
        currentData = self.comboBoxData.currentData()

        #get datavalues for u,t,y
        y = currentData.y
        u = currentData.u
        t = currentData.t

        #create trim data class
        _trista = trimDataClass()

        t0 = t[0]
        t1 = t[-1]
        # save initial value in  trimdata class as a variable
        _trista.initialName = currentText
        _trista.initialt1 = t0
        _trista.initialt2 = t1

        #fill in data value in trim data class
        _trista.lineEditNewName.setText(currentText)
        _trista.lineEditT1.setText(str(t0))
        _trista.lineEditT2.setText(str(t1))

        #exec the trimdata classs
        _trista.setFocus()
        _trista.exec_()

        #now trimdata is exited get new values
        newt1 = float(_trista.lineEditT1.text())
        newt2 = float(_trista.lineEditT2.text())
        newText = _trista.lineEditNewName.text()

        tlower = t >= newt1
        thigher = t <= newt2
        truncy = np.nonzero(tlower == thigher)

        newy = y[truncy]
        newu = u[truncy]
        newt = t[truncy]

        newdata = idData()
        newdata.y = newy
        newdata.u= newu
        newdata.t = newt

        if currentText == newText or newt.size == t.size:
            print("Data not Trimmed becasue there were no changes in Name or Time")
            self.statusbar.showMessage("Data not Trimmed becasue there were no changes in '{0}'".format(currentText), 10000)

        else:
            self.comboBoxData.addItem(newText, newdata)
            print("Trimmed '{}' to '{}' successfully".format(currentText,newText))
            self.statusbar.showMessage("Trimmed '{}' to '{}' successfully".format(currentText,newText), 10000)
            self.comboBoxData.setCurrentIndex(self.comboBoxData.count() - 1)

    def _stabilityCheck(self):
        try:

            _zero = self.textEdit_Zeros.toPlainText()
            _poles = self.textEdit_Poles.toPlainText()

            if self.checkBoxUseDelay.isChecked():
                _dt = float(self.lineEdit_Delay.text())
            else:
                _dt = 0

            identifiedSytem = newfotf(_zero,_poles,_dt)
            identifiedSytem.numberOfDecimal = int(self.lineEditLamda.text())
            isstabledata = identifiedSytem.isstable(True)
            num,nnum,den,nden,dt = fotfparam(identifiedSytem)
            newZero = poly2str(num,nnum,eps=identifiedSytem.numberOfDecimal)
            self.textEdit_Zeros.setText(newZero)

            newPole = poly2str(den,nden,eps=identifiedSytem.numberOfDecimal)
            self.textEdit_Poles.setText(newPole)

            if self.checkBoxUseDelay.isChecked():
                self.lineEdit_Delay.setText(_dt)
        except:
            pass


    def _roundOff(self):
        _zero = self.textEdit_Zeros.toPlainText()
        _poles = self.textEdit_Poles.toPlainText()
        epsi = int(self.lineEditLamda.text())

        semizero = str2poly(_zero)
        semipole = str2poly(_poles)

        newZero = poly2str(semizero[0], semizero[1], eps=epsi)
        self.textEdit_Zeros.setText(newZero)

        newPole = poly2str(semipole[0], semipole[1], eps=epsi)
        self.textEdit_Poles.setText(newPole)

    def _validate(self):
        #get verification Data
        currentIndex = self.comboBoxData.currentIndex()
        currentText = self.comboBoxData.currentText()
        verificationData = self.comboBoxData.currentData()

        #simulate verification data
        vy = verificationData.y
        vu = verificationData.u
        vt = verificationData.t


        #get identified system and step
        _zero = self.textEdit_Zeros.toPlainText()
        _poles = self.textEdit_Poles.toPlainText()

        if self.checkBoxUseDelay.isChecked():
            _dt = float(self.lineEdit_Delay.text())
        else:
            _dt = 0

        IdentifiedG = newfotf(_zero, _poles, _dt)
        IdentifiedG.numberOfDecimal = int(self.lineEditLamda.text())

        #simulate with the data from current combobox
        lsimG = lsim(IdentifiedG, vu, vt)

        # plot identified system output vs Data from initial system
        plt.figure(dpi=128)
        plt.subplot(2, 1, 1)
        plt.plot(vt, vy, 'b-', vt, lsimG, 'g-')
        plt.title("Validation Data '{0}' vs Identified System".format(currentText))
        plt.ylabel('output')
        plt.legend(['Valdata', 'idsystem'], loc='upper right')
        plt.grid(True, axis='both', which='both')

        # Fitness measure
        err = vy - lsimG
        fitness = 100 * (1 - (np.linalg.norm(err) / np.linalg.norm(vy - np.mean(lsimG))))

        plt.subplot(2, 1, 2)
        plt.plot(vt, vy - lsimG, 'r-')
        plt.title("Identified System error. fitness: {}%".format(round(fitness, 2)))
        plt.xlabel('time (sec)')
        plt.ylabel('error')
        plt.grid(True, axis='both', which='both')
        plt.show()

    def _identify(self):
        try:
            accu = int(self.lineEditLamda.text())
            #inital Guess
            initialGuess = newfotf(self.textEdit_Zeros.toPlainText(),self.textEdit_Poles.toPlainText())
            initialGuess.numberOfDecimal = accu
            #Similation method from combobox
            optimethod = self.comboBoxOptTypeMethod.currentData()

            if optimethod is optMethod.oustaloop:
                freqlower = int(self.lineEdit_StartFreq.text())
                freqhigher = int(self.lineEdit_StopFreq.text())
                freqOrder = int(self.lineEditOrder.text())
                oustoptions = np.array([freqlower, freqhigher, freqOrder],dtype=int)
            else:
                oustoptions = None

            #Identification algorithm from combobox
            optialg = self.comboBoxAlgorithm.currentData()

            #fix Coef or Expo or Free
            optiFix = self.comboBoxOptFix.currentData()

            #fix Zeros or Poles or None
            polyfix = [self.checkBoxFixZeros.isChecked(), self.checkBoxFixPoles.isChecked()]

            #optimize with delay?
            optiDelay = [self.checkBoxUseDelay.isChecked(), float(self.lineEdit_Delay.text())]

            findDelay = self.checkBoxUseDelay.isChecked()

            #generate optimization parameter class
            optiset = opt(initialGuess,optimethod,optialg,optiFix,polyfix, findDelay,oustaOption=oustoptions,accuracy= 10.0**-accu)

            #Get Currrent Data
            data = self.comboBoxData.currentData()

            #Get limits settings
            if self.checkBoxUseExpoLimits.isChecked():
                expLim = [int(self.lineEditExpLimitLower.text()),int(self.lineEditExpLimitUpper.text())]
            else:
                expLim = [10.0**-accu,10]

            if self.checkBoxUseCoefLimits.isChecked():
                coefLim = [int(self.lineEditCoefLimitLower.text()),int(self.lineEditCoefLimitUpper.text())]
            else:
                coefLim = [-10,10]

            # run Identification
            res = fid(data,optiset,[coefLim,expLim],plot=[False,False])
            res.numberOfDecimal = accu

            zero= res.poly2str(res.num[0][0], res.nnum[0][0])
            pole = res.poly2str(res.den[0][0], res.nden[0][0])
            dt = str(res.dt)


            #Update Text Boxes with identififed results
            self.textEdit_Zeros.setPlainText(zero)
            self.textEdit_Poles.setPlainText(pole)
            if self.lineEdit_Delay.isEnabled():
                self.lineEdit_Delay.setPlainText(dt)
        except:
            pass
            print("An exception occurred. Try using another limit/ initial guess settings")
            self.statusbar.showMessage("An exception occurred. Try using another settings",10000)


    def closeEvent(self, event):
        reply = QMessageBox.question(self, "Exit?", "Are you sure about Exit?", QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    fomcon = fotftidguiclass()
    app.exec_()



