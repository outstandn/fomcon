import sys
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtWidgets import QMessageBox
import numpy as np
from addict import Dict
from decimal import *

from fotf import *
from fopid_control import fopidcontrollerprototypemarkii as fofopdtTune, test_control
from pyGui.fomconoptimizegui import *
import asyncio

#gui
from pyGui import fopidopt, createnewfofopdtgui
#Constants
fofopdtModel = dict(K = 0, L = 0, T = 0, alpha = 0)
oustaloopModel = dict(wb = 0,wh = 0, N = 0,)
ALPHA_MIN = 0.001
ALPHA_MAX = 2
MAX_LAMBDA = 5
MIN_COEF = -1000
MAX_COEF = 1000
MIN_EXPO = 0.001
MAX_EXPO = 5
MAX_ITER = 50
MAX_SAMPLERATE = 1000
MIN_SAMPLERATE = 10
MAX_GAINMARGIN = 360
MAX_PHASEMARGIN = 1000
STATUSBAR_TIME = 7000
MIN_PORT = 49152
MAX_PORT = 65535

class fofopdtguiclass(QMainWindow, fopidopt.Ui_FOPIDOPT):
    def __init__(self):
        QMainWindow.__init__(self)
        fopidopt.Ui_FOPIDOPT.__init__(self)
        self.setupUi(self)
        self.setWindowIcon(QIcon('index.png'))
        self.IdentifiedModel = None

        self.comboBoxFOFOPDTSYS.currentIndexChanged.connect(self._comboBoxFOFOPDTSYSEmpty)
        # self.comboBoxApproxFilter.currentIndexChanged.connect(self._oustaloopSelected) #TODO :CHECK THIS
        self.lineEdit_StartFreq.textChanged.connect(self._OustaStartFreqChanged)
        self.lineEdit_StopFreq.textChanged.connect(self._OustaStopFreqChanged)
        self.lineEditOrder.textChanged.connect(self._OustaOrderChanged)
        self.pushButtonAddData.clicked.connect(self._addData)
        self.pushButtonEditData.clicked.connect(self._edit)
        self.pushButtonDeleteData.clicked.connect(self._deleteData)
        self.pushButtonTune.clicked.connect(self._tune)

        self.comboBoxTuneOption.currentIndexChanged.connect(self._TuneOptionChanged)
        # self.lineEdit_Kp.event.connect(self._KpChanged)
        self.lineEdit_Kp.textChanged.connect(self._KpChanged)
        self.lineEditConstMinKp.textChanged.connect(self._KpChanged)
        self.lineEditConstMaxKp.textChanged.connect(self._KpChanged)

        self.lineEdit_Ki.textChanged.connect(self._KiChanged)
        self.lineEditConstMinKi.textChanged.connect(self._KiChanged)
        self.lineEditConstMaxKi.textChanged.connect(self._KiChanged)

        self.lineEdit_Kd.textChanged.connect(self._KdChanged)
        self.lineEditConstMinKd.textChanged.connect(self._KdChanged)
        self.lineEditConstMaxKd.textChanged.connect(self._KdChanged)

        self.lineEdit_Lam.textChanged.connect(self._LamdaChanged)
        self.lineEditConstMinlam.textChanged.connect(self._LamdaChanged)
        self.lineEditConstMaxlam.textChanged.connect(self._LamdaChanged)

        self.lineEdit_Mu.textChanged.connect(self._MuChanged)
        self.lineEditConstMinMu.textChanged.connect(self._MuChanged)
        self.lineEditConstMaxMu.textChanged.connect(self._MuChanged)

        self.pushButtonSetParams.clicked.connect(self._setParams)
        self.lineEditParamsValue.textChanged.connect(self._paramValueChanged)

        self.lineEditIter.textChanged.connect(self._iterChanged)
        self.lineEditSamplRate.textChanged.connect(self._sampRateChanged)
        self.lineEditGainMargin.textChanged.connect(self._gainMargChanged)
        self.lineEditPhaseMargin.textChanged.connect(self._phaseMargChanged)


        self._isMuOk = self._isKdOk = self._isKiOk = self._isKpOk = True
        self._isOustaStartFreq = self._isOustaStopFreqOK = self._isOustaOrderOk = self._isLamdaOk = True
        self._isIterOk = self._isSampleRateOk = self._isGainMargOk = self._isPhaseMargOk = True
        self.comboFOFOPDTOk = False

        #TankControl Checkers and Buttons
        self.pushButtonTestConnection.clicked.connect(self._testControlConnection)
        self.pushButtonStartCom.clicked.connect(self.startNetContrl)
        self.pushButtonStopCom.clicked.connect(self.stopNetContrl)
        self.lineEditRecieveIP.textChanged.connect(self._ipRecChanged)
        self.lineEditRecievePort.textChanged.connect(self._recPortChanged)
        self.lineEditSendIP.textChanged.connect(self._ipSendChanged)
        self.lineEditSendPort.textChanged.connect(self._ipSendChanged)

        self._reloadAllFOTransFunc()
        self.show()

    #region Initialization Functions
    def _reloadAllFOTransFunc(self):
        # Startup Config
        fix = { 'Gains': [self.lineEdit_Kp, self.lineEdit_Ki, self.lineEdit_Kd],
                'Exponenets': [self.lineEdit_Lam, self.lineEdit_Mu],
                'Gains & Expo': [self.lineEdit_Kp, self.lineEdit_Ki, self.lineEdit_Kd, self.lineEdit_Lam, self.lineEdit_Mu]}

        algo = {"Tune Coefficient": optFix.Exp}#, "Tune All Parameters": optFix.Free,  "Tune Exponents": optFix.Coeff}
        # method = {"Grunwald Letnikov": optMethod.grunwaldLetnikov, "Oustaloop": optMethod.oustaloop}
        method = {"Oustaloop": optMethod.oustaloop}
        for i in fix:
            self.comboBoxParamSetOption.addItem(i, fix[i])

        for i in algo:
            self.comboBoxTuneOption.addItem(i, algo[i])

        for i in method:
            self.comboBoxApproxFilter.addItem(i, method[i])
        self.comboBoxParamSetOption.setCurrentIndex(0)
        self.groupBoxFOPIDParams.setChecked(False)
        self.groupBoxSimParams.setChecked(False)
        self.groupBoxGainPhaseMargin.setChecked(False)
        self.groupBoxTankControl.setChecked(False)
    #endregion

    #region Exception Handling, See fotf.py and pyfomcon.py for example
    def _ShowError(self, message, obj=None, obj2=None):
        try:
            if self.isDialogActive is False:
                self.isDialogActive = True
                if obj is not None:
                    obj.setCursorPosition(0)
                    obj.setSelection(0, len(obj.text()))

                self.statusbar.showMessage('Error: ' + message, STATUSBAR_TIME)
                raise ValueError(
                    QMessageBox.question(self, 'Error', message, QMessageBox.StandardButtons(QMessageBox.Ok)))

        except Exception as excep:
            self.isDialogActive = False
            print("\nError occured @ fofopdtguiclass._ShowError\n")
            print(type(excep))

    #region Button Checks
    # TODO :CHECK THIS _oustaloopSelected to be uysed lated for added aproxiamte methods features
    # def _oustaloopSelected(self):
    #     if self.comboBoxApproxFilter.currentData() is optMethod.grunwaldLetnikov:
    #         self.lineEdit_StartFreq.setEnabled(False)
    #         self.lineEdit_StopFreq.setEnabled(False)
    #         self.lineEditOrder.setEnabled(False)
    #         self._isOustaStartFreq = self._isOustaStopFreqOK = self._isOustaOrderOk = True
    #         try:
    #             self._ok2Tune()
    #         except:
    #             pass
    #     else:
    #         self.lineEdit_StartFreq.setEnabled(True)
    #         self.lineEdit_StopFreq.setEnabled(True)
    #         self.lineEditOrder.setEnabled(True)
    #         try:
    #             self._ok2Tune()
    #         except:
    #             pass

    def _OustaStartFreqChanged(self):
        try:
            x = int(self.lineEdit_StartFreq.text())
            if x < 0:
                self._isOustaStartFreq = True
            else:
                self._isOustaStartFreq = False
            self._ok2Tune()
        except:
            self._isOustaStartFreq = False

    def _OustaStopFreqChanged(self):
        try:
            x = int(self.lineEdit_StopFreq.text())
            if x > 0:
                self._isOustaStopFreqOK = True
            else:
                self._isOustaStopFreqOK = False
            self._ok2Tune()
        except:
            self._isOustaStopFreqOK = False

    def _OustaOrderChanged(self):
        try:
            x = int(self.lineEditOrder.text())
            if 0 < x < fofopdtTune.NMAX:
                self._isOustaOrderOk = True
            else:
                self._isOustaOrderOk = False
            self._ok2Tune()
        except:
            self._isOustaOrderOk = False

    def _paramValueChanged(self):
        try:
            x = self.lineEditParamsValue.text()
            if ALPHA_MIN <= float(x) <= ALPHA_MAX and self.groupBoxFOPIDParams.isChecked():
                self.pushButtonSetParams.setEnabled(True)
            else:
                self.pushButtonSetParams.setEnabled(False)
                self.statusbar.showMessage('{2} is out of range range: {0} <= x <= {1}'.format(ALPHA_MIN,ALPHA_MAX, x), STATUSBAR_TIME)
                print('{2} is out of range range: {0} <= x <= {1}'.format(ALPHA_MIN,ALPHA_MAX, x))
        except:
            self.pushButtonSetParams.setEnabled(False)
            self.statusbar.showMessage('fofopdtguiclass._qChanged:{2} is out of range range: {0} <= x <= {1}'.format(ALPHA_MIN, ALPHA_MAX, x), STATUSBAR_TIME)
            print('fofopdtguiclass._qChanged:{2} is out of range range: {0} <= x <= {1}'.format(ALPHA_MIN,ALPHA_MAX, x))

    def _TuneOptionChanged(self):
        currentTuneSettings = self.comboBoxTuneOption.currentData()
        if currentTuneSettings is optFix.Exp :
            self.lineEdit_Kp.setEnabled(True)
            self.lineEditConstMinKp.setEnabled(True)
            self.lineEditConstMaxKp.setEnabled(True)

            self.lineEdit_Ki.setEnabled(True)
            self.lineEditConstMinKi.setEnabled(True)
            self.lineEditConstMaxKi.setEnabled(True)

            self.lineEdit_Kd.setEnabled(True)
            self.lineEditConstMinKd.setEnabled(True)
            self.lineEditConstMaxKd.setEnabled(True)

            self.lineEdit_Lam.setEnabled(False)
            self.lineEditConstMinlam.setEnabled(False)
            self.lineEditConstMaxlam.setEnabled(False)

            self.lineEdit_Mu.setEnabled(False)
            self.lineEditConstMinMu.setEnabled(False)
            self.lineEditConstMaxMu.setEnabled(False)

        elif currentTuneSettings is optFix.Coeff:
            self.lineEdit_Kp.setEnabled(False)
            self.lineEditConstMinKp.setEnabled(False)
            self.lineEditConstMaxKp.setEnabled(False)

            self.lineEdit_Ki.setEnabled(False)
            self.lineEditConstMinKi.setEnabled(False)
            self.lineEditConstMaxKi.setEnabled(False)

            self.lineEdit_Kd.setEnabled(False)
            self.lineEditConstMinKd.setEnabled(False)
            self.lineEditConstMaxKd.setEnabled(False)

            self.lineEdit_Lam.setEnabled(True)
            self.lineEditConstMinlam.setEnabled(True)
            self.lineEditConstMaxlam.setEnabled(True)

            self.lineEdit_Mu.setEnabled(True)
            self.lineEditConstMinMu.setEnabled(True)
            self.lineEditConstMaxMu.setEnabled(True)

        elif currentTuneSettings is optFix.Free:
            self.lineEdit_Kp.setEnabled(False)
            self.lineEditConstMinKp.setEnabled(False)
            self.lineEditConstMaxKp.setEnabled(False)

            self.lineEdit_Ki.setEnabled(False)
            self.lineEditConstMinKi.setEnabled(False)
            self.lineEditConstMaxKi.setEnabled(False)

            self.lineEdit_Kd.setEnabled(False)
            self.lineEditConstMinKd.setEnabled(False)
            self.lineEditConstMaxKd.setEnabled(False)

            self.lineEdit_Lam.setEnabled(False)
            self.lineEditConstMinlam.setEnabled(False)
            self.lineEditConstMaxlam.setEnabled(False)

            self.lineEdit_Mu.setEnabled(False)
            self.lineEditConstMinMu.setEnabled(False)
            self.lineEditConstMaxMu.setEnabled(False)

    def _KpChanged(self):
        try:
            x = float(self.lineEdit_Kd.text())
            kpmin = float(self.lineEditConstMinKd.text())
            kpmax = float(self.lineEditConstMaxKd.text())
            if MIN_COEF <= kpmin <= x <= kpmax <= MAX_COEF:
                self._isKpOk = True
            else:
                self._isKpOk = False
        except:
            self._isKpOk = False
        finally:
            self._ok2Tune()

    def _KiChanged(self):
        try:
            x = float(self.lineEdit_Ki.text())
            kimin = float(self.lineEditConstMinKi.text())
            kimax = float(self.lineEditConstMaxKi.text())
            if MIN_COEF <= kimin <= x <= kimax <= MAX_COEF:
                self._isKiOk = True
            else:
                self._isKiOk = False
        except:
            self._isKiOk = False
        finally:
            self._ok2Tune()

    def _KdChanged(self):
        try:
            x = float(self.lineEdit_Kd.text())
            kdmin = float(self.lineEditConstMinKd.text())
            kdmax = float(self.lineEditConstMaxKd.text())
            if MIN_COEF <= kdmin <= x <= kdmax <= MAX_COEF:
                self._isKdOk = True
            else:
                self._isKdOk = False
        except:
            self._isKdOk = False
        finally:
            self._ok2Tune()

    def _LamdaChanged(self):
        try:
            x = float(self.lineEdit_Lam.text())
            lamdamin = float(self.lineEditConstMinlam.text())
            lamdamax = float(self.lineEditConstMaxlam.text())

            if MIN_EXPO < lamdamin <= x <= lamdamax <= MAX_EXPO:
                self._isLamdaOk = True
            else:
                self._isLamdaOk = False
        except:
            self._isLamdaOk = False
        finally:
            self._ok2Tune()

    def _MuChanged(self):
        try:
            x = float(self.lineEdit_Mu.text())
            mumin = float(self.lineEditConstMinMu.text())
            mumax = float(self.lineEditConstMaxMu.text())
            if MIN_EXPO < mumin <= x <= mumax <= MAX_EXPO:
                self._isMuOk = True
            else:
                self._isMuOk = False
        except:
            self._isMuOk = False
        finally:
            self._ok2Tune()

    def _ok2Tune(self):
        try:
            if self._isKpOk and self._isKiOk and self._isKdOk and  self._isLamdaOk and self._isMuOk and self.groupBoxSimParams.isChecked() \
                and self.comboFOFOPDTOk and self._isOustaStartFreq and self._isOustaStopFreqOK and self._isOustaOrderOk \
                    and self._isIterOk and self._isSampleRateOk and self.groupBoxSimParams.isChecked() \
                    and self._isGainMargOk and self._isPhaseMargOk and self.groupBoxGainPhaseMargin.isChecked():
                self.pushButtonTune.setEnabled(True)
            else:
                self.pushButtonTune.setEnabled(False)
        except:
            self.pushButtonTune.setEnabled(False)

    def _iterChanged(self):
        try:
            if 0 < int(self.lineEditIter.text()) <= MAX_ITER:
                self._isIterOk = True
            else:
                self._isIterOk = False
                print("0 < 'No. of Iteration' <= {0}".format(MAX_ITER))
                self.statusbar.showMessage("{0} < 'No. of Iteration' <= {1}".format(0,MAX_ITER), STATUSBAR_TIME)
        except:
            self._isIterOk = False
            print("'No. of Iteration' should be an integer")
            self.statusbar.showMessage("'No. of Iteration' should be an integer", STATUSBAR_TIME)
        finally:
            self._ok2Tune()

    def _sampRateChanged(self):
        try:
            if MIN_SAMPLERATE <= int(self.lineEditSamplRate.text()) <= MAX_SAMPLERATE:
                self._isSampleRateOk = True
            else:
                self._isSampleRateOk = False
                print("{0} < 'Sample Rate' <= {1}".format(MIN_SAMPLERATE, MAX_SAMPLERATE))
                self.statusbar.showMessage("{0} < 'Sample Rate' <= {1}".format(MIN_SAMPLERATE,MAX_SAMPLERATE), STATUSBAR_TIME)
        except:
            self._isSampleRateOk = False
            print("'Sample Rate' should be an integer")
            self.statusbar.showMessage("'Sample Rate' should be an integer", STATUSBAR_TIME)
        finally:
            self._ok2Tune()

    def _gainMargChanged(self):
        try:
            if 0 < float(self.lineEditGainMargin.text()) <= MAX_GAINMARGIN:
                self._isGainMargOk = True
            else:
                self._isGainMargOk = False
                print("{0} < 'Gain margin[dB]' <= {1}".format(0, MAX_GAINMARGIN))
                self.statusbar.showMessage("{0} < 'Gain margin[dB]' <= {1}".format(0,MAX_GAINMARGIN), STATUSBAR_TIME)
        except:
            self._isGainMargOk = False
            print("{0} < 'Gain margin[dB]' <= {1}".format(0, MAX_GAINMARGIN))
            self.statusbar.showMessage("{0} < 'Gain margin[dB]' <= {1}".format(0, MAX_GAINMARGIN), STATUSBAR_TIME)
        finally:
            self._ok2Tune()

    def _phaseMargChanged(self):
        try:
            if 0 < float(self.lineEditPhaseMargin.text()) <= MAX_PHASEMARGIN:
                self._isPhaseMargOk = True
            else:
                self._isPhaseMargOk = False
                print("{0} < Phase margin[deg] <= {1}".format(0, MAX_PHASEMARGIN))
                self.statusbar.showMessage("{0} < 'Phase margin[deg]' <= {1}".format(0,MAX_PHASEMARGIN), STATUSBAR_TIME)
        except:
            self._isPhaseMargOk = False
            print("'Sample Rate' should be an integer")
            self.statusbar.showMessage("'Sample Rate' should be an integer", STATUSBAR_TIME)
        finally:
            self._ok2Tune()

    def _ok2TestCon(self):
        if self._isIPRecOk and self._isIPSendOk and self._isPortRecvOk and self._isPortSendOk and self.groupBoxFOPIDParams.isChecked():
            self.pushButtonTestConnection.setEnabled(True)
        else:
            self.pushButtonTestConnection.setEnabled(False)

    def _ipRecChanged(self):
        try:
            ipaddress = self.lineEditRecieveIP.text().strip(" ").split(".")
            if ipaddress.count() == 4:
                for a in ipaddress:
                    if 0 <= int(a) <= 255:
                        pass
                    else:
                        self._isIPRecOk = False
                        print("{0} is Invalid in the IP address".format(a))
                        self.statusbar.showMessage("{0} is Invalid in the IP address".format(a), STATUSBAR_TIME)
                self._isIPRecOk = True
                print("'{0}' has a Valid IP Address format".format(ipaddress))
            else:
                self._isIPRecOk = False
                print("'{0}' needs to be formated with 3 '.'s i.e 'X.X.X.X'".format(ipaddress))
        except:
            self._isIPRecOk = False
            print("{0} is NOT an IP Address".format(ipaddress))
            self.statusbar.showMessage("{0} is NOT an IP Address".format(ipaddress), STATUSBAR_TIME)
        finally:
            self._ok2TestCon()

    def _ipSendChanged(self):
        try:
            ipaddress = self.lineEditSendIP.text().strip(" ").split(".")
            if ipaddress.count() == 4:
                for a in ipaddress:
                    if 0 <= int(a) <= 255:
                        pass
                    else:
                        self._isIPSendOk = False
                        print("{0} is Invalid in the IP address".format(a))
                        self.statusbar.showMessage("{0} is Invalid in the IP address".format(a), STATUSBAR_TIME)
                self._isIPSendOk = True
                print("'{0}' has a Valid IP Address format".format(ipaddress))
            else:
                self._isIPSendOk = False
                print("'{0}' needs to be formated with 3 '.'s i.e 'X.X.X.X'".format(ipaddress))
        except:
            self._isIPSendOk = False
            print("{0} is NOT an IP Address".format(ipaddress))
            self.statusbar.showMessage("{0} is NOT an IP Address".format(ipaddress), STATUSBAR_TIME)
        finally:
            self._ok2TestCon()

    def _recPortChanged(self):
        try:
            portnum = self.lineEditRecievePort.text()
            if MIN_PORT <= int(portnum) <= MAX_PORT:
                 self._isPortRecvOk = True
            else:
                self._isPortRecvOk = False
                print("Port:{0} is out of range. '{1} <= port <= {2}'".format(portnum,MIN_PORT, MAX_PORT))
                self.statusbar.showMessage(
                    "Port:{0} is out of range. '{1} <= port <= {2}'".format(portnum,MIN_PORT, MAX_PORT), STATUSBAR_TIME)
        except:
            self._isPortRecvOk = False
            print("ERROR: {0} is NOT an integer".format(portnum))
            self.statusbar.showMessage("ERROR: {0} is NOT an integer".format(portnum), STATUSBAR_TIME)
        finally:
            self._ok2TestCon()

    def _sendPortChanged(self):
        try:
            portnum = self.lineEditSendPort.text()
            if MIN_PORT <= int(portnum) <= MAX_PORT:
                self._isPortSendOk = True
            else:
                self._isPortSendOk = False
                print("Port:{0} is out of range. '{1} <= port <= {2}'".format(portnum, MIN_PORT, MAX_PORT))
                self.statusbar.showMessage(
                    "Port:{0} is out of range. '{1} <= port <= {2}'".format(portnum, MIN_PORT, MAX_PORT), STATUSBAR_TIME)
        except:
            self._isPortSendOk = False
            print("ERROR: {0} is NOT an integer".format(portnum))
            self.statusbar.showMessage("ERROR: {0} is NOT an integer".format(portnum), STATUSBAR_TIME)
        finally:
            self._ok2TestCon()
    #endregion

    #region Button Clicked Functions

    def _addData(self):
        _loadData = newfofopdtguiclass()
        _loadData.setFocus()
        _loadData.exec_()

        try:
            _sysname = _loadData.lineEditSysName.text()
            _pandasData = Dict(fofopdtModel)
            _pandasData.K = float(_loadData.lineEdit_GainK.text())
            _pandasData.L = float(_loadData.lineEdit_DelayText.text())
            _pandasData.T = float(_loadData.lineEdit_TimeConstant.text())
            _pandasData.alpha = float(_loadData.lineEdit_OrderAlpha.text())
            if _sysname:
                self.comboBoxFOFOPDTSYS.addItem(_sysname, _pandasData)
                self.comboBoxFOFOPDTSYS.setCurrentIndex(int(self.comboBoxFOFOPDTSYS.count()) - 1)
        except:
            self.statusbar.showMessage('fofopdtguiclass._addData: FOFOPDT Addition Failed', STATUSBAR_TIME)
            print('\nfofopdtguiclass._addData: FOFOPDT Addition Failed\n')

    def _deleteData(self):
        self.comboBoxFOFOPDTSYS.removeItem(self.comboBoxFOFOPDTSYS.currentIndex())

    def _comboBoxFOFOPDTSYSEmpty(self):
        if self.comboBoxFOFOPDTSYS.count() == 0:
            self.pushButtonDeleteData.setEnabled(False)
            self.pushButtonEditData.setEnabled(False)
            self.pushButtonTune.setEnabled(False)
            self.comboFOFOPDTOk = False
        else:
            self.comboFOFOPDTOk = True
            self.pushButtonDeleteData.setEnabled(True)
            self.pushButtonEditData.setEnabled(True)

        self._ok2Tune()

    def _tune(self):
        try:
            currentContext = getcontext()
            newContext = Context(prec=5, rounding=ROUND_HALF_EVEN, Emin=-999999, Emax=999999, capitals=1, clamp=0,
                                 flags=[], traps=[InvalidOperation, DivisionByZero, Overflow])
            setcontext(newContext)

            fofopdtTune.ACTIVATE_TUNING = True
            fofopdtTune.SAMPLE_RATE = int(self.lineEditSamplRate.text())   #very important to be first
            # fofopdtTune.NMAX = int(self.lineEditIter.text())                #Max allowed N for oustaloop approximation
            fofopdtTune.OPT_MAX_ITER = int(self.lineEditIter.text())        #Max iteration during tuning

            #get oustaloop model
            oustaloopModel = Dict(dict(wb=10**float(self.lineEdit_StartFreq.text()),  wh=10**float(self.lineEdit_StopFreq.text()),
                                        N=int(self.lineEditOrder.text()),  Ts= 1/float(self.lineEditSamplRate.text())))    #step = 1/samplerate
            #get fopidGuessModel
            fopidGuessModel = Dict(dict(Kp = float(self.lineEdit_Kp.text()), Ki = float(self.lineEdit_Ki.text()),
                                        Kd = float(self.lineEdit_Kd.text()),  lamda = float(self.lineEdit_Lam.text()),
                                        mu = float(self.lineEdit_Mu.text())))
            #get tuning/design parameter
            tunninParams = Dict(dict(wc = 0.1, pm = float(self.lineEditPhaseMargin.text()), opt_norm = fofopdtTune.OPT_NORM))

            #get Current FOFOPDT Model
            currentModel = self.comboBoxFOFOPDTSYS.currentData()

            #initialize the tuningProcess #TODO: if possible in another thread os as a task
            fofopdtTune.mainFOFOPIDOPT(currentModel, fopidGuessModel, oustaloopModel, tunninParams)

            # #change context for 4decimal place
            # precision4Context = Context(prec=4, rounding=ROUND_HALF_EVEN, Emin=-999999, Emax=999999, capitals=1, clamp=0,
            #                      flags=[], traps=[InvalidOperation, DivisionByZero, Overflow])
            # setcontext(precision4Context)

            self.lineEdit_Kp.setText(str(fofopdtTune.und_fopid.Kp))
            self.lineEdit_Ki.setText(str(fofopdtTune.und_fopid.Ki))
            self.lineEdit_Kd.setText(str(fofopdtTune.und_fopid.Kd))
            self.lineEdit_Lam.setText(str(fofopdtTune.und_fopid.lamda))
            self.lineEdit_Mu.setText(str(fofopdtTune.und_fopid.mu))

            setcontext(currentContext)
        except:
            print("An exception occurred. Try using another limit/ initial guess settings")
            self.statusbar.showMessage("An exception occurred. Try using another limit/ initial guess settings", 10000)

    async def startNetContrl(self):
        test_control.ALLOW_TO_SEND_RECV = True
        await test_control.startControl()

    def stopNetContrl(self):
        test_control.ALLOW_TO_SEND_RECV = False

    def _testControlConnection(self):
        recvIP, recvPort = self.lineEditRecieveIP.text(), int(self.lineEditRecievePort.text())
        sendIP, sendPort = self.lineEditSendIP.text(), int(self.lineEditSendPort.text())
        ok2StartContr = test_control.testUdpConnection(recvIP,recvPort, sendIP, sendPort)
        if ok2StartContr:
            self.pushButtonStartCom.setEnabled(True)
        else:
            self.pushButtonStartCom.setEnabled(False)

    def _edit(self):
        # read current name and Data
        currentText = self.comboBoxFOFOPDTSYS.currentText()
        currentData = self.comboBoxFOFOPDTSYS.currentData()

        # get datavalues for u,t,y
        K = currentData.K
        L = currentData.L
        T = currentData.T
        alpha = currentData.alpha

        # create trim data class
        _trista = newfofopdtguiclass()
        _trista.lineEditSysName.setText(currentText)
        _trista.lineEdit_GainK.setText(str(K))
        _trista.lineEdit_DelayText.setText(str(L))
        _trista.lineEdit_TimeConstant.setText(str(T))
        _trista.lineEdit_OrderAlpha.setText(str(alpha))

        # exec the trimdata classs
        _trista.setFocus()
        _trista.exec_()

        # now trimdata is exited get new values
        try:
            _sysname = _trista.lineEditSysName.text()
            _pandasData = Dict(fofopdtModel)
            _pandasData.K = float(_trista.lineEdit_GainK.text())
            _pandasData.L = float(_trista.lineEdit_DelayText.text())
            _pandasData.T = float(_trista.lineEdit_TimeConstant.text())
            _pandasData.alpha = float(_trista.lineEdit_OrderAlpha.text())
            if _sysname:
                self.comboBoxFOFOPDTSYS.addItem(_sysname, _pandasData)
                self.comboBoxFOFOPDTSYS.setCurrentIndex(int(self.comboBoxFOFOPDTSYS.count()) - 1)
        except:
            self.statusbar.showMessage('fofopdtguiclass._addData: FOFOPDT Addition Failed', STATUSBAR_TIME)
            print('\nfofopdtguiclass._addData: FOFOPDT Addition Failed\n')

    def _setParams(self):
        data = self.comboBoxParamSetOption.currentData()
        for i in data:
            i.setText(self.lineEditParamsValue.text())
        self._KpChanged()
        self._KiChanged()
        self._KdChanged()
        self._LamdaChanged()
        self._MuChanged()

    def closeEvent(self, event):
        reply = QMessageBox.question(self, "Exit?", "Are you sure about this Exit?", QMessageBox.Yes | QMessageBox.No,
                                     QMessageBox.Yes)
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()
    #endregion

class newfofopdtguiclass(QDialog, createnewfofopdtgui.Ui_dialogCreateNewFOTF):
    def __init__(self):
        QDialog.__init__(self)
        createnewfofopdtgui.Ui_dialogCreateNewFOTF.__init__(self)
        self.setWindowIcon(QIcon('index.png'))
        self.setupUi(self)

        self.lineEditSysName.textChanged.connect(self._checkSysName)
        self.lineEdit_GainK.textChanged.connect(self._checkGain)
        self.lineEdit_DelayText.textChanged.connect(self._checkDelay)
        self.lineEdit_TimeConstant.textChanged.connect(self._checkTimeConstant)
        self.lineEdit_OrderAlpha.textChanged.connect(self._checkAlpha)
        self.pushButtonOK.clicked.connect(self.close)
        self.pushButtonCancel.clicked.connect(self.close)
        self.sysnamecheck = self.gaincheck = self.delaycheck = self.timeconstantcheck = self.alphacheck =False
        self.show()

    #region Button OK Check
    def _checkOkButton(self):
        if self.sysnamecheck and self.gaincheck and self.delaycheck and self.timeconstantcheck and self.alphacheck:
            self.pushButtonOK.setEnabled(True)
        else:
            self.pushButtonOK.setEnabled(False)
    #endregion

    #region lineEdit Values Check
    def _checkSysName(self):
        try:
            self.sysnamecheck = len(self.lineEditSysName.text().strip(" ")) >= 1
            self._checkOkButton()
        except:
            self.sysnamecheck = False
            self._checkOkButton()

    def _checkGain(self):
        try:
            self.gaincheck = isinstance(float(self.lineEdit_GainK.text()), float)
            self._checkOkButton()
        except:
            self.gaincheck = False
            self._checkOkButton()

    def _checkDelay(self):
        try:
            self.delaycheck = float(self.lineEdit_DelayText.text()) >= 0
            self._checkOkButton()
        except:
            self.delaycheck = False
            self._checkOkButton()

    def _checkTimeConstant(self):
        try:
            self.timeconstantcheck = isinstance(float(self.lineEdit_TimeConstant.text()),float)
            self._checkOkButton()
        except:
            self.timeconstantcheck = False
            self._checkOkButton()

    def _checkAlpha(self):
        try:
            self.alphacheck = ALPHA_MIN < float(self.lineEdit_OrderAlpha.text()) < ALPHA_MAX
            self._checkOkButton()
        except:
            self.alphacheck = False
            self._checkOkButton()
    #endregion

    def closeEvent(self, event):
        sender = self.sender().text()

        close = QMessageBox.question(self, "{0}?".format(sender),
                                     "Are you sure you would like to '{0}' this form?".format(sender),
                                     QMessageBox.Yes | QMessageBox.No)
        if close == QMessageBox.Yes:
            if sender == "OK":
                pass
            else:
                self.lineEditSysName.clear()
                self.lineEdit_GainK.clear()
                self.lineEdit_DelayText.clear()
                self.lineEdit_TimeConstant.clear()
                self.lineEdit_OrderAlpha.clear()
            event.accept()
        else:
            event.ignore()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    fomcon = fofopdtguiclass()
    app.exec_()