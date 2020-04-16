import sys
import traceback
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import  *
from  multiprocessing import Process, Value, Array ,Pipe, Queue, Pool,Lock
from struct import pack, unpack

from addict import Dict
import socket
import time, datetime, traceback
from copy import deepcopy as copydeep

from fopid_control import fofopdtpidTuner as fofopdtTune
from pyGui.fomconoptimizegui import *
from fopid_control import controlServer

#gui
from pyGui import fopidoptgui, createnewfofopdtgui

#Constants
fofopdtModel = dict(K = 66.16, L = 1.93, T = 12.72, alpha = 0.5)
oustaloopModel = dict(wb = 0.0001,wh = 10000, N = 5,)
ALPHA_MIN,ALPHA_MAX,MAX_LAMBDA,MIN_COEF,MAX_COEF = 0.001,2,5,float('-inf'),float('inf')
MIN_EXPO,MAX_EXPO,MAX_ITER,MAX_DT,MIN_DT = 0.001,5,50,0.99, 0.01
MAX_GAINMARGIN,MAX_PHASEMARGIN,MAX_SIM_TIME,STATUSBAR_TIME = 20,359,3600,7000
MIN_PORT,MAX_PORT = 1081, 65535


class fofopdtguiclass(QMainWindow, fopidoptgui.Ui_FOPIDOPT):
    def __init__(self,guilock=None):
        QMainWindow.__init__(self)
        fopidoptgui.Ui_FOPIDOPT.__init__(self)
        self.setupUi(self)
        self.setWindowIcon(QIcon('index.png'))

        self.comboBoxFOFOPDTSYS.currentIndexChanged.connect(self._comboBoxFOFOPDTSYSEmpty)
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
        self.lineEditDt.textChanged.connect(self._dtChanged)
        self.lineEditTankSimTime.textChanged.connect(self._simTimeChanged)
        self.lineEditGainMargin.textChanged.connect(self._gainMargChanged)
        self.lineEditPhaseMargin.textChanged.connect(self._phaseMargChanged)

        self._isLamdaOk=self._isMuOk = self._isKdOk = self._isKiOk = self._isKpOk = True
        self._isOustaStartFreq = self._isOustaStopFreqOK = self._isOustaOrderOk =True
        self._isIterOk = self._isDtOk = self._isGainMargOk = self._isPhaseMargOk = self._isSimeTimeOk = True
        self.comboFOFOPDTOk = False

        #TankControl Checkers and Buttons
        self.pushButtonTestControl.clicked.connect(self._TestControl)
        self.pushButtonStartCom.clicked.connect(self.startNetContrl)
        self.pushButtonStopCom.clicked.connect(self.stopNetContrl)
        self.pushButtonUpdateTime.clicked.connect(self.updateTime)
        self.lineEditRecieveIP.textChanged.connect(self._ipRecChanged)
        self.lineEditRecievePort.textChanged.connect(self._recPortChanged)
        self.lineEditSendIP.textChanged.connect(self._ipSendChanged)
        self.lineEditSendPort.textChanged.connect(self._ipSendChanged)
        self.pushButtonToServer.clicked.connect(self._updateFOPIDServer)
        self.lineEditSendIPFOPIDContr.textChanged.connect(self._ipControlChanged)
        self.lineEditSendPortFOPIDContr.textChanged.connect(self._recvPortControlChanged)
        self.pushButtonExitServer.clicked.connect(self._exitServer)

        self._isPortRecvOk = self._isPortSendOk = self._isIPRecOk = \
        self._isIPSendOk = self._isIPControlOk = self._isPortControlOk =True

        self._reloadAllFOTransFunc()
        self.lock = guilock
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
                    QMessageBox.question(QWidget, 'Error', message, QMessageBox.StandardButtons(QMessageBox.Ok)))

        except Exception as excep:
            self.isDialogActive = False
            print("\nError occured @ fofopdtguiclass._ShowError\n")
            print(type(excep))

    #region Button Checks

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

            self.lineEdit_Lam.setEnabled(True)
            self.lineEditConstMinlam.setEnabled(False)
            self.lineEditConstMaxlam.setEnabled(False)

            self.lineEdit_Mu.setEnabled(True)
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

            self.lineEdit_Lam.setEnabled(True)
            self.lineEditConstMinlam.setEnabled(False)
            self.lineEditConstMaxlam.setEnabled(False)

            self.lineEdit_Mu.setEnabled(True)
            self.lineEditConstMinMu.setEnabled(False)
            self.lineEditConstMaxMu.setEnabled(False)

    def _KpChanged(self):
        try:
            x = float(self.lineEdit_Kp.text())
            kpmin = float(self.lineEditConstMinKp.text())
            kpmax = float(self.lineEditConstMaxKp.text())
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
            lammin = float(self.lineEditConstMinlam.text())
            lammax = float(self.lineEditConstMaxlam.text())

            if MIN_EXPO <= lammin <= x <= lammax <= MAX_EXPO:
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
            if MIN_EXPO <= mumin <= x <= mumax <= MAX_EXPO:
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
                and self._isIterOk and self._isDtOk and self.groupBoxSimParams.isChecked() \
                and self._isGainMargOk and self._isPhaseMargOk and self.groupBoxGainPhaseMargin.isChecked() \
                and self._isIPControlOk and self._isPortControlOk and self._isSimeTimeOk \
                and self.pushButtonTestControl.isEnabled()==False and self.groupBoxTankControl.isChecked():
                self.pushButtonTune.setEnabled(True)
                self.pushButtonToServer.setEnabled(True)
                self.pushButtonUpdateTime.setEnabled(True)
            else:
                self.pushButtonTune.setEnabled(False)
                self.pushButtonToServer.setEnabled(False)
                self.pushButtonUpdateTime.setEnabled(False)
        except:
            self.pushButtonTune.setEnabled(False)
            self.pushButtonToServer.setEnabled(False)
            self.pushButtonUpdateTime.setEnabled(False)

    def _iterChanged(self):
        iter = self.lineEditIter.text()
        try:
            if 1 <= int(iter) <= MAX_ITER:
                self._isIterOk = True
            else:
                self._isIterOk = False
                print("No. of Iteration: {0} is out of range. [{1} <= 'No. of Iteration' <= {2}]".format(iter,1,MAX_ITER))
                self.statusbar.showMessage("No. of Iteration: {0} is out of range. [{1} <= 'No. of Iteration' <= {2}]".format(iter,1,MAX_ITER), STATUSBAR_TIME)
        except:
            self._isIterOk = False
            print("ERROR!, 'No. of Iteration' Must be an integer")
            self.statusbar.showMessage("ERROR!, 'No. of Iteration' Must be an integer", STATUSBAR_TIME)
        finally:
            self._ok2Tune()

    def _dtChanged(self):
        step = self.lineEditDt.text()
        try:
            if MIN_DT <= float(step) <= MAX_DT:
                self._isDtOk = True
            else:
                self._isDtOk = False
                print("Step :{0}s is out of range. [{1} <= 'Step' <= {2}]".format(step,MIN_DT, MAX_DT))
                self.statusbar.showMessage("Step :{0}s is out of range. [{1} <= 'Step' <= {2}]".format(step,MIN_DT, MAX_DT), STATUSBAR_TIME)
        except:
            self._isDtOk = False
            print("ERROR!. Step :{0}s is NOT VALID".format(step))
            self.statusbar.showMessage("ERROR!. Step :{0}s is NOT VALID".format(step), STATUSBAR_TIME)
        finally:
            self._ok2Tune()

    def _gainMargChanged(self):
        gainMargin = self.lineEditGainMargin.text()
        try:
            if 0 < float(gainMargin) <= MAX_GAINMARGIN:
                self._isGainMargOk = True
            else:
                self._isGainMargOk = False
                print("Gain Margin: {0} dB, is out of range. [{1} < 'Gain margin[dB]' <= {2}]".format(gainMargin,0,MAX_GAINMARGIN))
                self.statusbar.showMessage("Gain Margin: {0} dB, is out of range. [{1} < 'Gain margin[dB]' <= {2}]".format(gainMargin,0,MAX_GAINMARGIN), STATUSBAR_TIME)
        except:
            self._isGainMargOk = False
            print("ERROR! Gain Margin: '{0}' dB is NOT VALID".format(gainMargin))
            self.statusbar.showMessage("ERROR! Gain Margin: '{0}' dB is NOT VALID".format(gainMargin), STATUSBAR_TIME)
        finally:
            self._ok2Tune()

    def _simTimeChanged(self):
        simtime = self.lineEditTankSimTime.text()
        try:
            if 0 < float(simtime) <= MAX_SIM_TIME and not self.pushButtonTestControl.isEnabled():
                self._isSimeTimeOk = True
            else:
                self._isSimeTimeOk = False

                print("Sim Time:{0}s is out of range. [{0} < 'Sim Time[s]' <= {1}]".format(simtime,0,MAX_SIM_TIME))
                self.statusbar.showMessage("Sim Time:{0}s is out of range. [{0} < 'Sim Time[s]' <= {1}]".format(simtime,0,MAX_SIM_TIME), STATUSBAR_TIME)
        except:
            self._isSimeTimeOk = False
            print("ERROR! Sim Time: '{0}'s is NOT VALID".format(simtime))
            self.statusbar.showMessage("ERROR! Sim Time: '{0}'s is NOT VALID".format(simtime), STATUSBAR_TIME)
        finally:
            self._ok2Tune()

    def _phaseMargChanged(self):
        phasemargin = self.lineEditPhaseMargin.text()
        try:
            if 1 <= float(phasemargin) <= MAX_PHASEMARGIN:
                self._isPhaseMargOk = True
            else:
                self._isPhaseMargOk = False
                print("Phase Margin:{0}[deg] is out opf range [{1} <= Phase margin[deg] <= {2}]".format(phasemargin,1, MAX_PHASEMARGIN))
                self.statusbar.showMessage("Phase Margin:{0}[deg] is out of range [{1} < Phase margin[deg] <= {2}]".format(phasemargin,1, MAX_PHASEMARGIN), STATUSBAR_TIME)
        except:
            self._isPhaseMargOk = False
            print("ERROR! 'Phase Margin':{0}[deg] is NOT VALID".format(phasemargin))
            self.statusbar.showMessage("ERROR! 'Phase Margin':{0}[deg] is NOT VALID".format(phasemargin), STATUSBAR_TIME)
        finally:
            self._ok2Tune()

    def _ok2TestCon(self):
        if self._isIPRecOk and self._isIPSendOk and self._isPortRecvOk and self._isPortSendOk and self.comboFOFOPDTOk \
                and self.groupBoxTankControl.isChecked() and self.groupBoxFOPIDParams.isChecked() \
                and self.lineEditRecieveIP.isEnabled() and self.lineEditSendIP.isEnabled()\
                and self.lineEditRecievePort.isEnabled()  and self.lineEditSendPort.isEnabled() and not self.pushButtonTestControl.isEnabled() :
            self.pushButtonStartCom.setEnabled(True)
        else:
            self.pushButtonStartCom.setEnabled(False)

    #region IP LINE EDIT CHECKERS
    def _ipRecChanged(self):
        ipaddress = self.lineEditRecieveIP.text().strip(" ").split(".")
        try:
            if len(ipaddress) == 4:
                for a in ipaddress:
                    if 0 <= int(a) < 255:
                        pass
                    else:
                        self._isIPRecOk = False
                        print("{0} is Invalid in the IP address".format(a))
                        self.statusbar.showMessage("{0} is Invalid in the IP address".format(a), STATUSBAR_TIME)
                self._isIPRecOk = True
                print("'{0}' VALID IP Address format".format(ipaddress))
                self.statusbar.showMessage("'{0}' VALID IP Address format".format(ipaddress), STATUSBAR_TIME)
            else:
                self._isIPRecOk = False
                print("'{0}' needs to be formated with 3 '.'s i.e 'X.X.X.X'".format(ipaddress))
        except:
            traceback.format_stack()
            self._isIPRecOk = False
            print("{0} is NOT an IP Address".format(ipaddress))
            self.statusbar.showMessage("{0} is NOT an IP Address".format(ipaddress), STATUSBAR_TIME)
        finally:
            self._ok2TestCon()

    def _ipSendChanged(self):
        ipaddress = self.lineEditSendIP.text().strip(" ").split(".")
        try:
            if len(ipaddress) == 4:
                for a in ipaddress:
                    if 0 <= int(a) < 255:
                        pass
                    else:
                        self._isIPSendOk = False
                        print("{0} is Invalid in the IP address".format(a))
                        self.statusbar.showMessage("{0} is Invalid in the IP address".format(a), STATUSBAR_TIME)
                self._isIPSendOk = True
                print("'{0}'  VALID IP Address format".format(ipaddress))
                self.statusbar.showMessage("'{0}' VALID IP Address format".format(ipaddress), STATUSBAR_TIME)
            else:
                self._isIPSendOk = False
                print("'{0}' needs to be formated with 3 '.'s i.e 'X.X.X.X'".format(ipaddress))
        except Exception as exce:
            traceback.format_stack()
            self._isIPSendOk = False
            print("{0} NOT an IP Address".format(ipaddress))
            self.statusbar.showMessage("{0} is NOT an IP Address".format(ipaddress), STATUSBAR_TIME)
        finally:
            self._ok2TestCon()

    def _ipControlChanged(self):
        ipaddress = self.lineEditSendIPFOPIDContr.text().strip(" ").split(".")
        try:
            if len(ipaddress) == 4:
                for a in ipaddress:
                    if 0 <= int(a) < 255:
                        pass
                    else:
                        self._isIPControlOk = False
                        print("{0} is Invalid in the IP address".format(a))
                        self.statusbar.showMessage("{0} is Invalid in the IP address".format(a), STATUSBAR_TIME)
                self._isIPControlOk = True
                print("'{0}'  VALID IP Address format".format(ipaddress))
                self.statusbar.showMessage("'{0}' VALID IP Address format".format(ipaddress), STATUSBAR_TIME)
            else:
                self._isIPControlOk = False
                print("'{0}' needs to be formated with 3 '.'s i.e 'X.X.X.X'".format(ipaddress))
        except Exception as exce:
            traceback.format_stack()
            self._isIPControlOk = False
            print("{0} NOT an IP Address".format(ipaddress))
            self.statusbar.showMessage("{0} is NOT an IP Address".format(ipaddress), STATUSBAR_TIME)
        finally:
            self._ok2Tune()
    #endregion

    #region PORT LINE EDIT CHECKER
    def _recPortChanged(self):
        portnum = self.lineEditRecievePort.text()
        try:

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
        portnum = self.lineEditSendPort.text()
        try:
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

    def _recvPortControlChanged(self):
        portnum = self.lineEditSendPortFOPIDContr.text()
        try:
            if MIN_PORT <= int(portnum) <= MAX_PORT:
                self._isPortControlOk = True
            else:
                self._isPortControlOk = False
                print("Port:{0} is out of range. '{1} <= port <= {2}'".format(portnum, MIN_PORT, MAX_PORT))
                self.statusbar.showMessage(
                    "Port:{0} is out of range. '{1} <= port <= {2}'".format(portnum, MIN_PORT, MAX_PORT),
                    STATUSBAR_TIME)
        except:
            self._isPortControlOk = False
            print("ERROR: {0} is NOT an integer".format(portnum))
            self.statusbar.showMessage("ERROR: {0} is NOT an integer".format(portnum), STATUSBAR_TIME)
        finally:
            self._ok2Tune()

    #endregion
    #endregion

    #region Button Clicked Functions
    def _addData(self):
        _loadData = newfofopdtguiclass()
        # _loadData.setFocus()
        # _loadData.foregroundRole()
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
        self._ok2TestCon()

    def _tune(self):
        try:
            fofopdtTune.ACTIVATE_TUNING = True
            fofopdtTune.DT = float(self.lineEditDt.text())   #very important to be first
            # fofopdtTune.NMAX = int(self.lineEditIter.text())                #Max allowed N for oustaloop approximation
            fofopdtTune.OPT_MAX_ITER = int(self.lineEditIter.text())        #Max iteration during tuning
            fofopdtTune.INPUT_MARGIN = float(self.lineEditGainMargin.text())

            #get oustaloop model
            oustalModel = Dict(dict(wb=10**float(self.lineEdit_StartFreq.text()), wh=10**float(self.lineEdit_StopFreq.text()),
                                       N=int(self.lineEditOrder.text()), Ts= 1/float(self.lineEditDt.text())))    #step = 1/samplerate
            #get fopidGuessModel
            fopidGuessModel = Dict(dict(Kp = float(self.lineEdit_Kp.text()), Ki = float(self.lineEdit_Ki.text()),
                                        Kd = float(self.lineEdit_Kd.text()),  lam = float(self.lineEdit_Lam.text()),
                                        mu = float(self.lineEdit_Mu.text())))
            #get tuning/design parameter
            tunninParams = Dict(dict(wc = 0.1, pm = float(self.lineEditPhaseMargin.text()), opt_norm = fofopdtTune.OPT_NORM))

            #get Current FOFOPDT Model
            currentModel = self.comboBoxFOFOPDTSYS.currentData()

            #initialize the tuningProcess #TODO: if possible in another thread os as a task
            fofopdtTune.mainFOFOPIDOPT(currentModel, fopidGuessModel, oustalModel, tunninParams)

            self.lineEdit_Kp.setText(str(fofopdtTune.und_fopid.Kp))
            self.lineEdit_Ki.setText(str(fofopdtTune.und_fopid.Ki))
            self.lineEdit_Kd.setText(str(fofopdtTune.und_fopid.Kd))
            self.lineEdit_Lam.setText(str(fofopdtTune.und_fopid.lam))
            self.lineEdit_Mu.setText(str(fofopdtTune.und_fopid.mu))

            self._ok2TestCon()
            fofopdtTune.ACTIVATE_TUNING = False
        except Exception as e:
            self._ok2TestCon()
            print("An exception occurred. Try using another limit/ initial guess settings")
            traceback.print_exc()
            self.statusbar.showMessage("An exception occurred. Try using another limit/ initial guess settings", STATUSBAR_TIME)

    def _updateFOPIDServer(self):
        try:
            confs = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

            # Set the parameters here to test and send them using the socket
            Kp = self.lineEdit_Kp.text()
            Ki = self.lineEdit_Ki.text()
            Kd = self.lineEdit_Kd.text()
            lam =self.lineEdit_Lam.text()
            mu = self.lineEdit_Mu.text()
            fopid = dict(Kp=Kp, Ki=Ki, Kd=Kd, lam=lam, mu=mu)
            controlfopid = dict(control = fopid)
            confs.sendto(bytes(repr(controlfopid),'utf-8'), (self.UDP_IP, self.UDP_PORT_CTRL))

        except Exception as e:
            traceback.print_exc()

    def startNetContrl(self):
        try:
            confs = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

            recvIp = self.lineEditRecieveIP.text()
            recvPort = self.lineEditRecievePort.text()
            sendPort = self.lineEditSendPort.text()
            sendIp = self.lineEditSendIP.text()
            simuTyme = self.lineEditTankSimTime.text()
            serverSettings = dict(recvIP=recvIp, recvPORT=recvPort, sendIP=sendIp, sendPORT=sendPort, time2Run=simuTyme)
            startControl = dict(start=serverSettings)

            if not self.isIPandPortExist(recvIp, recvPort):
                confs.sendto(bytes(repr(startControl), 'utf-8'), (self.UDP_IP, self.UDP_PORT_CTRL))
                self.pushButtonStartCom.setEnabled(False)
                self.pushButtonStopCom.setEnabled(True)
                self.pushButtonExitServer.setEnabled(True)

            else:
                confs.sendto(bytes(repr(startControl), 'utf-8'), (self.UDP_IP, self.UDP_PORT_CTRL))
                self.pushButtonStartCom.setEnabled(False)
                self.pushButtonStopCom.setEnabled(True)
                self.pushButtonExitServer.setEnabled(True)

            # print("Controller START Completed: Simulaion Time Set to {0}s",simuTyme)
        except Exception as e:
            traceback.print_exc()

            self.pushButtonStartCom.setEnabled(False)
            self.pushButtonStopCom.setEnabled(True)
            self.pushButtonExitServer.setEnabled(False)

    def stopNetContrl(self):
        try:
            confs = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            # Set the parameters here to test and send them using the socket
            stopControl = dict(stop = "")
            # NB! Parameter order is important! "s" means change start FOPID controller server parameters
            confs.sendto(bytes(repr(stopControl), 'utf-8'), (self.UDP_IP, self.UDP_PORT_CTRL))
        except Exception as e:
            traceback.print_exc()

        self.pushButtonStartCom.setEnabled(True)
        self.pushButtonStopCom.setEnabled(False)

    def _exitServer(self):
        try:
            confs = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            # Set the parameters here to test and send them using the socket
            exitServer = dict(exit = "")
            # NB! Parameter order is important! "s" means change start FOPID controller server parameters
            confs.sendto(bytes(repr(exitServer), 'utf-8'), (self.UDP_IP, self.UDP_PORT_CTRL))
            self.pushButtonStartCom.setEnabled(False)
            self.pushButtonStopCom.setEnabled(False)
            self.pushButtonTestControl.setEnabled(True)
            self.pushButtonExitServer.setEnabled(False)
            self.pushButtonToServer.setEnabled(False)

        except Exception as e:
            traceback.print_exc()

    def _TestControl(self):
        contrlPort = int(self.lineEditSendPortFOPIDContr.text())
        controlIP = self.lineEditSendIPFOPIDContr.text()
        try:
            if self.isIPandPortExist(controlIP,contrlPort):
                self.UDP_IP = self.lineEditSendIPFOPIDContr.text()
                self.UDP_PORT_CTRL = int(self.lineEditSendPortFOPIDContr.text())

                self.lineEditSendPortFOPIDContr.setEnabled(False)
                self.pushButtonTestControl.setEnabled(False)

                # self.lock.acquire()
                # print("Controller IP and Port Matched. Server FOUND")
                # self.lock.release()
            else:
                self.UDP_IP = 0
                self.UDP_PORT_CTRL = 0

                self.lineEditSendPortFOPIDContr.setEnabled(True)
                self.pushButtonTestControl.setEnabled(True)

            # self.lock.acquire()
            # print("ERROR!: Controller IP and Port Not Found. Run controlServer.py First")
            # self.lock.release()

        except Exception as e:
            traceback.print_exc()
            self.lineEditSendPortFOPIDContr.setEnabled(True)
            self.pushButtonTestControl.setEnabled(True)

        finally:
            self._ok2Tune()

    def updateTime(self):
        try:
            confs = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            # Set the parameters here to test and send them using the socket
            updatetimett = dict(time= dict(time2Run = self.lineEditTankSimTime.text()))
            # NB! Parameter order is important! "s" means change start FOPID controller server parameters
            confs.sendto(bytes(repr(updatetimett), 'utf-8'), (self.UDP_IP, self.UDP_PORT_CTRL))

        except:
            pass

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

        [z.setText(self.lineEditParamsValue.text()) for z in data]
        self._KpChanged()
        self._KiChanged()
        self._KdChanged()
        self._LamdaChanged()
        self._MuChanged()

    def closeEvent(self, event):
        reply = QMessageBox.question(self, "Exit?",
                                     "Are you sure you would like to 'Exit' this form?",
                                     QMessageBox.Yes | QMessageBox.No)
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    def isIPandPortExist(self, ip, port):
        contrlPort = int(port)
        controlIP = str(ip)
        open = None
        try:
            addressinfo = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
            addressinfo.bind((controlIP, contrlPort))
            addressinfo.shutdown(socket.SHUT_RDWR)
            addressinfo.close()
            return False

        except:
            return True
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
        self.sysnamecheck = self.gaincheck = self.delaycheck = self.timeconstantcheck = self.alphacheck = False
        self.lineEditSysName.setFocus()
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

# def fopidgui(lcok):
#     app = QApplication(sys.argv)
#     shareMem = fofopdtguiclass(lcok)
#     app.exec_()


# if __name__ == "__main__":
#     lock = Lock()
#     server = Process(target=controlServer.controlServerAutoStart, name="FOPIDServer", args = (controlServer.DURATION, lock,))
#     gui = Process(target= fopidgui, name="FOPIDGui", args= (lock,))
#
#     server.start()
#     gui.start()
#
#     server.join()
#     gui.join()
#
#     print("Gui and Server were 'Exited'")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    shareMem = fofopdtguiclass()
    app.exec_()