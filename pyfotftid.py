import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.uic import *
from PyQt5.QtWidgets import QMessageBox

from fotf import *
from fomconoptimize import  *
import fotftidgui
import  fotfviewergui

class fotftidguiclas(QMainWindow, fotftidgui.Ui_MainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        fotftidgui.Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.setWindowIcon(QIcon('index.png'))
        self.reloadAllFOTransFunc()
        self.show()


    def reloadAllFOTransFunc(self):
        #Startup Config
        for i in optMethod:
            self.comboBoxOptTypeMethod.addItem(str(i), i)

        for i in optAlgo:
            self.comboBoxAlgorithm.addItem(str(i), i)

        for i in optFix:
            self.comboBoxOptFix.addItem(str(i), i)

        self.comboBoxPolesOrZeros.addItems(['Zero','Pole'])
        self.lineEditCoefLimitLower.setEnabled(False)
        self.lineEditCoefLimitUpper.setEnabled(False)
        self.lineEditExpLimitLower.setEnabled(False)
        self.lineEditExpLimitUpper.setEnabled(False)

    def Exit(self):
        reply = QMessageBox.question(self, "Exit?", "Would you like to exit?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            sys.exit()



if __name__ == "__main__":
    app = QApplication(sys.argv)
    fomcon = fotftidguiclas()
    app.exec_()