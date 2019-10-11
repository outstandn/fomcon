import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.uic import *
from PyQt5.QtWidgets import QMessageBox

from fotf import *
import fotftidgui
import  fotfviewergui

class FomconForm(QMainWindow, fotftidgui.Ui_MainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        fotftidgui.Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.setWindowIcon(QIcon('index.png'))
        self.show()


    def Exit(self):
        reply = QMessageBox.question(self, "Exit?", "Would you like to exit?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            sys.exit()



if __name__ == "__main__":
    app = QApplication(sys.argv)
    fomcon = FomconForm()
    app.exec_()