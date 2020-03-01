# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'D:\CompSys\4th Semeter THESIS\Code\fomcon\loaddatagui.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_LoadDataForm(object):
    def setupUi(self, LoadDataForm):
        LoadDataForm.setObjectName("LoadDataForm")
        LoadDataForm.setWindowModality(QtCore.Qt.WindowModal)
        LoadDataForm.resize(402, 134)
        LoadDataForm.setMinimumSize(QtCore.QSize(0, 0))
        LoadDataForm.setMaximumSize(QtCore.QSize(16777215, 16777215))
        LoadDataForm.setModal(True)
        self.gridLayout = QtWidgets.QGridLayout(LoadDataForm)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(LoadDataForm)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 2)
        self.lineEditSysName = QtWidgets.QLineEdit(LoadDataForm)
        self.lineEditSysName.setObjectName("lineEditSysName")
        self.gridLayout.addWidget(self.lineEditSysName, 1, 0, 1, 2)
        self.label_4 = QtWidgets.QLabel(LoadDataForm)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 2, 0, 1, 1)
        self.lineEditDataPath = QtWidgets.QLineEdit(LoadDataForm)
        self.lineEditDataPath.setReadOnly(True)
        self.lineEditDataPath.setObjectName("lineEditDataPath")
        self.gridLayout.addWidget(self.lineEditDataPath, 3, 0, 1, 2)
        self.pushButtonBrowse = QtWidgets.QPushButton(LoadDataForm)
        self.pushButtonBrowse.setMinimumSize(QtCore.QSize(125, 23))
        self.pushButtonBrowse.setMaximumSize(QtCore.QSize(125, 23))
        self.pushButtonBrowse.setObjectName("pushButtonBrowse")
        self.gridLayout.addWidget(self.pushButtonBrowse, 3, 2, 1, 1)
        self.pushButtonOK = QtWidgets.QPushButton(LoadDataForm)
        self.pushButtonOK.setEnabled(False)
        self.pushButtonOK.setObjectName("pushButtonOK")
        self.gridLayout.addWidget(self.pushButtonOK, 4, 1, 1, 1)
        self.pushButtonCancel = QtWidgets.QPushButton(LoadDataForm)
        self.pushButtonCancel.setObjectName("pushButtonCancel")
        self.gridLayout.addWidget(self.pushButtonCancel, 4, 2, 1, 1)

        self.retranslateUi(LoadDataForm)
        QtCore.QMetaObject.connectSlotsByName(LoadDataForm)

    def retranslateUi(self, LoadDataForm):
        _translate = QtCore.QCoreApplication.translate
        LoadDataForm.setWindowTitle(_translate("LoadDataForm", "Load Data"))
        self.label.setText(_translate("LoadDataForm", "Data Name: (e.g. Data_1):"))
        self.label_4.setText(_translate("LoadDataForm", "Data Path:"))
        self.pushButtonBrowse.setText(_translate("LoadDataForm", "Browse"))
        self.pushButtonOK.setText(_translate("LoadDataForm", "OK"))
        self.pushButtonCancel.setText(_translate("LoadDataForm", "Cancel"))

