# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'D:\CompSys\4th Semeter THESIS\Code\fomcon\trimdatagui.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_TrimDataForm(object):
    def setupUi(self, TrimDataForm):
        TrimDataForm.setObjectName("TrimDataForm")
        TrimDataForm.setWindowModality(QtCore.Qt.WindowModal)
        TrimDataForm.resize(174, 131)
        TrimDataForm.setMinimumSize(QtCore.QSize(0, 0))
        TrimDataForm.setMaximumSize(QtCore.QSize(16777215, 16777215))
        TrimDataForm.setModal(True)
        self.gridLayout = QtWidgets.QGridLayout(TrimDataForm)
        self.gridLayout.setObjectName("gridLayout")
        self.lineEditT1 = QtWidgets.QLineEdit(TrimDataForm)
        self.lineEditT1.setReadOnly(False)
        self.lineEditT1.setObjectName("lineEditT1")
        self.gridLayout.addWidget(self.lineEditT1, 3, 0, 1, 2)
        self.label_4 = QtWidgets.QLabel(TrimDataForm)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 2, 0, 1, 2)
        self.label = QtWidgets.QLabel(TrimDataForm)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 2)
        self.lineEditNewName = QtWidgets.QLineEdit(TrimDataForm)
        self.lineEditNewName.setObjectName("lineEditNewName")
        self.gridLayout.addWidget(self.lineEditNewName, 1, 0, 1, 2)
        self.pushButtonOK = QtWidgets.QPushButton(TrimDataForm)
        self.pushButtonOK.setEnabled(False)
        self.pushButtonOK.setObjectName("pushButtonOK")
        self.gridLayout.addWidget(self.pushButtonOK, 4, 0, 1, 1)
        self.pushButtonCancel = QtWidgets.QPushButton(TrimDataForm)
        self.pushButtonCancel.setObjectName("pushButtonCancel")
        self.gridLayout.addWidget(self.pushButtonCancel, 4, 1, 1, 1)

        self.retranslateUi(TrimDataForm)
        QtCore.QMetaObject.connectSlotsByName(TrimDataForm)

    def retranslateUi(self, TrimDataForm):
        _translate = QtCore.QCoreApplication.translate
        TrimDataForm.setWindowTitle(_translate("TrimDataForm", "Trim Data"))
        self.label_4.setText(_translate("TrimDataForm", "Delay Ends at T1[s]"))
        self.label.setText(_translate("TrimDataForm", "New Data Name: (e.g. Data_1):"))
        self.pushButtonOK.setText(_translate("TrimDataForm", "OK"))
        self.pushButtonCancel.setText(_translate("TrimDataForm", "Cancel"))

