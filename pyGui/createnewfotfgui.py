# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\Set\Documents\CompSys\4th Semeter THESIS\Code\fomcon\guipyQtui\createnewfotfgui.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_dialogCreateNewFOTF(object):
    def setupUi(self, dialogCreateNewFOTF):
        dialogCreateNewFOTF.setObjectName("dialogCreateNewFOTF")
        dialogCreateNewFOTF.setWindowModality(QtCore.Qt.WindowModal)
        dialogCreateNewFOTF.resize(320, 240)
        dialogCreateNewFOTF.setMinimumSize(QtCore.QSize(320, 240))
        dialogCreateNewFOTF.setMaximumSize(QtCore.QSize(480, 360))
        dialogCreateNewFOTF.setModal(True)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(dialogCreateNewFOTF)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setContentsMargins(5, 5, 5, 5)
        self.verticalLayout.setSpacing(5)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(dialogCreateNewFOTF)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.lineEditSysName = QtWidgets.QLineEdit(dialogCreateNewFOTF)
        self.lineEditSysName.setObjectName("lineEditSysName")
        self.verticalLayout.addWidget(self.lineEditSysName)
        self.label_4 = QtWidgets.QLabel(dialogCreateNewFOTF)
        self.label_4.setObjectName("label_4")
        self.verticalLayout.addWidget(self.label_4)
        self.lineEdit_ZeroPoly = QtWidgets.QLineEdit(dialogCreateNewFOTF)
        self.lineEdit_ZeroPoly.setObjectName("lineEdit_ZeroPoly")
        self.verticalLayout.addWidget(self.lineEdit_ZeroPoly)
        self.label_3 = QtWidgets.QLabel(dialogCreateNewFOTF)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)
        self.lineEdit_PolePoly = QtWidgets.QLineEdit(dialogCreateNewFOTF)
        self.lineEdit_PolePoly.setObjectName("lineEdit_PolePoly")
        self.verticalLayout.addWidget(self.lineEdit_PolePoly)
        self.label_2 = QtWidgets.QLabel(dialogCreateNewFOTF)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.lineEdit_DelayText = QtWidgets.QLineEdit(dialogCreateNewFOTF)
        self.lineEdit_DelayText.setObjectName("lineEdit_DelayText")
        self.verticalLayout.addWidget(self.lineEdit_DelayText)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setContentsMargins(5, 5, 5, 5)
        self.horizontalLayout.setSpacing(5)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.pushButtonOK = QtWidgets.QPushButton(dialogCreateNewFOTF)
        self.pushButtonOK.setEnabled(False)
        self.pushButtonOK.setObjectName("pushButtonOK")
        self.horizontalLayout.addWidget(self.pushButtonOK)
        self.pushButtonCancel = QtWidgets.QPushButton(dialogCreateNewFOTF)
        self.pushButtonCancel.setObjectName("pushButtonCancel")
        self.horizontalLayout.addWidget(self.pushButtonCancel)
        self.verticalLayout_2.addLayout(self.horizontalLayout)

        self.retranslateUi(dialogCreateNewFOTF)
        QtCore.QMetaObject.connectSlotsByName(dialogCreateNewFOTF)

    def retranslateUi(self, dialogCreateNewFOTF):
        _translate = QtCore.QCoreApplication.translate
        dialogCreateNewFOTF.setWindowTitle(_translate("dialogCreateNewFOTF", "Create new FOTF"))
        self.label.setText(_translate("dialogCreateNewFOTF", "System Name: (e.g. System_1):"))
        self.label_4.setText(_translate("dialogCreateNewFOTF", "Zero Polynomial (e.g. s+11):"))
        self.label_3.setText(_translate("dialogCreateNewFOTF", "Pole Polynomial (e.g. -s^1.5-1):"))
        self.label_2.setText(_translate("dialogCreateNewFOTF", "Delay [sec]"))
        self.lineEdit_DelayText.setText(_translate("dialogCreateNewFOTF", "0"))
        self.pushButtonOK.setText(_translate("dialogCreateNewFOTF", "OK"))
        self.pushButtonCancel.setText(_translate("dialogCreateNewFOTF", "Cancel"))

