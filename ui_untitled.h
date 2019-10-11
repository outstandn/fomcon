/********************************************************************************
** Form generated from reading UI file 'untitledYa5920.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UNTITLEDYA5920_H
#define UNTITLEDYA5920_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_Dialog
{
public:
    QVBoxLayout *verticalLayout;
    QGroupBox *groupBox_10;
    QGridLayout *gridLayout_8;
    QComboBox *comboBoxFOTF_2;
    QPushButton *pushButton_AddFotf_2;
    QPushButton *pushButton_EditFOTF_2;
    QPushButton *pushBottonRefresh_2;
    QPushButton *pushButton_DeleteFOTF_2;
    QPushButton *pushButton_GetOustaloop_2;
    QGroupBox *groupBox_8;
    QVBoxLayout *verticalLayout_3;
    QGroupBox *groupBox_Freq_2;
    QGridLayout *gridLayout_6;
    QLabel *label_23;
    QLabel *label_24;
    QPushButton *bottonBodePlot_2;
    QLineEdit *lineEdit_LowerFreq_2;
    QLabel *label_25;
    QPushButton *pushButtonViewInConsole_2;
    QLineEdit *lineEdit_FreqDataPoints_2;
    QPushButton *pushButton_StabilityTest_2;
    QLineEdit *lineEdit_HigherFreq;
    QGroupBox *groupBox_9;
    QGridLayout *gridLayout_7;
    QLabel *label_26;
    QLabel *label_27;
    QLabel *label_28;
    QLineEdit *lineEdit_StartTime_2;
    QLineEdit *lineEdit_StepTime_2;
    QLineEdit *lineEdit_StopTime_2;
    QLabel *label_29;
    QLabel *label_TimeVector_2;
    QComboBox *comboBoxtimeType_2;
    QLineEdit *lineEdit_inputTime_2;
    QLineEdit *lineEdit_TimeVextorDisplay_2;
    QPushButton *buttoSimulate_2;

    void setupUi(QDialog *Dialog)
    {
        if (Dialog->objectName().isEmpty())
            Dialog->setObjectName(QStringLiteral("Dialog"));
        Dialog->resize(754, 557);
        verticalLayout = new QVBoxLayout(Dialog);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        groupBox_10 = new QGroupBox(Dialog);
        groupBox_10->setObjectName(QStringLiteral("groupBox_10"));
        groupBox_10->setMinimumSize(QSize(0, 0));
        groupBox_10->setMaximumSize(QSize(16777215, 16777215));
        QFont font;
        font.setPointSize(12);
        font.setBold(false);
        font.setItalic(false);
        font.setWeight(50);
        groupBox_10->setFont(font);
        gridLayout_8 = new QGridLayout(groupBox_10);
        gridLayout_8->setObjectName(QStringLiteral("gridLayout_8"));
        comboBoxFOTF_2 = new QComboBox(groupBox_10);
        comboBoxFOTF_2->setObjectName(QStringLiteral("comboBoxFOTF_2"));
        QFont font1;
        font1.setPointSize(12);
        comboBoxFOTF_2->setFont(font1);

        gridLayout_8->addWidget(comboBoxFOTF_2, 1, 0, 2, 3);

        pushButton_AddFotf_2 = new QPushButton(groupBox_10);
        pushButton_AddFotf_2->setObjectName(QStringLiteral("pushButton_AddFotf_2"));
        QFont font2;
        font2.setPointSize(12);
        font2.setBold(false);
        font2.setWeight(50);
        font2.setKerning(false);
        pushButton_AddFotf_2->setFont(font2);

        gridLayout_8->addWidget(pushButton_AddFotf_2, 0, 0, 1, 1);

        pushButton_EditFOTF_2 = new QPushButton(groupBox_10);
        pushButton_EditFOTF_2->setObjectName(QStringLiteral("pushButton_EditFOTF_2"));
        pushButton_EditFOTF_2->setFont(font2);

        gridLayout_8->addWidget(pushButton_EditFOTF_2, 0, 1, 1, 1);

        pushBottonRefresh_2 = new QPushButton(groupBox_10);
        pushBottonRefresh_2->setObjectName(QStringLiteral("pushBottonRefresh_2"));
        QFont font3;
        font3.setPointSize(12);
        font3.setBold(false);
        font3.setWeight(50);
        pushBottonRefresh_2->setFont(font3);

        gridLayout_8->addWidget(pushBottonRefresh_2, 3, 0, 1, 1);

        pushButton_DeleteFOTF_2 = new QPushButton(groupBox_10);
        pushButton_DeleteFOTF_2->setObjectName(QStringLiteral("pushButton_DeleteFOTF_2"));
        pushButton_DeleteFOTF_2->setFont(font2);

        gridLayout_8->addWidget(pushButton_DeleteFOTF_2, 0, 2, 1, 1);

        pushButton_GetOustaloop_2 = new QPushButton(groupBox_10);
        pushButton_GetOustaloop_2->setObjectName(QStringLiteral("pushButton_GetOustaloop_2"));
        pushButton_GetOustaloop_2->setFont(font2);

        gridLayout_8->addWidget(pushButton_GetOustaloop_2, 3, 2, 1, 1);


        verticalLayout->addWidget(groupBox_10);

        groupBox_8 = new QGroupBox(Dialog);
        groupBox_8->setObjectName(QStringLiteral("groupBox_8"));
        groupBox_8->setMinimumSize(QSize(0, 0));
        groupBox_8->setMaximumSize(QSize(16777215, 16777215));
        groupBox_8->setFont(font3);
        verticalLayout_3 = new QVBoxLayout(groupBox_8);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        groupBox_Freq_2 = new QGroupBox(groupBox_8);
        groupBox_Freq_2->setObjectName(QStringLiteral("groupBox_Freq_2"));
        QFont font4;
        font4.setPointSize(12);
        font4.setBold(false);
        font4.setItalic(false);
        font4.setWeight(50);
        font4.setKerning(false);
        groupBox_Freq_2->setFont(font4);
        gridLayout_6 = new QGridLayout(groupBox_Freq_2);
        gridLayout_6->setObjectName(QStringLiteral("gridLayout_6"));
        label_23 = new QLabel(groupBox_Freq_2);
        label_23->setObjectName(QStringLiteral("label_23"));
        label_23->setFont(font1);

        gridLayout_6->addWidget(label_23, 0, 1, 1, 1);

        label_24 = new QLabel(groupBox_Freq_2);
        label_24->setObjectName(QStringLiteral("label_24"));
        label_24->setFont(font1);
        label_24->setAlignment(Qt::AlignCenter);

        gridLayout_6->addWidget(label_24, 0, 0, 1, 1);

        bottonBodePlot_2 = new QPushButton(groupBox_Freq_2);
        bottonBodePlot_2->setObjectName(QStringLiteral("bottonBodePlot_2"));
        bottonBodePlot_2->setFont(font3);

        gridLayout_6->addWidget(bottonBodePlot_2, 3, 2, 1, 1);

        lineEdit_LowerFreq_2 = new QLineEdit(groupBox_Freq_2);
        lineEdit_LowerFreq_2->setObjectName(QStringLiteral("lineEdit_LowerFreq_2"));
        lineEdit_LowerFreq_2->setMaximumSize(QSize(16777215, 16777215));
        lineEdit_LowerFreq_2->setFont(font3);
        lineEdit_LowerFreq_2->setAlignment(Qt::AlignCenter);

        gridLayout_6->addWidget(lineEdit_LowerFreq_2, 1, 0, 1, 1);

        label_25 = new QLabel(groupBox_Freq_2);
        label_25->setObjectName(QStringLiteral("label_25"));
        label_25->setFont(font1);

        gridLayout_6->addWidget(label_25, 0, 2, 1, 1);

        pushButtonViewInConsole_2 = new QPushButton(groupBox_Freq_2);
        pushButtonViewInConsole_2->setObjectName(QStringLiteral("pushButtonViewInConsole_2"));
        pushButtonViewInConsole_2->setFont(font3);

        gridLayout_6->addWidget(pushButtonViewInConsole_2, 3, 0, 1, 1);

        lineEdit_FreqDataPoints_2 = new QLineEdit(groupBox_Freq_2);
        lineEdit_FreqDataPoints_2->setObjectName(QStringLiteral("lineEdit_FreqDataPoints_2"));
        lineEdit_FreqDataPoints_2->setMaximumSize(QSize(16777215, 16777215));
        lineEdit_FreqDataPoints_2->setFont(font3);
        lineEdit_FreqDataPoints_2->setAlignment(Qt::AlignCenter);

        gridLayout_6->addWidget(lineEdit_FreqDataPoints_2, 1, 2, 1, 1);

        pushButton_StabilityTest_2 = new QPushButton(groupBox_Freq_2);
        pushButton_StabilityTest_2->setObjectName(QStringLiteral("pushButton_StabilityTest_2"));
        pushButton_StabilityTest_2->setFont(font3);

        gridLayout_6->addWidget(pushButton_StabilityTest_2, 3, 1, 1, 1);

        lineEdit_HigherFreq = new QLineEdit(groupBox_Freq_2);
        lineEdit_HigherFreq->setObjectName(QStringLiteral("lineEdit_HigherFreq"));

        gridLayout_6->addWidget(lineEdit_HigherFreq, 1, 1, 1, 1);


        verticalLayout_3->addWidget(groupBox_Freq_2);

        groupBox_9 = new QGroupBox(groupBox_8);
        groupBox_9->setObjectName(QStringLiteral("groupBox_9"));
        QFont font5;
        font5.setPointSize(12);
        font5.setBold(false);
        font5.setItalic(false);
        font5.setUnderline(false);
        font5.setWeight(50);
        font5.setKerning(false);
        groupBox_9->setFont(font5);
        gridLayout_7 = new QGridLayout(groupBox_9);
        gridLayout_7->setObjectName(QStringLiteral("gridLayout_7"));
        label_26 = new QLabel(groupBox_9);
        label_26->setObjectName(QStringLiteral("label_26"));
        label_26->setFont(font1);
        label_26->setAlignment(Qt::AlignCenter);

        gridLayout_7->addWidget(label_26, 0, 0, 1, 1);

        label_27 = new QLabel(groupBox_9);
        label_27->setObjectName(QStringLiteral("label_27"));
        QFont font6;
        font6.setStrikeOut(false);
        label_27->setFont(font6);
        label_27->setAlignment(Qt::AlignCenter);

        gridLayout_7->addWidget(label_27, 0, 1, 1, 1);

        label_28 = new QLabel(groupBox_9);
        label_28->setObjectName(QStringLiteral("label_28"));
        label_28->setAlignment(Qt::AlignCenter);

        gridLayout_7->addWidget(label_28, 0, 2, 1, 1);

        lineEdit_StartTime_2 = new QLineEdit(groupBox_9);
        lineEdit_StartTime_2->setObjectName(QStringLiteral("lineEdit_StartTime_2"));
        lineEdit_StartTime_2->setMaximumSize(QSize(16777215, 16777215));
        lineEdit_StartTime_2->setFont(font3);
        lineEdit_StartTime_2->setAlignment(Qt::AlignCenter);

        gridLayout_7->addWidget(lineEdit_StartTime_2, 1, 0, 1, 1);

        lineEdit_StepTime_2 = new QLineEdit(groupBox_9);
        lineEdit_StepTime_2->setObjectName(QStringLiteral("lineEdit_StepTime_2"));
        lineEdit_StepTime_2->setMaximumSize(QSize(16777215, 16777215));
        lineEdit_StepTime_2->setFont(font3);
        lineEdit_StepTime_2->setAlignment(Qt::AlignCenter);

        gridLayout_7->addWidget(lineEdit_StepTime_2, 1, 1, 1, 1);

        lineEdit_StopTime_2 = new QLineEdit(groupBox_9);
        lineEdit_StopTime_2->setObjectName(QStringLiteral("lineEdit_StopTime_2"));
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(lineEdit_StopTime_2->sizePolicy().hasHeightForWidth());
        lineEdit_StopTime_2->setSizePolicy(sizePolicy);
        lineEdit_StopTime_2->setMaximumSize(QSize(16777215, 16777215));
        lineEdit_StopTime_2->setFont(font3);
        lineEdit_StopTime_2->setAlignment(Qt::AlignCenter);

        gridLayout_7->addWidget(lineEdit_StopTime_2, 1, 2, 1, 1);

        label_29 = new QLabel(groupBox_9);
        label_29->setObjectName(QStringLiteral("label_29"));
        label_29->setFont(font1);
        label_29->setAlignment(Qt::AlignCenter);

        gridLayout_7->addWidget(label_29, 2, 0, 1, 1);

        label_TimeVector_2 = new QLabel(groupBox_9);
        label_TimeVector_2->setObjectName(QStringLiteral("label_TimeVector_2"));
        label_TimeVector_2->setMaximumSize(QSize(16777215, 16777215));
        label_TimeVector_2->setFont(font1);
        label_TimeVector_2->setAlignment(Qt::AlignCenter);

        gridLayout_7->addWidget(label_TimeVector_2, 2, 1, 1, 1);

        comboBoxtimeType_2 = new QComboBox(groupBox_9);
        comboBoxtimeType_2->setObjectName(QStringLiteral("comboBoxtimeType_2"));
        sizePolicy.setHeightForWidth(comboBoxtimeType_2->sizePolicy().hasHeightForWidth());
        comboBoxtimeType_2->setSizePolicy(sizePolicy);
        comboBoxtimeType_2->setFont(font3);
        comboBoxtimeType_2->setEditable(false);
        comboBoxtimeType_2->setMaxVisibleItems(2);
        comboBoxtimeType_2->setModelColumn(0);

        gridLayout_7->addWidget(comboBoxtimeType_2, 2, 2, 1, 1);

        lineEdit_inputTime_2 = new QLineEdit(groupBox_9);
        lineEdit_inputTime_2->setObjectName(QStringLiteral("lineEdit_inputTime_2"));
        lineEdit_inputTime_2->setMaximumSize(QSize(16777215, 16777215));
        lineEdit_inputTime_2->setFont(font3);
        lineEdit_inputTime_2->setAlignment(Qt::AlignCenter);

        gridLayout_7->addWidget(lineEdit_inputTime_2, 3, 0, 1, 1);

        lineEdit_TimeVextorDisplay_2 = new QLineEdit(groupBox_9);
        lineEdit_TimeVextorDisplay_2->setObjectName(QStringLiteral("lineEdit_TimeVextorDisplay_2"));
        lineEdit_TimeVextorDisplay_2->setMaximumSize(QSize(16777215, 16777215));
        lineEdit_TimeVextorDisplay_2->setAlignment(Qt::AlignCenter);
        lineEdit_TimeVextorDisplay_2->setReadOnly(true);

        gridLayout_7->addWidget(lineEdit_TimeVextorDisplay_2, 3, 1, 1, 1);

        buttoSimulate_2 = new QPushButton(groupBox_9);
        buttoSimulate_2->setObjectName(QStringLiteral("buttoSimulate_2"));
        sizePolicy.setHeightForWidth(buttoSimulate_2->sizePolicy().hasHeightForWidth());
        buttoSimulate_2->setSizePolicy(sizePolicy);
        buttoSimulate_2->setFont(font3);

        gridLayout_7->addWidget(buttoSimulate_2, 3, 2, 1, 1);


        verticalLayout_3->addWidget(groupBox_9);


        verticalLayout->addWidget(groupBox_8);

#ifndef QT_NO_SHORTCUT
#endif // QT_NO_SHORTCUT

        retranslateUi(Dialog);

        comboBoxtimeType_2->setCurrentIndex(-1);


        QMetaObject::connectSlotsByName(Dialog);
    } // setupUi

    void retranslateUi(QDialog *Dialog)
    {
        Dialog->setWindowTitle(QApplication::translate("Dialog", "Dialog", Q_NULLPTR));
        groupBox_10->setTitle(QApplication::translate("Dialog", "FO Transfer Functions", Q_NULLPTR));
        pushButton_AddFotf_2->setText(QApplication::translate("Dialog", "Add", Q_NULLPTR));
        pushButton_EditFOTF_2->setText(QApplication::translate("Dialog", "Edit", Q_NULLPTR));
        pushBottonRefresh_2->setText(QApplication::translate("Dialog", "Refresh List", Q_NULLPTR));
        pushButton_DeleteFOTF_2->setText(QApplication::translate("Dialog", "Delete", Q_NULLPTR));
        pushButton_GetOustaloop_2->setText(QApplication::translate("Dialog", "Get Oustaloop Model", Q_NULLPTR));
        groupBox_8->setTitle(QApplication::translate("Dialog", "System Analysis", Q_NULLPTR));
        groupBox_Freq_2->setTitle(QApplication::translate("Dialog", "Frequency Domain", Q_NULLPTR));
        label_23->setText(QApplication::translate("Dialog", "Higher Freq (rad/s):", Q_NULLPTR));
        label_24->setText(QApplication::translate("Dialog", "Lower freq (rad/s):", Q_NULLPTR));
        bottonBodePlot_2->setText(QApplication::translate("Dialog", "Bode Plot", Q_NULLPTR));
        lineEdit_LowerFreq_2->setText(QApplication::translate("Dialog", "0.001", Q_NULLPTR));
        label_25->setText(QApplication::translate("Dialog", "Data points (int):", Q_NULLPTR));
        pushButtonViewInConsole_2->setText(QApplication::translate("Dialog", "View in console", Q_NULLPTR));
        lineEdit_FreqDataPoints_2->setText(QApplication::translate("Dialog", "5000", Q_NULLPTR));
        pushButton_StabilityTest_2->setText(QApplication::translate("Dialog", "Stability test", Q_NULLPTR));
        groupBox_9->setTitle(QApplication::translate("Dialog", "Time Domain", Q_NULLPTR));
        label_26->setText(QApplication::translate("Dialog", "Start (s):", Q_NULLPTR));
        label_27->setText(QApplication::translate("Dialog", "Step (s):", Q_NULLPTR));
        label_28->setText(QApplication::translate("Dialog", "Stop (s):", Q_NULLPTR));
        lineEdit_StartTime_2->setText(QApplication::translate("Dialog", "0", Q_NULLPTR));
        lineEdit_StepTime_2->setText(QApplication::translate("Dialog", "0.1", Q_NULLPTR));
        lineEdit_StopTime_2->setText(QApplication::translate("Dialog", "30", Q_NULLPTR));
        label_29->setText(QApplication::translate("Dialog", "Input:", Q_NULLPTR));
        label_TimeVector_2->setText(QApplication::translate("Dialog", "TimeVector:", Q_NULLPTR));
        comboBoxtimeType_2->setCurrentText(QString());
        lineEdit_inputTime_2->setText(QApplication::translate("Dialog", "1", Q_NULLPTR));
        buttoSimulate_2->setText(QApplication::translate("Dialog", "Simulate", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class Dialog: public Ui_Dialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UNTITLEDYA5920_H
