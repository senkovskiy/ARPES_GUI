# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'main5.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1124, 557)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setBold(False)
        font.setWeight(50)
        font.setKerning(True)
        MainWindow.setFont(font)
        MainWindow.setStyleSheet("")
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setBold(False)
        font.setItalic(False)
        font.setUnderline(False)
        font.setWeight(50)
        font.setKerning(True)
        self.centralwidget.setFont(font)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setBold(False)
        font.setItalic(False)
        font.setUnderline(False)
        font.setWeight(50)
        font.setKerning(False)
        self.groupBox.setFont(font)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.checkBox_angle = QtWidgets.QCheckBox(self.groupBox)
        self.checkBox_angle.setObjectName("checkBox_angle")
        self.gridLayout_2.addWidget(self.checkBox_angle, 0, 0, 1, 2)
        self.lcd = QtWidgets.QLCDNumber(self.groupBox)
        self.lcd.setEnabled(True)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setBold(False)
        font.setItalic(False)
        font.setUnderline(False)
        font.setWeight(50)
        font.setKerning(True)
        self.lcd.setFont(font)
        self.lcd.setStatusTip("")
        self.lcd.setAutoFillBackground(True)
        self.lcd.setStyleSheet("")
        self.lcd.setSmallDecimalPoint(False)
        self.lcd.setObjectName("lcd")
        self.gridLayout_2.addWidget(self.lcd, 1, 0, 1, 2)
        self.Slider = QtWidgets.QSlider(self.groupBox)
        self.Slider.setOrientation(QtCore.Qt.Horizontal)
        self.Slider.setObjectName("Slider")
        self.gridLayout_2.addWidget(self.Slider, 2, 0, 1, 2)
        self.label = QtWidgets.QLabel(self.groupBox)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 3, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        self.label_2.setObjectName("label_2")
        self.gridLayout_2.addWidget(self.label_2, 3, 1, 1, 1)
        self.checkBox_k = QtWidgets.QCheckBox(self.groupBox)
        self.checkBox_k.setObjectName("checkBox_k")
        self.gridLayout_2.addWidget(self.checkBox_k, 4, 0, 1, 1)
        self.gridLayout.addWidget(self.groupBox, 1, 0, 1, 1)
        self.MainGraphView_2 = GraphicsLayoutWidget(self.centralwidget)
        self.MainGraphView_2.setMouseTracking(False)
        self.MainGraphView_2.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.MainGraphView_2.setObjectName("MainGraphView_2")
        self.gridLayout.addWidget(self.MainGraphView_2, 0, 2, 5, 1)
        self.groupBox_3 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_3.setObjectName("groupBox_3")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.groupBox_3)
        self.verticalLayout.setObjectName("verticalLayout")
        self.LoadButton = QtWidgets.QPushButton(self.groupBox_3)
        self.LoadButton.setObjectName("LoadButton")
        self.verticalLayout.addWidget(self.LoadButton)
        self.SaveImageButton = QtWidgets.QPushButton(self.groupBox_3)
        self.SaveImageButton.setObjectName("SaveImageButton")
        self.verticalLayout.addWidget(self.SaveImageButton)
        self.SaveARPESButton = QtWidgets.QPushButton(self.groupBox_3)
        self.SaveARPESButton.setObjectName("SaveARPESButton")
        self.verticalLayout.addWidget(self.SaveARPESButton)
        self.gridLayout.addWidget(self.groupBox_3, 0, 0, 1, 1)
        self.MainGraphView = GraphicsLayoutWidget(self.centralwidget)
        self.MainGraphView.setAutoFillBackground(False)
        self.MainGraphView.setInteractive(True)
        self.MainGraphView.setObjectName("MainGraphView")
        self.gridLayout.addWidget(self.MainGraphView, 0, 1, 5, 1)
        self.groupBox_2 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_2.setObjectName("groupBox_2")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.groupBox_2)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.label_3 = QtWidgets.QLabel(self.groupBox_2)
        self.label_3.setObjectName("label_3")
        self.gridLayout_3.addWidget(self.label_3, 2, 0, 1, 1)
        self.YSpinBox = QtWidgets.QSpinBox(self.groupBox_2)
        self.YSpinBox.setMinimum(1)
        self.YSpinBox.setProperty("value", 1)
        self.YSpinBox.setObjectName("YSpinBox")
        self.gridLayout_3.addWidget(self.YSpinBox, 3, 1, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.groupBox_2)
        self.label_4.setObjectName("label_4")
        self.gridLayout_3.addWidget(self.label_4, 3, 0, 1, 1)
        self.IntegrateCheckBox = QtWidgets.QCheckBox(self.groupBox_2)
        self.IntegrateCheckBox.setObjectName("IntegrateCheckBox")
        self.gridLayout_3.addWidget(self.IntegrateCheckBox, 1, 1, 1, 1)
        self.XSpinBox = QtWidgets.QSpinBox(self.groupBox_2)
        self.XSpinBox.setMinimum(1)
        self.XSpinBox.setProperty("value", 1)
        self.XSpinBox.setObjectName("XSpinBox")
        self.gridLayout_3.addWidget(self.XSpinBox, 2, 1, 1, 1)
        self.LockAutoRangeBox = QtWidgets.QCheckBox(self.groupBox_2)
        self.LockAutoRangeBox.setObjectName("LockAutoRangeBox")
        self.gridLayout_3.addWidget(self.LockAutoRangeBox, 0, 1, 1, 1)
        self.gridLayout.addWidget(self.groupBox_2, 2, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1124, 19))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "ANTARES Data Browser"))
        self.groupBox.setTitle(_translate("MainWindow", "k-space data"))
        self.checkBox_angle.setText(_translate("MainWindow", "Set zero:"))
        self.label.setText(_translate("MainWindow", "min"))
        self.label_2.setText(_translate("MainWindow", "max"))
        self.checkBox_k.setText(_translate("MainWindow", "k-space"))
        self.groupBox_3.setTitle(_translate("MainWindow", "Load / Safe data"))
        self.LoadButton.setText(_translate("MainWindow", "Load"))
        self.SaveImageButton.setText(_translate("MainWindow", "Save Image"))
        self.SaveARPESButton.setText(_translate("MainWindow", "Save ARPES"))
        self.groupBox_2.setTitle(_translate("MainWindow", "GroupBox"))
        self.label_3.setText(_translate("MainWindow", "X Size"))
        self.label_4.setText(_translate("MainWindow", "Y Size"))
        self.IntegrateCheckBox.setText(_translate("MainWindow", "Integrate"))
        self.LockAutoRangeBox.setText(_translate("MainWindow", "LockAutoRange"))
from pyqtgraph import GraphicsLayoutWidget
