import sys
import h5py
import numpy as np
from decimal import *
from scipy import interpolate

import pyqtgraph as pg
from PyQt5 import QtCore, uic, QtGui
from PyQt5.QtWidgets import (
    QMainWindow,
    QApplication,
    QLabel,
    QErrorMessage,
    QMessageBox,
    QFileDialog,
    QSlider,
    QHBoxLayout,
    QLCDNumber,
    )

from PyQt5.QtCore import (QDir)
from main4 import Ui_MainWindow

# This code is used if main4.ui file is converted to main4.py file 
class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.UI = Ui_MainWindow()
        self.UI.setupUi(self)
        self.init_buttons()
        self.check_box_angle()
        self.check_box_k()
        self.init_graph_spatial()
        self.init_graph_ARPES()

# This code is used if main4.ui file is NOT converted to main4.py file 
# (by the comand #pyuic5 main4.ui > main4.py):
#     
# class MainApp(QMainWindow):
#    def __init__(self):
#        QMainWindow.__init__(self)
#        self.UI = uic.loadUi('main4.ui', self)
#        self.init_buttons()
#        self.check_box_angle()
#        self.check_box_k()
#        self.init_graph_spatial()
#        self.init_graph_ARPES()

    def check_box_angle(self):
        ''' Connect check box 'Set zero' to the set_zero function (chosing the zero angle)
        '''
        self.UI.checkBox_angle.stateChanged.connect(self.set_zero)   

    def check_box_k(self):
        ''' Connect check box 'k-scale' to the update_ARPES_to_k function
        '''
        self.UI.checkBox_k.stateChanged.connect(self.update_ARPES_to_k)

    def slider(self):
        ''' Determine values of Slider regarding to the chosen zero_angle
        '''
        self.UI.Slider.setTickPosition(QSlider.TicksBelow)
        self.UI.Slider.valueChanged.connect(self.zero_angle) 

    def init_buttons(self):
        ''' Initialize load and safe buttons
        '''
        self.UI.LoadButton.released.connect(self.load_data)
        self.UI.SaveImageButton.released.connect(self.save_image)
        self.UI.SaveARPESButton.released.connect(self.save_arpes)
 
    def save_image(self):
        ''' Save the spatial image.
            ROI - region of interest. 
            if ROI is determined (update_roi is running), pass the ROI values to the string, othervise - integrated image
        '''
        try:                                                             
            A1 = self.angle_scale[self.Angle_range_items[0]]           
            A2 = self.angle_scale[self.Angle_range_items[1]]
            E1 = self.energy_scale[self.Energy_range_items[0]]
            E2 = self.energy_scale[self.Energy_range_items[1]]
          
            # Pass angles (determined by ROI) to the string_save_image
            string_save_image = f"Ang_{A1:.1f}_{A2:.1f}__En_{E1:.2f}_{E2:.2f}"
           
            # If checkBox_k is checked, pass the k-values to the string_save_image
            if self.UI.checkBox_k.isChecked():
                k1 = self.k_scale[self.k_roi_index_0]           
                k2 = self.k_scale[self.k_roi_index_1]
                string_save_image = f"k_{k1:.2f}_{k2:.2f}__En_{E1:.2f}_{E2:.2f}" 
        
        # If ROI is not determined (update_roi not running)        
        except:
            string_save_image = "integrated_image"
       
        # Chose the path for saving and save
        path, _ = QFileDialog.getSaveFileName(self, 'Save Spatial Image', QDir.homePath() + 
            "/Spatial_" + string_save_image + ".csv", "CSV Files(*.csv *.txt)") 
        if path:
            try: 
                np.savetxt(path, self.roi_spatial_image(), delimiter=",")
            
            # If spacial image is not yet selected (update_roi), save the integrated image
            except:
                np.savetxt(path, self.integrated_image, delimiter=",")
    
    def save_arpes(self): 
        ''' Save the ARPES spatial image together with X, Y coordinates
        '''
        path, _ = QFileDialog.getSaveFileName(self, 'Save ARPES', QDir.homePath() + 
            f"/ARPES_X_{self.X_coord}__Y_{self.Y_coord}.csv", "CSV Files(*.csv *.txt)") 
        if path:
            np.savetxt(path, self.ARPES_image(), delimiter=",")

    def load_data(self):
        ''' Load NXS file
            if some file was opened before, we should delite the lines, ROI, and unckeck the checkBox_k
        '''
        try:
            self.UI.checkBox_k.setChecked(False)
            self.UI.checkBox_angle.setChecked(False)
            self.plotSpatial.removeItem(self.vLine)
            self.plotSpatial.removeItem(self.hLine)
            self.plotARPES.removeItem(self.roi)
        except:
            pass
       
        # Suggest only NXS files for opening
        file = QFileDialog.getOpenFileName(self, "Open NXS File", "*.nxs")
        self.path = file[0]
        self.f = h5py.File(self.path, "r")
        self.ARPES_data = np.array(self.f['salsaentry_1/scan_data/data_12'])
        
        # Determine the integrated image and show it in the left window  
        self.integrated_image = np.rot90(np.sum(self.ARPES_data, axis=(2,3)))
        self.update_spatial_image(self.integrated_image)
       
        # Determine data parameters (X, Y, energy, angle scales, offsets), show ARPES image from the middle position (X/Y), add lines and ROI, set real scales
        # param (0) is the real zero angle, which initially = 0
        self.determine_data_parameters(0)

    def update_spatial_image(self, image):
        self.imgSpatial.setImage(image)

    def determine_data_parameters(self, offset):
        ''' Determine parameters of the spatial and ARPES images

        :param offset: the chosen offset in zero angle to convert the ARPES image from angular to the k-space 

                       Data Parameters
                        ----------
                X_scale - X scale of the spatial image 
                Y_scale - Y scale of the spatial image 
                precision_energy - decimal places in the energy scale
                precision_angle -  decimal places in the angle scale
                energy_offset - lowest energy
                energy_step - step in enrgy scale
                energy_scale - energy scale of the ARPES image
                angle_offset - smallest angle            
                angle_step - step in angle scale
                angle_scale - angle scale of the ARPES image
                        ----------
        The scale_img function is used to set the real scales in spatial and ARPES windows;
        add lines, which crossing point determines location of the ARPES image;
        show the ARPES image determined by the X/Y position of the lines' crossing point.
        '''
        self.X_scale = np.array(self.f['salsaentry_1/scan_data/actuator_1_1'][0]) 
        self.Y_scale = np.array(self.f['salsaentry_1/scan_data/actuator_2_1'])    

        precision_energy = 7 
        precision_angle = 5

        self.energy_offset = np.around(self.f.get('salsaentry_1/scan_data/data_04')[0][0], precision_energy)
        self.energy_step = np.around(self.f.get('salsaentry_1/scan_data/data_05')[0][0], precision_energy)
        energy_end = np.around(self.energy_offset + self.energy_step * (self.ARPES_data.shape[3] - 1), precision_energy)
        self.energy_scale = np.linspace(self.energy_offset, energy_end, self.ARPES_data.shape[3]) 

        self.angle_offset = np.around(self.f.get('salsaentry_1/scan_data/data_07')[0][0], precision_angle) 
       
        # if check box 'Set zero' is selected, we should shift the angle offset according to the chosen zero angle:
        if self.UI.checkBox_angle.isChecked():
            self.angle_offset -= offset
        self.angle_step = np.around(self.f.get('salsaentry_1/scan_data/data_08')[0][0], precision_angle)
        angle_end = np.around(self.angle_offset + self.angle_step * (self.ARPES_data.shape[2] - 1), precision_angle)
        self.angle_scale = np.linspace(self.angle_offset, angle_end, self.ARPES_data.shape[2])
              
        # set real scales to the spatial image:
        rect = self.scale_img(self.X_scale, self.Y_scale)
        self.imgSpatial.setRect(rect)   
        
        # add horizontal and vertical lines in the middle of the spatial image
        self.add_lines(self.X_scale[int(self.X_scale.shape[0]/2)], self.Y_scale[int(self.Y_scale.shape[0]/2)])
        
        # show the ARPES image from the X/Y point
        self.update_ARPES_image()
       
        # set real scale (in angles or in k, if checkBox_k is selected) to the ARPES image:
        if self.UI.checkBox_k.isChecked():
            rect = self.scale_img(self.k_scale, self.energy_scale) 
        else:
            rect = self.scale_img(self.angle_scale, self.energy_scale)
        self.imgARPES.setRect(rect)
        
        # after real scale is shown, add ROI,
        # self.angle_step - means that the ROI should be added in the angle scale (if checkBox_k is not checked)
        self.add_roi(self.angle_step)

    def update_ARPES_image(self):
        ''' Update ARPES image when the lines are moved, and puts actual coordinates in the title, 
            self.ARPES_image - is the ARPES image data determined from the data file (self.ARPES_data)
        '''
        self.imgARPES.setImage(self.ARPES_image())
        self.plotSpatial.setTitle(f'X: {self.X_coord:.1f}, Y: {self.Y_coord:.1f}')

    def ARPES_image(self):
        ''' Determine data of the ARPES image from the chosed X/Y coordinates
        :return: ARPES image at selected coordinates
        '''
        self.X_coord = self.find_position(self.X_scale, self.vLine.getPos()[0]) 
        self.Y_coord = self.find_position(self.Y_scale, self.hLine.getPos()[1]) 
        X_index = self.X_scale.shape[0] - 1 - np.where(self.X_scale == self.X_coord)[0][0]
        Y_index = np.where(self.Y_scale == self.Y_coord)[0][0] 
        
        # if check box 'k-scale' is selected, convert angles to k
        if self.UI.checkBox_k.isChecked():
            image_at_coord = self.ang_to_k(self.ARPES_data[Y_index, X_index, :, :], self.energy_scale, self.angle_scale)
        else: 
            image_at_coord = self.ARPES_data[Y_index, X_index, :, :]
        
        return image_at_coord
        
    def add_roi(self, angle_step):
        ''' Add ROI to the ARPES image - from zero to the middle in energy and angular/k-scale;
            if checkBox_k is checked, ROI is added in the k-scale
        :param angle_step: step in the angular scale (if 'k-selected', ROI is displaied in the k-scale)
        '''

        if self.UI.checkBox_k.isChecked():
            angle_step = self.convertAngletoK(angle_step,self.energy_scale[0])
        self.roi = pg.ROI([0, self.energy_offset + self.energy_step], 
            [angle_step * self.ARPES_data.shape[2]/2, self.energy_step * self.ARPES_data.shape[3]/2])
        self.roi.addScaleHandle([0.0, 0.0], [1, 1])
        
        # Signal is emitted any time when the position of ROI changes (ROI is dragged by the user)
        self.roi.sigRegionChanged.connect(self.update_roi) 
       
        # Always have ROi on top of the image
        self.roi.setZValue(10) 
        self.plotARPES.addItem(self.roi)

    def add_lines(self, pos_X, pos_Y):
        ''' Add horizontal (hLine) and vertical (vLine) lines to the spatial image
            and connects the chosen position with the update_ARPES_image method.
        :param pos_X: position along the X scale of spatial image
        :param pos_Y: position along the Y scale of spatial image
        '''
        self.vLine = pg.InfiniteLine(pos=pos_X, angle=90, movable=True, bounds=[self.X_scale[0], self.X_scale[-1]])
        self.hLine = pg.InfiniteLine(pos=pos_Y, angle=0, movable=True, bounds=[self.Y_scale[0], self.Y_scale[-1]])
        self.plotSpatial.addItem(self.vLine, ignoreBounds=True)
        self.plotSpatial.addItem(self.hLine, ignoreBounds=True)
        self.vLine.sigPositionChanged.connect(self.update_ARPES_image) 
        self.hLine.sigPositionChanged.connect(self.update_ARPES_image)

    def init_graph_spatial(self): 
        ''' Initialize spatial image in the MainGraphView window.
        '''
        self.winI = self.UI.MainGraphView   
        
        # Add 2D Image View
        self.plotSpatial = self.winI.addPlot(title="Spatial image")  # plot area (ViewBox + axes) for displaying the image
        self.imgSpatial = pg.ImageItem() # item for displaying image data: either a 2D numpy array or a 3D array
        self.plotSpatial.addItem(self.imgSpatial)
       
        # Add histrogramm
        hist = pg.HistogramLUTItem()
        hist.setImageItem(self.imgSpatial)
        hist.setMaximumWidth(150)
        hist.gradient.loadPreset('viridis')
        self.winI.addItem(hist)

    def init_graph_ARPES(self):
        ''' Initialize ARPES image in the MainGraphView_2 window.
        '''
        self.winA= self.UI.MainGraphView_2
        
        # Add 2D Image View
        self.plotARPES = self.winA.addPlot(title="ARPES image")
        self.imgARPES = pg.ImageItem() 
        self.plotARPES.addItem(self.imgARPES)   
        
        # Add histrogramm
        hist = pg.HistogramLUTItem()
        hist.setImageItem(self.imgARPES)
        hist.setMaximumWidth(150)
        hist.gradient.loadPreset('flame')
        self.winA.addItem(hist) 

    def update_roi(self):
        ''' Update angle and energy for the selected ROI (ROI is selected from the ARPES image)
        '''
        selected = self.roi.getArrayRegion(self.ARPES_image(), self.imgARPES)
        shape_selected = np.shape(selected)
        handle_pos = self.roi.pos()  

        # set the k_scale ROI and the k-limits if k-scale is selected
        if self.UI.checkBox_k.isChecked():

            if handle_pos[0] > (self.k_scale[-1]+(self.k_scale[1]-self.k_scale[0])/2): 
                self.k_roi_index_0 = self.k_scale.shape[0]-1
            else:
                self.k_roi_index_0  = np.where(self.k_scale == self.find_position(self.k_scale, handle_pos[0]))[0][0]

            handle_pos_k_right = handle_pos[0] + shape_selected[0]*(self.k_scale[1]-self.k_scale[0])

            if handle_pos_k_right > (self.k_scale[-1] + (self.k_scale[1]-self.k_scale[0])/2): 
                self.k_roi_index_1 = self.k_scale.shape[0]-1
            else:
                self.k_roi_index_1 = np.where(self.k_scale == self.find_position(self.k_scale, handle_pos_k_right))[0][0]
            # convert handle position to the angle
            handle_pos[0] = self.convertKtoAngle(handle_pos[0], handle_pos[1])
        
        # set the angle_scale ROI and the limits
        if handle_pos[0] > (self.angle_scale[-1]+self.angle_step/2): 
            angle_roi_index_0 = self.angle_scale.shape[0]-1
        else:
            angle_roi_index_0 = np.where(self.angle_scale == self.find_position(self.angle_scale, handle_pos[0]))[0][0]
     
        if self.UI.checkBox_k.isChecked():
            handle_pos_angle_right = self.convertKtoAngle(handle_pos_k_right, handle_pos[1])
        else:
            handle_pos_angle_right = handle_pos[0] + shape_selected[0]*self.angle_step
        
        if handle_pos_angle_right > (self.angle_scale[-1] + self.angle_step/2): 
            angle_roi_index_1 = self.angle_scale.shape[0]-1
        else:
            angle_roi_index_1 = np.where(self.angle_scale == self.find_position(self.angle_scale, handle_pos_angle_right))[0][0]
        
        self.Angle_range_items = (angle_roi_index_0, angle_roi_index_1)
        
        # set the energy_scale ROI and the limits
        if handle_pos[1] > (self.energy_scale[-1]+self.energy_step/2): 
            energy_roi_index_0 = self.energy_scale.shape[0]-1
        else:
            energy_roi_index_0 = np.where(self.energy_scale == self.find_position(self.energy_scale, handle_pos[1]))[0][0]
       
        handle_pos_energy_top = handle_pos[1] + shape_selected[1]*self.energy_step
        
        if handle_pos_energy_top > (self.energy_scale[-1]+self.energy_step/2): 
            energy_roi_index_1 = self.energy_scale.shape[0]-1
        else:
            energy_roi_index_1 = np.where(self.energy_scale == self.find_position(self.energy_scale, handle_pos_energy_top))[0][0]
    
        self.Energy_range_items = (energy_roi_index_0, energy_roi_index_1)
        
        # diasplay spatial image integrated in the ROI
        image = self.roi_spatial_image()
        self.update_spatial_image(image)
        self.plotARPES.setTitle(f'A ({self.angle_scale[self.Angle_range_items[0]]:.1f} : {self.angle_scale[self.Angle_range_items[1]]:.1f}), \
                                  E ({self.energy_scale[self.Energy_range_items[0]]:.2f} : {self.energy_scale[self.Energy_range_items[1]]:.2f})')

        if self.UI.checkBox_k.isChecked():
            self.plotARPES.setTitle(f'k ({self.k_scale[self.k_roi_index_0]:.2f} : {self.k_scale[self.k_roi_index_1]:.2f}), \
                                  E ({self.energy_scale[self.Energy_range_items[0]]:.2f} : {self.energy_scale[self.Energy_range_items[1]]:.2f})')

    def roi_spatial_image(self):
        ''' Determine the spatial image based on the selected ROI
        :return: spatial image
        '''
        lowerA = self.Angle_range_items[0]
        upperA = self.Angle_range_items[1]
        lowerE = self.Energy_range_items[0]
        upperE = self.Energy_range_items[1]
        selection = self.ARPES_data[:, :, lowerA:upperA, lowerE:upperE]
        image = np.rot90(np.sum(selection, axis=(2,3)))
        return image
      
    def scale_img(self, X_scale, Y_scale):
        ''' Transform image to real x and y scales by setting the rectangle (rect) in the proper scale
        :param X_scale: spatial X scale
        :param Y_scale: spatial Y scale
        :return: rectangle in the X/Y scalae
        '''
        left = X_scale[0] - (X_scale[1]-X_scale[0])/2
        bottom = Y_scale[0] - (Y_scale[1]-Y_scale[0])/2
        # this is needed to set the precice scale, and avoid the problem with the way how the float number is stored in binary
        width = np.around(X_scale[-1] - X_scale[0], - Decimal(str(X_scale[0])).as_tuple().exponent) + (X_scale[1]-X_scale[0])
        height = np.around(Y_scale[-1] - Y_scale[0], - Decimal(str(Y_scale[0])).as_tuple().exponent) + (Y_scale[1]-Y_scale[0])
        rect = QtCore.QRectF(left , bottom , width, height)
        return rect

    def find_position(self, array, value):
        ''' Find position in the array (e.g. energy sclae) closest to the chosen value
        '''
        idx = np.argmax(array + (array[1] - array[0])/2 - value > 0)
        return array[idx]

    def update_ARPES_to_k(self):
        ''' Update the ARPES image if the checkbox 'k-scale' is selected
        ''' 
        # remove ROI and horizontal and vertical lines    
        self.plotSpatial.removeItem(self.vLine)
        self.plotSpatial.removeItem(self.hLine)
        self.plotARPES.removeItem(self.roi)
        
        # determine new parameters with respect to the chosen zero angle and update ROI with respect to the new parameters
        # if 'k-scale' is selected, the ARPES image is converyted to k-space
        self.determine_data_parameters(self.angle_0)
        self.update_roi()
        
        # remove ZeroLine (which choses the zero angle) if it is already shown
        try:
            self.plotARPES.removeItem(self.ZeroLine)
        except:
            pass

    def set_zero(self, state): 
        ''' When the check box 'Set zero' is selected, initiate slider and chose zero angle by zero_angle function 
        :param state: state of the check box
        '''      
        if state == QtCore.Qt.Checked:
            self.slider() 
            self.zero_angle()
        else:
            self.angle_0 = 0
            self.UI.lcd.display(self.angle_0)
            self.UI.Slider.setSliderPosition(self.angle_0)
            self.update_ARPES_to_k()
            self.UI.lcd.setPalette(QApplication.style().standardPalette())
            self.plotARPES.removeItem(self.ZeroLine)

    def zero_angle(self):
        ''' Select zero angle (self.angle_0) by the slider
        ''' 
        self.angle_0 = self.UI.Slider.value()/10
        palette = self.UI.lcd.palette()
        minimum = self.angle_scale[0]
        maximum = self.angle_scale[-1]
        self.UI.Slider.setMinimum(int(round(minimum*10)))
        self.UI.Slider.setMaximum(int(round(maximum*10)))
        palette.setColor(palette.Background, QtGui.QColor(204, 255, 204))
        palette.setColor(palette.WindowText, QtGui.QColor(0,0,0))
        self.UI.lcd.setSegmentStyle(QtGui.QLCDNumber.Flat)
        self.UI.lcd.setPalette(palette)
        self.UI.lcd.display(self.angle_0)
        # remove ZeroLine (which choses the zero angle) if it is already shown
        try:
            self.plotARPES.removeItem(self.ZeroLine)
        except:
            pass
        # draw ZeroLine as a vertical line at chosen angle_0
        self.ZeroLine = pg.InfiniteLine(pos=self.angle_0, angle=90, movable=False)
        self.plotARPES.addItem(self.ZeroLine, ignoreBounds=True)

    def ang_to_k(self, array, energy_scale, angle_scale):
        ''' Convert angles (angle_acale) to the k-space (k_scale) for chosen ARPES image (array)
        :param array: ARPES image in angular space
        :param energy_scale: float
        :param angle_scale: float
        :return: ARPES image in k-space
        ''' 
        # interpolate angle-energy data 
        f = interpolate.interp2d(angle_scale, energy_scale, array.T, fill_value=0.0)                           
        Emax = energy_scale[-1]
        ang_min = angle_scale[0] 
        ang_max = angle_scale[-1]
        # determine k-scale
        kmin    = self.convertAngletoK(ang_min, Emax)
        kmax    = self.convertAngletoK(ang_max, Emax)
        self.k_scale = np.linspace(kmin, kmax, angle_scale.shape[0])
        kdelta = self.k_scale[1] - self.k_scale[0]
        # generate temporary ARPES image, which will be rescaled in the k-scale
        temp = np.zeros((len(energy_scale), len(angle_scale)))
        
        # convert each energy slice from angles to the k-space
        for count, E in enumerate(energy_scale):
            cutoff_lower = self.convertAngletoK(ang_min, energy_scale[count]) # = 0.5124*np.sqrt(energy_scale[count])*np.sin(np.pi/180*(ang_min))
            cutoff_upper = self.convertAngletoK(ang_max, energy_scale[count]) # = 0.5124*np.sqrt(energy_scale[count])*np.sin(np.pi/180*(ang_max))
            k_range = np.where(np.logical_and(cutoff_lower<self.k_scale, self.k_scale<cutoff_upper))[0]
            lower_k = int(k_range[0])
            upper_k = int(k_range[-1])
            k_scale_temp = self.k_scale[k_range]
            ang = self.convertKtoAngle(k_scale_temp, E)
            temp[count, lower_k:upper_k+1] = f(ang, E)

        return temp.T

    def convertKtoAngle(self, k, Energy):
        ''' Convert a k-posion to an angle
        :param k: k-posion
        :param Energy: energy
        :return: angle
        '''
        # const = sqrt(2 * m_e)/h-bar
        const = 0.5124
        try:
            angle = 180/np.pi*np.arcsin(k / (const*np.sqrt(Energy)))
        except Exception as e:
            print(e)
        return angle

    def convertAngletoK(self, ang, Energy):
        ''' Convert an angle to a k-position
        :param k: angle
        :param Energy: energy
        :return: momentum k
        ''' 
        # const = sqrt(2 * m_e)/h-bar
        const = 0.5124
        k = const*np.sqrt(Energy)*np.sin(np.pi/180*(ang))
        return k

if __name__ == "__main__":
    qApp = QApplication(sys.argv)
    qApp.setFont(QtGui.QFont('Helvetica', 12))
    aw = MainWindow()
    aw.show()
    sys.exit(qApp.exec_())

