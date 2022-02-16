# ANTARES data browser for micro-ARPES spatial maps

This is an interactive software allowing to visualize and treat the micro-ARPES spatial maps from ANTARES beamline at syncrotron radiation facility SOLEIL (France).
The micro-ARPES data are stored in Nexus files with HDF5.

Main features:
- visualize the ARPES image from the selected point in spatial map,
- visualize integrated ARPES image from the selected area,
- convert ARPES data from real (E vs. angle) to momentum (E vs. k) space
- select the region of interest in real and momentum space and integrate the spatial map inside the selected region,
- save images (cvg, jpg etc.) and data (ARPES or spatial) in 'csv' format. 

![ARPES_GUI](https://user-images.githubusercontent.com/81705695/149383129-2ddd80c9-31a8-4cee-94b3-6241c02d6483.png)

Installation: 
1. the main files are in the folder 'src' 
2. the packages to be installed:
     - pip install h5py 
     - pip install scipy 
     - pip install pyqtgraph 
     - pip install PyQt5 
3. file to run 'ARPES_GUI_4.py'
      