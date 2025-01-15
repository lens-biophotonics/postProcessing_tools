## Scripts/Tools

### crop_zoom.py
The **crop_zoom.py** is used for cropping the OME-TIF file to the desired x dimensions and then 
zoom it back to the original X shape. A **crop_zoom.ipynb** jupyter notebook is provided for testing. The **crop_zoom.py** takes **crop_zoom-p.yaml** config yaml file as a argument. One can run the script by the following commonds:
```
python3 crop_zoom.py -p crop_zoom-p.yaml
```

> **Note**: To run the python scripts one has to create a python venv or conda with the right python libraries installed. 
