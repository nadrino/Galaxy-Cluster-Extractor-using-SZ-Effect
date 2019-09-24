# README #

### What is this repository for? ###

* Galaxy cluster extraction using SZ-effect with Planck data.

### Dependencies ###

* Install Python (worked on 2.7.11) with apt-get (linux) or homebrew (mac / linux)
* Install Python packages :

```
pip install numpy
pip install matplotlib
pip install healpy
pip install scipy
```

* The following files are required:

```
HFI_SkyMap_100_2048_R2.02_full.fits
HFI_SkyMap_143_2048_R2.02_full.fits
HFI_SkyMap_217_2048_R2.02_full.fits
HFI_SkyMap_353_2048_R2.02_full.fits
HFI_SkyMap_545_2048_R2.02_full.fits
HFI_SkyMap_857_2048_R2.02_full.fits
```

You can download them here: [https://irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/HFI_Products.html](https://irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/HFI_Products.html)

Or go into your input folder (you will need to set it up in the next section) and type the commands: 

```
wget https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/HFI_SkyMap_100_2048_R2.02_full.fits
wget https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/HFI_SkyMap_143_2048_R2.02_full.fits
wget https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/HFI_SkyMap_217_2048_R2.02_full.fits
wget https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/HFI_SkyMap_353_2048_R2.02_full.fits
wget https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/HFI_SkyMap_545_2048_R2.02_full.fits
wget https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/HFI_SkyMap_857_2048_R2.02_full.fits
```

### How do I get set up? ###

* Go to the folder you want the 'ias_project' folder to be.
* Use ```git clone https://bitbucket.org/aazlm/ias_project.git```
* Go into  ```ias_project/src/ ``` and open ```setup.py```
* You have to indicate the paths to the required file ```in_files_path``` and the path for the output files ```out_files_path``` for example:
```python
in_files_path = './input/'
out_files_path = './'
# Don't forget the '/' at the end of the path
```
* The program is ready to run !

### How to run the cluster finder? ###

* Convolute the maps with ```./convolute_degrade_maps.py```
* Convert all maps in MJy ```./convert_map_units.py```
* Create the maps ```./create_mask.py```
* Create the SZ map ```./ILC.py```
* Find clusters ```./find_clusters.py```
* All the results can be found in the folder ```out_files_path``` you just specified.
