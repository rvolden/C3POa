# C3POa
Computational pipeline for calling consensi on R2C2 nanopore data.

### Dependencies ###
- Python 3.6
- NumPy 1.13.3
- poa v1.0.0 Revision: 1.2.2.9
- EMBOSS water: watHerON v8
- minimap2 2.7-r654
- racon

To install dependencies, use setup.sh.  
setup.sh will download and make all of the packages that you need to run C3POa (except Python and NumPy).
```
chmod +x setup.sh
./setup.sh
```

To install NumPy, you can go [here](https://scipy.org/install.html).  
Otherwise, you can use your computer's package manager (apt-get, dnf, brew, etc.) to install.  
Pip3 is another option for NumPy installation.  
Example:
```
sudo dnf install python3-numpy
```
or
```
pip3 install numpy
```
