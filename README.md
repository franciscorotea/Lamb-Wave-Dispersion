# Lamb Wave Dispersion

Lamb waves are a type of ultrasonic elastic waves that propagate guided between the two parallel surfaces of solid plates. Lamb waves have recently gained popularity from researchers and engineers in the non-destructive evaluation community for damage identification, due to their relatively low attenuation ratio, strong penetration capability, convenience of generation and collection, and high sensitivity to structural damage even of small size. 

Lamb waves propagate as two infinite sets of modes: ***symmetric modes***, in which the displacement is symmetrical about the midplane of the plate, and ***anti-symmetric modes***, with displacements anti-symmetric about the midplane. Another important characteristic of Lamb waves is their dispersive behavior: Lamb wave velocity depends on both the excitation frequency (f) and the thickness (d) of the plate combined in a frequency–thickness (f·d) product. 

This Python package presents tools to calculate and plot Lamb wave dispersion curves and particle displacement in traction-free, homogeneous and isotropic plates.

## Getting Started

The code is tested with Python 3.7. Next section provides the prerequisites to run the program.

### Prerequisites

The code is dependant on the following external libraries: Numpy, Scipy, Matplotlib. These can be installed with Python's inbuilt package management system, [pip](https://pip.pypa.io/en/stable/). See Python's tutorial on [installing packages](https://packaging.python.org/tutorials/installing-packages/#id17) for information about this issue. In short, the installation can be made as:

```
pip install numpy
pip install scipy
pip install matplotlib
```

## Usage:

First, you need to import the `Lamb` class from the `lamb` module, and create an instanciate it. For this example, we are going to use a 10 mm Aluminum plate.

```python
from lamb import Lamb

alum = Lamb(thickness=10, 
            nmodes_sym=5, 
            nmodes_antisym=5, 
            fd_max=10000, 
            vp_max=15000, 
            c_L=6420, 
            c_S=3040)
```

#### Parameters:
    
`thickness`: Thickness of the plate, in mm.  
`nmodes_sym`: Number of symmetric modes to calculate.  
`nmodes_antisym`: Number of antisymmetric modes to calculate.  
`fd_max`: Maximum value of frequency × thickness to calculate, in kHz × mm.  
`vp_max`: Maximum value of phase velocity to calculate, in m/s.  
`c_L`: Longitudinal wave velocity of the material, in m/s.  
`c_S`: Shear wave velocity of the material, in m/s.   
    
The following parameters are optional:
        
`c_R`: Rayleigh wave velocity of the material, in m/s. Defaults to None.  
`fd_points`: Number of frequency × thickness points. Defaults to 100.  
`vp_step`: Increment between phase velocity intervals. Defaults to 100.  
`material`: Name of the material being analyzed. Defaults to ''.  

### Methods

* ***Phase Velocity***

To generate a plot of phase velocity as a function of frequency × thickness, you can use:

```python
alum.plot_phase_velocity()
```

This method produces the following plot:

   ![alt text](https://i.imgur.com/yrHQj9L.png)
   
#### Parameters:

You can use the following optional parameters with this method:  

`modes`: Which family of modes to plot. Can be 'symmetric', 'antisymmetric' or 'both'. Defaults to 'both'.  
`cutoff_frequencies`: Add cutoff frequencies to the plot. Defaults to True.  
`material_velocities`: Add material velocities (longitudinal, shear and Rayleigh) to the plot. Defaults to True.  
`save_img`: Save the result image as png. Defaults to False.  
`sym_style`: A dictionary with [matplotlib kwargs](https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html) to modify the symmetric curves (to change color, linewidth, linestyle, etc.).  
`antisym_style`: A dictionary with [matplotlib kwargs](https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html) to modify the antisymmetric curves (to change color, linewidth, linestyle, etc.).  

* ***Group Velocity***

To generate a plot of group velocity as a function of frequency × thickness, you can use:

```python
alum.plot_group_velocity()
```

This method produces the following plot:

   ![alt text](https://i.imgur.com/HfcJfJI.png)
   
#### Parameters:

You can use the following optional parameters with this method:  

`modes`: Which family of modes to plot. Can be 'symmetric', 'antisymmetric' or 'both'. Defaults to 'both'.  
`cutoff_frequencies`: Add cutoff frequencies to the plot. Defaults to True.  
`save_img`: Save the result image as png. Defaults to False.  
`sym_style`: A dictionary with [matplotlib kwargs](https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html) to modify the symmetric curves (to change color, linewidth, linestyle, etc.).  
`antisym_style`: A dictionary with [matplotlib kwargs](https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html) to modify the antisymmetric curves (to change color, linewidth, linestyle, etc.).  

* ***Wave Number***

To generate a plot of wave number as a function of frequency × thickness, you can use:

```python
alum.plot_wave_number()
```

This method produces the following plot:

   ![alt text](https://i.imgur.com/uLitFVR.png)
   
#### Parameters:

You can use the following optional parameters with this method:  

`modes`: Which family of modes to plot. Can be 'symmetric', 'antisymmetric' or 'both'. Defaults to 'both'.   
`save_img`: Save the result image as png. Defaults to False.  
`sym_style`: A dictionary with [matplotlib kwargs](https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html) to modify the symmetric curves (to change color, linewidth, linestyle, etc.).  
`antisym_style`: A dictionary with [matplotlib kwargs](https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html) to modify the antisymmetric curves (to change color, linewidth, linestyle, etc.).  

* ***Wave Structure***

To generate a plot of the wave structure (i.e., the displacement profile across the thickness of the plate), you can use:

```python
alum.plot_wave_structure(mode='A0', nrows=3, ncols=2, fd=[500, 1000, 1500, 2000, 2500, 3000])
```

This method produces the following plot:

   ![alt text](https://i.imgur.com/F3fNEvL.png)
   
#### Parameters:

This method has to be used with the following parameters:  

`mode`: Mode to be analyzed. Can be "A0", "A1", "A2", ..., "An" or "S0", "S1", "S2", ..., "Sn", with 'n' being the order of the corresponding mode.  
`nrows`: Number of rows in the subplot.  
`ncols`: Number of columns in the subplot.  
`fd`: Array with the frequency × thickness values to analyze. The length of the array must be equal to `nrows` x `ncols`.  

The following parameters are optional:

`save_img`: Number of rows in the subplot. Defaults to False.  
`inplane_style`: A dictionary with [matplotlib kwargs](https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html) to modify the in-plane curves (to change color, linewidth, linestyle, etc.).  
`outofplane_style`: A dictionary with [matplotlib kwargs](https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html) to modify the out-of-plane curves (to change color, linewidth, linestyle, etc.).  

* ***Particle Displacement Field***

To generate an animation of the particle displacement field, you can use:

```python
alum.animate_displacement(mode='A0', fd=1000)
```

This method produces the following animation:

   ![alt text](https://thumbs.gfycat.com/OrderlyEnviousBoilweevil-size_restricted.gif)
   
#### Parameters:

This method has to be used with the following parameters:  

`mode`: Mode to be animated. Can be "A0", "A1", "A2", ..., "An" or "S0", "S1", "S2", ..., "Sn", with 'n' being the order of the corresponding mode.  
`fd`: Frequency × thickness value to animate.

The following parameters are optional:

`speed`: Delay between frames in milliseconds. It can be used to control the speed of the rotating vectors in the animation (a smaller value produces a faster animation). Defaults to 30.  
`save_gif`: Set to True if you want to save the result animation as a gif. Defaults to False.  
`save_video`: Choose a video format if you want to save the result animation as a video. Can be 'mp4', 'mov' or 'avi'. Defaults to False.  

***Note***: If you want to save the animation as a gif, you should install [ImageMagick](http://www.imagemagick.org/script/download.php) and specify the full path to magick.exe like this before using the `animate_displacement` method:

```python
Lamb.magick_path = 'C:/Program Files/ImageMagick-7.0.10-Q16/magick.exe'
```

If you want to save the animation as .mp4, .avi or .mov, you should specify the full path to the ffmpeg executable in ImageMagick installation folder:

```python
Lamb.ffmpeg_path = 'C:/Program Files/ImageMagick-7.0.10-Q16/ffmpeg.exe'
```    

If you are using some flavor of Unix, chances are ImageMagick is already installed on your computer.  

### Attributes

* ***Phase Velocity***

You can use the attributes `vp_sym` and `vp_antisym` to find the phase velocity at a particular `fd` value or an array of `fd` values. They are dictionaries with interpolators at each mode, where the keys are "A0", "A1", "A2", ..., "An" (for `vp_antisym`) and "S0", "S1", "S2", ..., "Sn" (for `vp_sym`), with 'n' being the order of the corresponding mode.  

For example, if you need the phase velocity for the S0 mode at 1000 kHz × mm, you can do:

```python
alum.vp_sym['S0'](1000)
```      

And this should return 5265.14 m/s. Always make sure that the fd values are within the valid range for the corresponding mode (i. e., above the cutoff frequency and below the `fd_max` you chose). Also, make sure the mode selected is within the selected `nmodes`. For example, if you chose `nmodes_sym = 5`, you can use 'S0', 'S1', 'S2', 'S3' or 'S4'.

* ***Group Velocity***

You can use the attributes `vg_sym` and `vg_antisym` to find the group velocity at a particular `fd` value or an array of `fd` values. They are dictionaries with interpolators at each mode, where the keys are "A0", "A1", "A2", ..., "An" (for `vg_antisym`) and "S0", "S1", "S2", ..., "Sn" (for `vg_sym`), with 'n' being the order of the corresponding mode.  

For example, if you need the group velocity for the A1 mode at 2000, 3000, and 4000 kHz × mm, you can do:

```python
alum.vg_antisym['A1']([2000,3000,4000])
```      

And this should return 3241.72, 3577.26, and 2486.33 m/s. Always make sure that the fd values are within the valid range for the corresponding mode (i. e., above the cutoff frequency and below the `fd_max` you chose). Also, make sure the mode selected is within the selected `nmodes`. For example, if you chose `nmodes_antisym = 5`, you can use 'A0', 'A1', 'A2', 'A3' or 'A4'.

* ***Wave Number***

You can use the attributes `k_sym` and `k_antisym` to find the wave number at a particular `fd` value or an array of `fd` values. They are dictionaries with interpolators at each mode, where the keys are "A0", "A1", "A2", ..., "An" (for `k_antisym`) and "S0", "S1", "S2", ..., "Sn" (for `k_sym`), with 'n' being the order of the corresponding mode.  

For example, if you need the wave number for the S3 mode at 8000 kHz × mm, you can do:

```python
alum.k_sym['S0'](8000)
```      

And this should return 726.38 m-1. Always make sure that the fd values are within the valid range for the corresponding mode (i. e., below the `fd_max` you chose). Also, make sure the mode selected is within the selected `nmodes`. For example, if you chose `nmodes_sym = 5`, you can use 'S0', 'S1', 'S2', 'S3' or 'S4'.

Run `example_code.py` for a ... continuara.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
