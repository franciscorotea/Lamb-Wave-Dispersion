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

First, you need to import the `Lamb` class from the `lamb` module, and create an instance of this class.

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

Run `example_code.py` for a ... continuara.

2. Select a subject from the CIPIC database. You should select a subject with similar anthropometric measurements as yourself for the best experience.

      ***Note:** Due to storage limitations, the repository has only 4 subjects of the database to choose from. The full database is ~170MB and has 45 subjects. It can be [downloaded for free](https://www.ece.ucdavis.edu/cipic/spatial-sound/hrtf-data/) at the CIPIC webpage. Make sure to download the MATLAB version of the database. In order to make it work, you should simply replace the folder `CIPIC_hrtf_database` with the one you downloaded.*

   ![alt text](https://i.imgur.com/wgbHujh.png)

3. Press `Play` to start playing the default sound file. Make sure to be using headphones as your audio output. Move the Azimuth and        Elevation sliders to position the sound in the 3D space. You can load your own audio file in File/Load audio file. Also, there are      other sound samples in the folder resources/sound.
   
   **IMPORTANT: For now, the only working format is a mono WAV file at 44100 Hz sample rate and 16 bit depth.**
   
   You can save the file at the specified pair of Azimuth/Elevation in `File/Save audio file`.
   Lastly, you can choose to use a crossover in order not to spatialize low frequencies, since low frequencies are non-directional in      nature. Go to `Settings/Change cutoff frequency` to set the desired frequency. By default, crossover is set at 200 Hz.

   ![alt text](https://i.imgur.com/xmcz00n.png)

## Implementation details

Before sound arrives to the auditory system, it is filtered by the diffraction and reflection properties of the head, pinna, and torso. This information is captured in the head related transfer function (HRTF), a pair of functions (one for each ear) that characterizes how an ear receives a sound from a point in space. The HRTF is highly dependent on the location of the sound source relative to the listener, which is the main reason we are able to locate the sound source. A pair of HRTFs for two ears can be used to synthesize a binaural sound that seems to come from a particular point in space.

The CIPIC database provides head related impulse responses (HRIR), that is, the inverse Fourier transform of the HRTF. The process of positioning a sound source in a virtual space using HRIRs consists on the convolution of the mono signal with the HRIR for left and right ear:

<a href="https://www.codecogs.com/eqnedit.php?latex=x_{L,&space;R}(t)&space;=&space;h_{L,&space;R}(t)&space;*&space;x(t)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x_{L,&space;R}(t)&space;=&space;h_{L,&space;R}(t)&space;*&space;x(t)" title="x_{L, R}(t) = h_{L, R}(t) * x(t)" /></a>

![alt text](https://i.imgur.com/WCnl0mG.png)

### HRIR Interpolation

The CIPIC database has HRIRs for a finite set of points in space. Nevertheless, it would be desirable to pan the audio source in space 
smoothly, without audible jumps from one point to another. Therefore, an interpolation must be made for the points in space in which there are no HRIRs recorded. For this panner, the interpolation approach outlined by H. Gamper [in this paper](https://asa.scitation.org/doi/full/10.1121/1.4828983) (available for free) is used.

In short, the method consists on performing a triangulation of the measured set of points (every combination of azimuth, elevation). These points become then vertices for every triangle in the triangulation. Lastly, the interpolated HRIR at any point X inside a triangle can be represented as the weighted sum of the HRIRs measured at each vertex of the triangle, using its barycentric coordinates as interpolation weights.

### Real-time convolution

It is also desirable that filtering is performed in real time, so that the audio changes as the user move the azimuth/elevation sliders. This procedure is not straightforward: since the length of the convolution between signals with length L and M would be L+M-1, it is not possible to simply concatenate the output blocks. Another issue that must be taken into account is the aliasing introduced by the cyclic FFT convolution.

To overcome these issues, the [overlap-save method](https://en.wikipedia.org/wiki/Overlap%E2%80%93save_method) is implemented. In short, this method consists on breaking the input audio signal into chunks of size L, transform the chunks into the frequency domain with the FFT and multiply it by the impulse response's DFT (i.e. convolution in time domain), transform back to the time domain and lop on the last L samples from the resulting L+M-1 chunk.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
