# Lamb Wave Dispersion

Lamb waves are a type of ultrasonic elastic waves that propagate guided between the two parallel surfaces of solid plates. Lamb waves have recently gained popularity from researchers and engineers in the non-destructive evaluation community for damage identification, due to their relatively low attenuation ratio, strong penetration capability, convenience of generation and collection, and high sensitivity to structural damage even of small size. 

Lamb waves propagate as two infinite sets of modes: ***symmetric modes***, in which the displacement is symmetrical about the midplane of the plate, and ***anti-symmetric modes***, with displacements anti-symmetric about the midplane. Another important characteristic of Lamb waves is their dispersive behavior: Lamb wave velocity depends on both the excitation frequency (f) and the thickness (d) of the plate combined in a frequency–thickness (f·d) product. 

This Python package presents tools to calculate and plot Lamb wave dispersion curves and particle displacement in traction-free, homogeneous and isotropic plates.

# 3D Audio Panner

This code provides a Python implementation of a 3D Audio Panner with a GUI, using the head-related impulse responses (HRIR) recorded in the [CIPIC HRTF Database](https://www.ece.ucdavis.edu/cipic/spatial-sound/hrtf-data/) by the CIPIC Interface Laboratory at UC Davis.

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

1. Run `example_code.py`.

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
