"""A module with tools to calculate and plot Lamb wave dispersion 
curves.

Usage:
    
First, you need to create an instance of the Lamb class:
    
    mat = Lamb(thickness, nmodes_sym, nmodes_antisym, fd_max, vp_max, 
               c_L, c_S [, c_R=None][, fd_points=100][, vp_step=100]
               [, material=''])

Then, you can use this instance with the following methods:
    
    plot_phase_velocity(modes, cutoff_frequencies, material_velocities,
                        save_img, sym_style, antisym_style):
        Plot phase velocity as a function of frequency × thickness.
    plot_group_velocity(modes, cutoff_frequencies, save_img, sym_style,
                        antisym_style):
        Plot group velocity as a function of frequency × thickness.
    plot_wave_number(modes, save_img, sym_style, antisym_style):
        Plot wavenumber as a function of frequency × thickness.
    plot_wave_structure(mode, nrows, ncols, fd, save_img, inplane_style,
                        outofplane_style):
        Plot particle displacement across the thickness of the plate.
    animate_displacement(mode, fd, speed, save_gif, save_video):
        Generate an animation of the displacement vector field.
    save_results()
        Save all results to a txt file.

You can also use the following attributes:
    
    vp_sym:
        Phase velocity interpolators for symmetric modes.
    vg_sym:
        Group velocity interpolators for symmetric modes.
    k_sym:
        Wavenumber interpolators for symmetric modes.
    vp_antisym:
        Phase velocity interpolators for antisymmetric modes.
    vg_antisym:
        Group velocity interpolators for antisymmetric modes.
    k_antisym:
       Wavenumber interpolators for antisymmetric modes.

For example, if you need the phase velocity for the S0 mode at 1000 
kHz × mm, you can do:
    
    mat.vp_sym['S0'](1000)

You can also use a `np.array` instead of a single fd value. Always make
sure that the fd values are within the valid range for the corresponding 
mode (i. e., above the cutoff frequency and below the fd_max you chose).
Also, make sure the mode selected is within the selected `nmodes`. For 
example, if you chose `nmodes_sym = 4`, you can use 'S0', 'S1', 'S2' or 
'S3'.

For information about the equations implemented, please refer to:

Rose, J. L., Ultrasonic Guided Waves in Solid Media, Chapter 6: Waves in 
Plates, Cambridge University Press, 2014.
        
Graff, K. F., Wave Motion in Elastic Solids, Chapter 8: Wave Propagation 
in Plates and Rods, Dover Publications, 1975.

Author:         Francisco Rotea
                (Buenos Aires, Argentina)
Repository:     https://github.com/franciscorotea
Email:          francisco.rotea@gmail.com

"""

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation

import scipy.optimize

from .plot_utils import add_plot, add_cutoff_freqs, add_velocities
from .utils import interpolate, correct_instability, write_txt, find_max

class Lamb:
    """A class used to calculate and plot Lamb wave dispersion curves 
    for traction-free, homogeneous and isotropic plates. It also allows 
    to generate an animation of the displacement vector field.
    
    Methods
    -------
    plot_phase_velocity(modes, cutoff_frequencies, material_velocities,
                        save_img, sym_style, antisym_style):
        Plot phase velocity as a function of frequency × thickness.
    plot_group_velocity(modes, cutoff_frequencies, save_img, sym_style,
                        antisym_style):
        Plot group velocity as a function of frequency × thickness.
    plot_wave_number(modes, save_img, sym_style, antisym_style):
        Plot wavenumber as a function of frequency × thickness.
    plot_wave_structure(mode, nrows, ncols, fd, save_img, inplane_style,
                        outofplane_style):
        Plot particle displacement across the thickness of the plate.
    animate_displacement(mode, fd, speed, save_gif, save_video):
        Generate an animation of the displacement vector field.
    save_results()
        Save all results to a txt file.

    Attributes
    ----------
    vp_sym:
        Dictionary with phase velocity interpolators for symmetric 
        modes.
    vg_sym:
        Dictionary with group velocity interpolators for symmetric 
        modes.
    k_sym:
        Dictionary with wavenumber interpolators for symmetric 
        modes.
    vp_antisym:
        Dictionary with phase velocity interpolators for antisymmetric 
        modes.
    vg_antisym:
        Dictionary with group velocity interpolators for antisymmetric 
        modes.
    k_antisym:
       Dictionary with wavenumber interpolators for antisymmetric 
       modes.
       
    """
    
    # If you want to save the animation as a gif, you should install 
    # ImageMagick from http://www.imagemagick.org/script/download.php 
    # and specify the full path to magick.exe:
    
    magick_path = 'C:/Program Files/ImageMagick-7.0.10-Q16/magick.exe'
    
    # If you want to save the animation as .mp4, .avi or .mov, you 
    # should specify the full path to the ffmpeg executable in 
    # ImageMagick installation folder:
    
    ffmpeg_path = 'C:/Program Files/ImageMagick-7.0.10-Q16/ffmpeg.exe'
    
    def __init__(self, thickness, nmodes_sym, nmodes_antisym, fd_max, vp_max, 
                 c_L, c_S, c_R = None, fd_points=100, vp_step=100, 
                 material=''):
        """"
        Parameters
        ----------
        thickness : float or int
            Thickness of the plate, in mm.
        nmodes_sym : int
            Number of symmetric modes to calculate.
        nmodes_antisym : int
            Number of antisymmetric modes to calculate.
        fd_max : float or int
            Maximum value of frequency × thickness to calculate.
        vp_max : float or int
            Maximum value of phase velocity to calculate, in m/s.
        c_L : float or int
            Longitudinal wave velocity of the material, in m/s.
        c_S: float or int
            Shear wave velocity of the material, in m/s.
        c_R: float or int, optional
            Rayleigh wave velocity of the material, in m/s.     
        fd_points : int, optional
            Number of frequency × thickness points.
        vp_step : int, optional
            Increment between phase velocity intervals.
        material : str, optional
            Name of the material being analyzed.
            
        """
        
        self.d = thickness/1e3
        self.h = (thickness/2)/1e3
        self.nmodes_sym = nmodes_sym
        self.nmodes_antisym = nmodes_antisym
        self.fd_max = fd_max
        self.vp_max = vp_max
        self.c_L = c_L
        self.c_S = c_S
        self.c_R = c_R
        self.fd_points = fd_points
        self.vp_step = vp_step
        self.material = material
        
        # Solve the dispersion equations.
        
        sym = self._solve_disp_eqn(function=self._symmetric, 
                                   nmodes=nmodes_sym, 
                                   c=c_S, 
                                   label='S')
        
        antisym = self._solve_disp_eqn(function=self._antisymmetric, 
                                       nmodes=nmodes_antisym, 
                                       c=c_L, 
                                       label='A')
        
        # Calculate group velocity (vg) and wavenumber (k) from phase 
        # velocity (vp) and interpolate all results.
        
        self.vp_sym, self.vg_sym, self.k_sym = interpolate(sym, self.d)
        self.vp_antisym, self.vg_antisym, self.k_antisym = interpolate(antisym, 
                                                                       self.d)

    def _calc_constants(self, vp, fd):
        """Calculate the constants p and q (defined to simplify the 
        dispersion equations) and wavenumber from a pair of phase 
        velocity and frequency × thickness product.
        
        Parameters
        ----------
        vp : float or int
            Phase velocity.
        fd : float or int
            Frequency × thickness product.
        
        Returns
        -------
        k : float
            Wavenumber.
        p, q : float
            A pair of constants introduced to simplify the dispersion 
            relations.
            
        """
        
        omega = 2*np.pi*(fd/self.d)
        
        k = omega/vp
    
        p = np.sqrt((omega/self.c_L)**2 - k**2, dtype=np.complex128)
        q = np.sqrt((omega/self.c_S)**2 - k**2, dtype=np.complex128)

        return k, p, q

    def _symmetric(self, vp, fd):
        """Rayleigh-Lamb frequency relation for symmetric modes, used to 
        determine the velocity at which a wave of a particular frequency 
        will propagate within the plate. The roots of this equation are 
        used to generate the dispersion curves.

        Parameters
        ----------
        vp : float or int
            Phase velocity.
        fd : float or int
            Frequency × thickness product.
        
        Returns
        -------
        symmetric : float
            Dispersion relation for symmetric modes.
        
        """
        
        k, p, q = self._calc_constants(vp, fd)
    
        symmetric = (np.tan(q*self.h)/q
                     + (4*(k**2)*p*np.tan(p*self.h))/(q**2 - k**2)**2)

        return np.real(symmetric)
        
    def _antisymmetric(self, vp, fd):
        """Rayleigh-Lamb frequency relation for antisymmetric modes, 
        used to determine the velocity at which a wave of a particular 
        frequency will propagate within the plate. The roots of this 
        equation are used to generate the dispersion curves.

        Parameters
        ----------
        vp : float or int
            Phase velocity.
        fd : float or int
            Frequency × thickness product.
            
        Returns
        -------
        antisymmetric : float
            Dispersion relation for antisymmetric modes.
            
        """
        
        k, p, q = self._calc_constants(vp, fd)

        antisymmetric = (q * np.tan(q*self.h)
                         + (((q**2 - k**2)**2)*np.tan(p*self.h))/(4*(k**2)*p))

        return np.real(antisymmetric)
    
    def _calc_wave_structure(self, modes, vp, fd, y):
        """Calculate the wave structure across the thickness of the 
        plate.
        
        Parameters
        ----------
        modes : {'A', 'S'}
            Family of modes to analyze. Can be 'A' (antisymmetric modes) 
            or 'S' (symmetric modes).
        vp : float or int
            Phase velocity.
        fd : float or int
            Frequency × thickness product.
        y : array
            Array representing thickness values to calculate wave 
            structure, from -d/2 to d/2.
            
        Returns
        -------
        u : array
            In plane displacement profile.
        w : array
            Out of plane plane displacement profile.
            
        """
        
        k, p, q = self._calc_constants(vp, fd)
        
        if modes == 'S':
            C = 1
            B = -2*k*q*np.cos(q*self.h) / ((k**2 - q**2) * np.cos(p*self.h))   
            u = 1j*(k*B*np.cos(p*y) + q*C*np.cos(q*y))
            w = -p*B*np.sin(p*y) + k*C*np.sin(q*y)
        elif modes == 'A':
            D = 1
            A = 2*k*q*np.sin(q*self.h) / ((k**2 - q**2) * np.sin(p*self.h))
            u = 1j*(k*A*np.sin(p*y) - q*D*np.sin(q*y))
            w = p*A*np.cos(p*y) + k*D*np.cos(q*y)

        return u, w        

    def _solve_disp_eqn(self, function, nmodes, c, label):
        """Function to calculate the numerical solution to the 
        dispersion equations.
        
        The algorithm works as follows:
            
            1) Fix a value of frequency × thickness product.
            2) Evaluate the function at two values of phase velocity 
               (vp and vp+step) and check their signs.
            3) Since the function is continuous, if the sign changes  
               in the interval under analysis, a root exists in this 
               interval. Use the bisection method to locate it 
               precisely.
            4) Continue searching for other roots at this value of 
               frequency × thickness.
            5) Change the value of frequency × thickness and repeat 
               steps 2 to 4.

        Parameters
        ----------
        function : {self._symmetric, self._antisymmetric}
            Family of modes to solve. Can be `self._symmetric` (to solve 
            symmetric modes) or `self._antisymmetric` (to solve 
            antisymmetric modes).
            
        Returns
        -------
        result_dict : dict
            A dictionary, where the keys are the corresponding mode 
            (e.g., 'A0', 'A1', 'A2', ..., 'An' for antisymmetric modes 
             or 'S0', 'S1', 'S2', ..., 'Sn' for symmetric modes) and the 
            values are numpy arrays of dimensions (2, fd_points), where 
            the first row has the fd values and the second row has the 
            phase velocity values calculated.
            
        """
                        
        fd_arr = np.linspace(0, self.fd_max, self.fd_points)      
        result = np.zeros((len(fd_arr), nmodes + 1))
        
        print(f'\nCalculating {function.__name__[1:]} modes..\n')
        
        for i, fd in enumerate(fd_arr):
            
            print(f'{i}/{self.fd_points} - {np.around(fd, 1)} kHz × mm')
        
            result[i][0] = fd
        
            j = 1
            
            vp_1 = 0
            vp_2 = self.vp_step

            while vp_2 < self.vp_max:
                x_1 = function(vp_1, fd)
                x_2 = function(vp_2, fd)

                if j < nmodes + 1:
                    if not np.isnan(x_1) and not np.isnan(x_2):
                        if np.sign(x_1) != np.sign(x_2):
                            bisection = scipy.optimize.bisect(f=function, 
                                                              a=vp_1, 
                                                              b=vp_2, 
                                                              args=(fd,))
                            
                            # TO FIX: I don't know why at some points 
                            # the function changes sign, but the roots
                            # found by the bisect method don't evaluate 
                            # to zero.
                            
                            # For now, these values are ignored (only
                            # take into account those values that 
                            # evaluate to 0.01 or less).
                            
                            if (np.abs(function(bisection, fd)) < 1e-2 and not
                                    np.isclose(bisection, c)):
                                
                                result[i][j] = bisection
                                j += 1
                                
                vp_1 = vp_2
                vp_2 = vp_2 + self.vp_step
        
        # Correct some instabilities and replace zeros with NaN, so it
        # is easier to filter.
           
        result = correct_instability(result, function)
        result[result == 0] = np.nan
        
        result_dict = {}
        
        for nmode in range(nmodes):
            
            # Filter all NaN values.
            
            mode_result = np.vstack((result[:, 0], result[:, nmode + 1]))
            mode_result = mode_result[:, ~np.isnan(mode_result).any(axis=0)]
            
            # Append to a dictionary with keys 'An' or 'Sn'.
            
            result_dict[label + str(nmode)] = mode_result
        
        return result_dict
    
    def animate_displacement(self, mode, fd, speed=30, 
                             save_gif=False, save_video=False): 
        """Generate an animation of the displacement vector field across 
        the plate. The mesh grid created cover a full wavelength of the 
        current selected wave mode and fd value.
        
        Parameters
        ----------
        mode : str
            Mode to be animated. Can be "A0", "A1", "A2", ..., "An" or 
            "S0", "S1", "S2", ..., "Sn", with 'n' being the order of the 
            corresponding mode.
        fd : float or int
            Frequency × thickness product.
        speed : int
            Delay between frames in milliseconds. It can be used to 
            control the speed of the rotating vectors in the animation
            (a smaller value produces a faster animation). Default to 30.
        save_gif : bool
            Set to True if you want to save the result animation as a 
            gif. Defaults to False.
        save_video : {'mp4', 'mov', 'avi'}
            Choose a video format if you want to save the result 
            animation as a video. Can be 'mp4', 'mov' or 'avi'. 
            Defaults to False.
        
        Returns
        -------
        fig, ax : matplotlib objects
            The figure and the axes of the generated plot.
        
        """
        
        if mode[0] == 'S' and int(mode[1:]) < self.nmodes_sym:
            vp = self.vp_sym[mode](fd)
        elif mode[0] == 'A' and int(mode[1:]) < self.nmodes_antisym:
            vp = self.vp_antisym[mode](fd)
        else:
            raise Exception('mode not recognized. Mode must be "Sn" or '
                            '"An", where n is an integer greater or equal '
                            'than 0. For example: "S0", "S1", "A0", "A1", '
                            'etc. Make sure the mode order selected is within '
                            'the number of modes requested when setting up the'
                            ' Lamb class.')
            
        # Generate the mesh grid, with the x-values covering a full
        # wavelength and the y-values covering the thickness of the 
        # plate (from -thickness/2 to +thickness/2).
        
        wavelength = vp/(fd/self.d)
        
        xx = np.linspace(0, wavelength, 40)
        yy = np.linspace(-self.h, self.h, 40)
        
        x, y = np.meshgrid(xx, yy)
        u, w = np.zeros_like(x), np.zeros_like(y)
        
        # Generate the time vector necessary to complete one cycle 
        # (i.e., wave period).
        
        time = np.linspace(0, 1/(fd/self.d), 30)
        
        # Calculate angular frequency and wavenumber.
        
        omega = 2*np.pi*(fd/self.d)  
        k = omega/vp
                 
        def compute_displacement(t):
            """Calculate particle displacement as a function of time."""
            
            u, w = self._calc_wave_structure(mode[0], vp, fd, y)

            u = u * np.exp(1j*(k*x-omega*t))
            w = w * np.exp(1j*(k*x-omega*t))

            return np.real(u), np.real(w)

        # Find the largest displacement vector to use for normalization.

        max_disp_arr = []
        
        for t in time:
            u, w = compute_displacement(t)
            max_disp_arr.append(np.amax(np.sqrt(u**2 + w**2)))

        max_disp = max(max_disp_arr)
        
        # Generate the quiver plot animation.

        fig, ax = plt.subplots(figsize=(8, 5))
        fig.canvas.set_window_title(f'Displacement Field (mode {mode})')
        
        quiver = ax.quiver(x, y, u, w, scale=5*max_disp, scale_units='inches')
        
        ax.set_title('Mode $\mathregular{' + mode[0] + '_' + mode[1:] + '}$')
        ax.text(0.5, 0.05, f'fd = {np.around(fd, 1)} kHz × mm', ha='center', 
                va='center', transform = ax.transAxes)

        ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        
        ax.set_yticks([-self.h, 0, self.h])
        ax.set_yticklabels(['-d/2', '0', 'd/2'])
        
        ax.set_ylabel('Thickness')
        
        ax.set_xlim([0 - wavelength/4, wavelength + wavelength/4])
        ax.set_ylim([-self.d, self.d])
    
        def init():
            return quiver,
    
        def animate(t):
            u, w = compute_displacement(t)
            quiver.set_UVC(u, w)
            
            return quiver,
        
        anim = matplotlib.animation.FuncAnimation(fig, animate, init_func=init, 
                                                  frames=time, interval=speed, 
                                                  blit=True)
        
        if save_gif:
            plt.rcParams['animation.convert_path'] = Lamb.magick_path
            anim.save(f'results/Mode_{mode}_fd_{int(fd)}_animation.gif', 
                      writer='imagemagick', extra_args='convert')
            
        if save_video:
            plt.rcParams['animation.ffmpeg_path'] = Lamb.ffmpeg_path
            anim.save(f'results/Mode_{mode}_fd_{int(fd)}_animation.'
                      f'{save_video}', writer='imagemagick')
        
        return fig, ax
    
    def plot(self, ax, result, y_max, cutoff_frequencies=False, 
             arrow_dir=None, material_velocities=False, plt_kwargs={}):
        """Generate a dispersion plot for a family of modes (symmetric
        or antisymmetric).
        
        Parameters
        ----------
        ax : axes
            Matplotlib axes in which the plot will be added.
        result : dict
            A dictionary with a result (vp, vg or k) interpolator at 
            each mode.
        y_max : float or int
            Maximum y value in the plot. 
        cutoff_frequencies : bool, optional
            Set to True to add cutoff frequencies to the plot.
        arrow_dir : {'up', 'down'}, optional
            Set arrows direction of cutoff frequencies. Can be 'up' (for 
            group velocity plots) or 'down' (for phase velocity plots).
        material_velocities : bool, optional
            Add material velocities (longitudinal, shear and Rayleigh) 
            to the plot. Defaults to True.       
        plt_kwargs : dict, optional
            Matplotlib kwargs (to change color, linewidth, linestyle, 
            etc.).
                    
        """
        
        for mode, arr in result.items():
            
            # Generate an fd array for each mode and add the 
            # corresponding mode plot.
            
            fd = np.arange(np.amin(arr.x), np.amax(arr.x), 0.1)
            add_plot(ax, result, mode, fd, **plt_kwargs)   
            
            if cutoff_frequencies:
                add_cutoff_freqs(ax, mode, arrow_dir, y_max, 
                                 self.c_L, self.c_S)
        
        if material_velocities:
            add_velocities(ax, self.c_L, self.c_S, self.c_R, self.fd_max)
                
        ax.set_xlim([0, self.fd_max])
        ax.set_ylim([0, y_max])
        
        ax.set_xlabel('Frequency × thickness [KHz × mm]')
        
    def plot_phase_velocity(self, modes='both', cutoff_frequencies=True, 
                            material_velocities=True, save_img=False,
                            sym_style={'color': 'blue'}, 
                            antisym_style={'color': 'red'}):
        """Generate a plot of phase velocity as a function of frequency 
        × thickness.
        
        Parameters
        ----------
        modes : {'both', 'symmetric', 'antisymmetric'}, optional
            Which family of modes to plot. Can be 'symmetric', 
            'antisymmetric' or 'both'. Defaults to 'both'.
        cutoff_frequencies : bool, optional
            Add cutoff frequencies to the plot. Defaults to True.
        material_velocities : bool, optional
            Add material velocities (longitudinal, shear and Rayleigh) 
            to the plot. Defaults to True.
        save_img : bool, optional
            Save the result image as png. Defaults to False.
        sym_style : dict, optional
            A dictionary with matplotlib kwargs to modify the symmetric 
            curves (to change color, linewidth, linestyle, etc.).
        antisym_style : dict, optional
            A dictionary with matplotlib kwargs to modify the 
            antisymmetric curves (to change color, linewidth, linestyle, 
            etc.).
            
        Returns
        -------
        fig, ax : matplotlib objects
            The figure and the axes of the generated plot.
            
        """
        
        fig, ax = plt.subplots(figsize=(7, 4))
        fig.canvas.set_window_title('Phase Velocity')
        
        # Calculate the maximum value to scale the ylim of the axes.
        
        max_sym, max_antisym = find_max(self.vp_sym), find_max(self.vp_antisym)
        
        if modes == 'symmetric':
            self.plot(ax, self.vp_sym, max_sym, cutoff_frequencies, 'down', 
                      material_velocities, plt_kwargs=sym_style)
        elif modes == 'antisymmetric':
            self.plot(ax, self.vp_antisym, max_antisym, cutoff_frequencies, 
                      'down', material_velocities, plt_kwargs=antisym_style)
        elif modes == 'both':
            max_ = max(max_sym, max_antisym)
            self.plot(ax, self.vp_sym, max_, cutoff_frequencies, 
                      'down', material_velocities, plt_kwargs=sym_style)
            self.plot(ax, self.vp_antisym, max_, cutoff_frequencies, 
                      'down', material_velocities, plt_kwargs=antisym_style)
        else:
            raise Exception('modes must be "symmetric", "antisymmetric"'
                            'or "both".') 
            
        ax.legend(loc='lower right')
        ax.set_ylabel('Phase Velocity [m/s]')
        
        if save_img:
            fig.savefig(f'results/Phase Velocity - {self.d*1e3} mm '
                        f'{self.material} plate.png', 
                        bbox_inches='tight')
        
        return fig, ax

    def plot_group_velocity(self, modes='both', cutoff_frequencies=True, 
                            save_img=False, sym_style={'color': 'blue'}, 
                            antisym_style={'color': 'red'}):
        """Generate a plot of group velocity as a function of frequency 
        × thickness.
        
        Parameters
        ----------
        modes : {'both', 'symmetric', 'antisymmetric'}, optional
            Which family of modes to plot. Can be 'symmetric', 
            'antisymmetric' or 'both'. Defaults to 'both'.
        cutoff_frequencies : bool, optional
            Add cutoff frequencies to the plot. Defaults to True.
        save_img : bool, optional
            Save the result image as png. Defaults to False.            
        sym_style : dict, optional
            A dictionary with matplotlib kwargs to modify the symmetric 
            curves (to change color, linewidth, linestyle, etc.).
        antisym_style : dict, optional
            A dictionary with matplotlib kwargs to modify the 
            antisymmetric curves (to change color, linewidth, linestyle, 
            etc.).
            
        Returns
        -------
        fig, ax : matplotlib objects
            The figure and the axes of the generated plot.
            
        """
              
        fig, ax = plt.subplots(figsize=(7, 4))
        fig.canvas.set_window_title('Group Velocity')
        
        # Calculate the maximum value to scale the ylim of the axes.
        
        max_sym, max_antisym = find_max(self.vg_sym), find_max(self.vg_antisym)
    
        if modes == 'symmetric':
            self.plot(ax, self.vg_sym, max_sym, cutoff_frequencies, 
                      'up', plt_kwargs=sym_style)
        elif modes == 'antisymmetric':
            self.plot(ax, self.vg_antisym, max_antisym, cutoff_frequencies, 
                      'up', plt_kwargs=antisym_style)
        elif modes == 'both':
            max_ = max(max_sym, max_antisym)
            self.plot(ax, self.vg_sym, max_, cutoff_frequencies, 
                      'up', plt_kwargs=sym_style)
            self.plot(ax, self.vg_antisym, max_, cutoff_frequencies, 
                      'up', plt_kwargs=antisym_style)
        else:
            raise Exception('modes must be "symmetric", "antisymmetric"'
                            'or "both".')  
             
        ax.legend(loc='lower right')
        ax.set_ylabel('Group Velocity [m/s]')
        
        if save_img:
            fig.savefig(f'results/Group Velocity - {self.d*1e3} mm '
                        f'{self.material} plate.png', 
                        bbox_inches='tight')    
            
        return fig, ax
    
    def plot_wave_number(self, modes='both', save_img=False,
                         sym_style={'color': 'blue'}, 
                         antisym_style={'color': 'red'}):       
        """Generate a plot of wavenumber as a function of frequency × 
        thickness.
        
        Parameters
        ----------
        modes : {'both', 'symmetric', 'antisymmetric'}, optional
            Which family of modes to plot. Can be 'symmetric', 
            'antisymmetric' or 'both'. Defaults to 'both'.
        save_img : bool, optional
            Save the result image as png. Defaults to False.                  
        sym_style : dict, optional
            A dictionary with matplotlib kwargs to modify the symmetric 
            curves (to change color, linewidth, linestyle, etc.).
        antisym_style : dict, optional
            A dictionary with matplotlib kwargs to modify the 
            antisymmetric curves (to change color, linewidth, linestyle, 
            etc.).
            
        Returns
        -------
        fig, ax : matplotlib objects
            The figure and the axes of the generated plot.
            
        """
            
        fig, ax = plt.subplots(figsize=(7, 4))
        fig.canvas.set_window_title('Wave Number')
        
        # Calculate the maximum value to scale the ylim of the axes.
        
        max_sym, max_antisym = find_max(self.k_sym), find_max(self.k_antisym)

        if modes == 'symmetric':
            self.plot(ax, self.k_sym, max_sym, plt_kwargs=sym_style)
        elif modes == 'antisymmetric':
            self.plot(ax, self.k_antisym, max_antisym, plt_kwargs=antisym_style)
        elif modes == 'both':
            max_ = max(max_sym, max_antisym)
            self.plot(ax, self.k_sym, max_, plt_kwargs=sym_style)
            self.plot(ax, self.k_antisym, max_, plt_kwargs=antisym_style)
        else:
            raise Exception('modes must be "symmetric", "antisymmetric"'
                            'or "both".') 
            
        ax.legend(loc='upper left')
        ax.set_ylabel('Wave Number [1/m]')
        
        if save_img:
            fig.savefig(f'results/Wave Number - {self.d*1e3} mm '
                        f'{self.material} plate.png',
                        bbox_inches='tight')
            
        return fig, ax
    
    def plot_wave_structure(self, mode, nrows, ncols, fd, save_img=False,
                            inplane_style={'color': 'blue'}, 
                            outofplane_style={'color': 'red'}):
        """Generate a plot of the wave structure, i.e., the in-plane and 
        out-of-plane displacement profiles across the thickness of the 
        plate.
        
        Parameters
        ----------
        mode : str
            Mode to be analyzed. Can be "A0", "A1", "A2", ..., "An" or 
            "S0", "S1", "S2", ..., "Sn", with 'n' being the order of the 
            corresponding mode.
        nrows : int
            Number of rows in the subplot.
        ncols : int
            Number of columns in the subplot.
        fd : array
            Array with the frequency × thickness points to analyze. The 
            length of the array must be equal to nrows x ncols.
        save_img : bool, optional
            Save the result image as png. Defaults to False.                  
        inplane_style : dict, optional
            A dictionary with matplotlib kwargs to modify the in-plane 
            curves (to change color, linewidth, linestyle, etc.).
        outofplane_style : dict, optional
            A dictionary with matplotlib kwargs to modify the 
            out-of-plane curves (to change color, linewidth, linestyle, 
            etc.).
            
        Returns
        -------
        fig, axs : matplotlib objects
            The figure and the axes of the generated plot.
        
        """
        
        y = np.linspace(-self.h, self.h, 100) 
           
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols)        
        fig.canvas.set_window_title(f'Wave Structure (mode {mode})')
        
        fig.suptitle('Mode $\mathregular{' + mode[0] + '_' + mode[1:] + '}$')
        
        for ax, freq in zip(axs.flatten(), fd):
            if mode[0] == 'S' and int(mode[1:]) < self.nmodes_sym:
                vp = self.vp_sym[mode](freq)
            elif mode[0] == 'A' and int(mode[1:]) < self.nmodes_antisym:
                vp = self.vp_antisym[mode](freq)
            else:
                raise Exception('mode not recognized. Mode must be "Sn" or '
                                '"An", where n is an integer greater or equal '
                                'than 0. For example: "S0", "S1", "A0", "A1", '
                                'etc. Make sure the mode order selected is '
                                'within the number of modes requested when '
                                'setting up the Lamb class.')
            
            u, w = self._calc_wave_structure(mode[0], vp, freq, y)
            
            # All values of u, w are purely real or purely imaginary.
            
            if np.all(np.iscomplex(u)):
                ax.plot(np.imag(u), y, label='In plane', **inplane_style)
            else:
                ax.plot(np.real(u), y, label='In plane', **inplane_style)
                
            if np.all(np.isreal(w)):
                ax.plot(np.real(w), y, label='Out of plane', **outofplane_style)
            else:
                ax.plot(np.imag(w), y, label='Out of plane', **outofplane_style)                
            
            ax.set_title(f'fd: {np.around(freq, 1)} KHz × mm')
    
            ax.set_ylim([-self.h, self.h])
            
            ax.set_yticks([-self.h, 0, self.h]) 
            ax.set_yticklabels(['-d/2', '0', 'd/2'])
            
            ax.spines['left'].set_position('zero')
            ax.spines['right'].set_color('none')
            ax.spines['bottom'].set_position('zero')
            ax.spines['top'].set_color('none')
            
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
        
        # TO FIX: For some reason tight_layout() isn't working with some
        # subplot configurations, producing overlaping plots (e.g.
        # nrows=2 and ncols=4). This happens even if I remove the
        # fig.suptitle() and fig.legend() (not considered by
        # tight_layout())
        
        fig.tight_layout()
        
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', ncol=2)  

        if save_img:
            fig.savefig(f'results/Wave Structure - {self.d*1e3} mm '
                        f'{self.material} plate - Mode {mode}.png',
                        bbox_inches='tight')
            
        return fig, axs
    
    def save_results(self):
        """Save all results to a txt file."""
        
        if self.material:
            filename = f'{self.material} plate - {self.d*1e3} mm.txt'
        else:
            filename = f'{self.d*1e3} mm plate.txt'
            
        header = (f'Material: {self.material}\n'
                  f'Thickness: {str(self.d*1e3)} mm\n'
                  f'Longitudinal wave velocity: {str(self.c_L)} m/s\n'
                  f'Shear wave velocity: {str(self.c_S)} m/s\n\n')
                
        write_txt(self.vp_sym, self.vp_antisym, 'Phase Velocity', 
                  filename, header)
        write_txt(self.vg_sym, self.vg_antisym, 'Group Velocity', 
                  filename, header)
        write_txt(self.k_sym, self.k_antisym, 'Wavenumber', 
                  filename, header)