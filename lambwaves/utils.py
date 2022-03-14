"""Some helper functions used to calculate dispersion curves.

Author:         Francisco Rotea
                (Buenos Aires, Argentina)
Repository:     https://github.com/franciscorotea
Email:          francisco.rotea@gmail.com

"""

import itertools

import numpy as np
import scipy.interpolate

def interpolate(result, d, kind='cubic'):
    """Interpolate the results for phase velocity, group velocity and
    wave number.
    
    Parameters
    ----------
    result : dict
        Dictionary with the phase velocity values obtained by solving
        the dispersion equations.
    kind : str, optional
        Specifies the kind of interpolation as a string. Can be
        ‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, 
        ‘previous’, ‘next’. Defaults to 'cubic'.
        
    Returns
    -------
    interp_vp : dict
        Dictionary with phase velocity interpolator at each mode.
    interp_vg : dict
        Dictionary with group velocity interpolator at each mode.
    interp_k : dict
        Dictionary with wave number interpolator at each mode.
        
    """
    
    interp_vp = {}
    interp_vg = {}
    interp_k = {}
    
    for mode, arr in result.items():
            
            if arr[1].size > 3:
                
                fd = arr[0]
                vp = arr[1]             

                interp_vp[mode] = scipy.interpolate.interp1d(fd, vp, kind=kind)
                
                k = (fd*2*np.pi/d)/vp
                
                interp_k[mode] = scipy.interpolate.interp1d(fd, k, kind=kind)
                
                # Find the derivative of phase velocity using a 
                # interpolating spline.
                
                univ_s = scipy.interpolate.InterpolatedUnivariateSpline(fd, vp)
                vp_prime = univ_s.derivative()
                
                vg = np.square(vp) * (1/(vp - vp_prime(fd)*fd))

                interp_vg[mode] = scipy.interpolate.interp1d(fd, vg, kind=kind)
    
    return interp_vp, interp_vg, interp_k

def correct_instability(result, function):
    """A function to correct the instability produced when two roots are 
    in close proximity, making the function change sign twice or more in 
    the phase velocity interval under analysis. Since these values of 
    phase velocity are not computed, it ultimately produces a wrong mode 
    assignation, e.g., a phase velocity value corresponding to the S1 
    mode is wrongly assigned to S0.
                    
    Since phase velocity is strictly decreasing for each mode order 
    (except for A0), the function works by looping through each mode, 
    detecting if a phase velocity value is not decreasing. If this is 
    the case, the value is appended to the next mode, and replaced by 0.
    
    Parameters
    ----------
    result : array
        An array of shape (fd_points, nmodes+1), where the first column
        contains the fd values and the following columns are the phase 
        velocity values for the requested modes (S0, S1, S2, etc., or 
        A0, A1, A2, etc.)
    function : object
        Family of modes to solve (symmetric or antisymmetric).
        
    Returns
    -------
    corr : array
        The corrected result array.
        
    """
    
    # Compensate for antisymmetric mode (A0 is strictly increasing).
    
    n = 1 if function.__name__ == '_symmetric' else 2
    nmodes = result.shape[1] - 1

    corr = np.copy(result)
    
    for idx, col in enumerate(corr.T[n:,:]):
        if np.any(col):
            i = 0
            while col[i] == 0 and i < len(col)-1:
                i += 1
            if idx < nmodes-n:
                corr[i][idx+n+1] = 0
            
    for idx, col in enumerate(corr.T[n:,:]):
        for i in range(len(col)-1):
            if i == len(col)-2:
                corr[i+1][idx+n] = 0
            if col[i] != 0:
                j = i + 1
                if col[j] == 0:
                    while col[j] == 0 and j < len(col)-1:
                        j += 1
                if j < len(col)-1:
                    if col[i] < col[j] or col[j] == 0:
                        while (col[i] < col[j] or col[j] == 0) and j < len(col)-1:
                            if col[j] == 0:
                                    j += 1
                            else:
                                for idx2 in range(nmodes):
                                    if idx == idx2:
                                        corr[j][idx+n] = 0
                                        p = n + 1
                                        while p <= nmodes - idx2:
                                            corr[j][idx+p] = result[j][idx+p-1]
                                            p += 1
                                j += 1
    
    return corr

def write_txt(data_sym, data_antisym, kind, filename, header):
    """Function to write the results to a txt file.
    
    Parameters
    ----------
    data_sym : dict
        A dictionary consisting of interpolators for the specified 
        symmetric modes.
    data_antisym : dict
        A dictionary consisting of interpolators for the specified 
        antisymmetric modes.
    kind : {'Phase Velocity', 'Group Velocity', 'Wavenumber'}
        The type of results to write. Can be 'Phase Velocity', 'Group
        Velocity' or 'Wavenumber'.
    filename : str
        The filename of the txt file.
    header : str
        The header of the txt file (to include material information for 
        example)

    """
    
    if kind == 'Phase Velocity':
        label = 'vp [m/s]'
    elif kind == 'Group Velocity':
        label = 'vg [m/s]'
    else:
        label = 'k  [1/m]'
        
    # Get the calculated (non-interpolated) data.
    
    results = []
    
    for n in range(len(data_sym)):
        x_vals = data_sym[f'S{n}'].x
        y_vals = data_sym[f'S{n}'](x_vals)
        results.append(np.around(x_vals, 1))
        results.append(np.around(y_vals, 1))
        
    for n in range(len(data_antisym)):
        x_vals = data_antisym[f'A{n}'].x
        y_vals = data_antisym[f'A{n}'](x_vals)
        results.append(np.around(x_vals, 1))
        results.append(np.around(y_vals, 1))
    
    # Write the results in a txt file.
    
    with open('results/' + kind + ' - ' + filename, 'w') as f:
        f.write(header)
        
        f.write('\t\t\t\t'.join(data_sym.keys()) + '\t\t\t\t')
        f.write('\t\t\t\t'.join(data_antisym.keys()) + '\n')
        
        for _ in range(len(data_sym) + len(data_antisym)):
            f.write('fd [kHz mm]\t' + label + '\t')
            
        f.write('\n')
        
        for k in itertools.zip_longest(*results, fillvalue=''):
             f.write('\t\t'.join(map(str, k)) + '\n')
        
def find_max(result):
    """Find the maximum value in all modes analyzed. Used to limit the 
    scale of the dispersion plots.
    
    Parameters
    ---------- 
    result : dict
        A dictionary with a result (vp, vg or k) interpolator at each 
        mode.
        
    """
    
    max_val_arr = []
    
    for mode, arr in result.items():
        fd = np.arange(np.amin(arr.x), np.amax(arr.x), 0.1)
        max_val_arr.append(np.amax(result[mode](fd)))  
    
    return max(max_val_arr)