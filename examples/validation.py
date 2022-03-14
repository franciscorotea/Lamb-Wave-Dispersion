import matplotlib.pyplot as plt
import pandas as pd
from lambwaves import Lamb
    
# Load phase and group velocity data from an Excel file exported from 
# Dispersion software (1 mm aluminum plate)

df_vp = pd.read_excel('alum1mm.xlsx', skiprows = 9, usecols = 'U:AN', 
                      header=None)
df_vg = pd.read_excel('alum1mm.xlsx', skiprows = 9, usecols = 'A:T', 
                      header=None)

# Create an instance of the same material using the Lamb class.

alum = Lamb(thickness=1, 
            nmodes_sym=5, 
            nmodes_antisym=5, 
            fd_max=10000, 
            vp_max=15000, 
            c_L=6420,
            c_S=3040)

# Plot phase velocity using the Lamb class.

fig1, ax1 = alum.plot_phase_velocity(material_velocities=False,
                                     cutoff_frequencies=False,
                                     sym_style={'color' : 'black'}, 
                                     antisym_style={'color' : 'black'})

# Remove the legend that labels Symmetric and Antisymmetric modes
# (we are interested in labeling only Lamb module and Dispersion).

ax1.get_legend().remove()  
 
# Plot phase velocity obtained by Dispersion.

line1 = ax1.lines[0]

for mode in df_vp.columns[::2]:
    ax1.plot(df_vp[mode]*1e3, df_vp[mode+1]*1e3, 
             color = 'orange', 
             linestyle='--')
    
line2 = ax1.lines[-1]    

ax1.legend((line1, line2), ('Lamb module', 'Dispersion'))

# Plot group velocity using the Lamb class.

fig2, ax2 = alum.plot_group_velocity(cutoff_frequencies=False,
                                     sym_style={'color' : 'black'},
                                     antisym_style={'color' : 'black'})

# Remove the legend that labels Symmetric and Antisymmetric modes
# (we are interested in labeling only Lamb module and Dispersion).

ax2.get_legend().remove()

# Plot group velocity obtained by Dispersion.
   
line1 = ax2.lines[0]

for mode in df_vg.columns[::2]:
    ax2.plot(df_vg[mode]*1e3, df_vg[mode+1]*1e3, 
             color = 'orange', 
             linestyle='--')

line2 = ax2.lines[-1]       

ax2.legend((line1, line2), ('Lamb module', 'Dispersion'))

plt.show()
