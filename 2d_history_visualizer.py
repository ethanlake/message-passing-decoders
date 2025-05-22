import numpy as np
import matplotlib.pyplot as plt
import h5py 
import argparse
from matplotlib.widgets import Slider
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation

parser = argparse.ArgumentParser() 
parser.add_argument('-fin',default='no') # spacetime history of fields and spins 
parser.add_argument('-save_animation',default='no') # saves animation for max_animation_time (defined below) if not equal to "no" 
args = parser.parse_args() 

def loaddata(fin):
    with h5py.File(fin, "r") as f:
        key_list = [key for key in f.keys() if not key.startswith("_")] 
        data_dict = dict.fromkeys(key_list)
        for key in key_list:
            data_dict[key] = f[key][()]  
    return data_dict 

# load field history data 
hist_data = loaddata(args.fin)
L = hist_data["L"]
hist = hist_data["hist"].T 
anyon_hist = hist_data["anyon_hist"].T
field_hist = hist_data["field_hist"].T
field_hist[field_hist == 0] = 100000
maxtime = np.shape(hist)[0]

if "power" in args.fin:
    force_hist = hist_data["force_hist"].T 
    maxforce = np.max(force_hist[:,:,:,1])

sz = 6.5
fig, ax = plt.subplots(figsize=(sz,sz))
ax.axis('off')  

if args.save_animation == "no":
    plt.subplots_adjust(bottom=0.25) 
    slider_ax = plt.axes([0.2, 0.1, 0.65, 0.03])  
else: 
    slider_ax = plt.axes([-.5, -.5, 0, 0])  

plot_slider = Slider(slider_ax,r'',0,maxtime-1,valinit=0, valstep=1)

def update(val):
    ax.clear()
    ax.axis('off') 
    time = int(plot_slider.val)+1

    if "power" in args.fin: 
        mags = force_hist[time-1,:,:,0] / maxforce 
        ax.imshow(mags.T, origin='lower',cmap='gist_heat_r',vmax=1,vmin=0)
    else:
        minfields = (1/np.min(field_hist[time-1,:,:,:,:],axis=(2,3)))**.75 # to make the color gradients look nice 
        ax.imshow(minfields.T, origin='lower',cmap='Purples',vmin=0,vmax=1)

    fig.canvas.draw_idle()


if args.save_animation != 'no':
    fps = 20
    dpi = 150 

    def animate(frame):
        plot_slider.set_val(frame)  
        update(frame)
    
    max_animation_time = maxtime
    animation = FuncAnimation(fig, animate, frames=np.arange(0, max_animation_time))
    
    animation.save(args.save_animation, writer='ffmpeg', dpi=dpi,fps=fps)

else: 
    def on_key(event):
        if event.key == 'left':
            plot_slider.set_val(max(plot_slider.val - 1, plot_slider.valmin))
        elif event.key == 'right':
            plot_slider.set_val(min(plot_slider.val + 1, plot_slider.valmax))
    fig.canvas.mpl_connect('key_press_event', on_key)

    # connect the slider to the update function
    plot_slider.on_changed(update)

    plt.show()

