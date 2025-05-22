import numpy as np
import matplotlib.pyplot as plt
import h5py 
import argparse
import matplotlib.cm as cm 

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman'] + plt.rcParams['font.serif']
pi = np.pi 

# load data from h5 files 
def loaddata(fin): 
    with h5py.File(fin, "r") as f:
        key_list = [key for key in f.keys() if not key.startswith("_")] # saving dictionaries to jld2 files can sometimes cause the creation of additional keys describing the data structure that we don't care about
        data_dict = dict.fromkeys(key_list)
        for key in key_list:
            data_dict[key] = f[key][()]  
    return data_dict 

parser = argparse.ArgumentParser() 
parser.add_argument('-fin',nargs='+')  
parser.add_argument('-save',default='') 
parser.add_argument('-full_stats',action='store_true') # if true, plots variance etc. of erode times
parser.add_argument('-just_te',action='store_true') # if true, only plots erosion times
parser.add_argument('-just_f',action='store_true') # if true, only plots erosion fidelity
parser.add_argument('-error_rate',action='store_true') # if true, plots logical error rate
parser.add_argument('-extend_fit',action='store_true') # if true, extends the fit to the trel data to small p 
parser.add_argument('-max_te',action='store_true') # if true, plots the max erosion time 
args = parser.parse_args() 

### aesthetics ### 
ms = 7 # marker size 
plotsize = 3.5 # size of plots (in.)
lw = 1.2 # line width 
legend_title_fontsize = 13.5; legend_fontsize = 12 

### load data ### 
fins = args.fin   
data = np.zeros(len(fins),dtype=object)
for (f,fin) in enumerate(fins): 
    data[f] = loaddata(fin)

synch = data[0]["synch"]
mode = data[0]["mode"].decode('utf-8') 
vary_L = data[0]["vary_L"]
if vary_L: 
    cmap = cm.Spectral 
else: 
    cmap = cm.coolwarm 

### plot labels etc. ### 
title = r"${\sf TC}$" if "2d" in args.fin[0] else r"${\sf rep.\, code}$"
# title = r"${\sf TC\, \, |\,\,{\sf uncoord}}$" if "2d" in args.fin[0] else r"${\sf rep.\, code}$"
xlab = r"$L$" if vary_L else r"$p$"
if mode == "erode": 
    ylab = r"$T_{\sf dec}$" if args.just_te else r"$p_{\sf log}$"
    # ylab = r"$T_{\sf dec}$" if args.just_te else r"$p_{\sf log}^{\sf un}$"
elif mode == "Ft":
    ylab = r"$F_T$" if not args.error_rate else r"$\varepsilon$"
elif mode == "trel":
    ylab = r"$t_{\sf mem}$"    
leg_title = r'$L$' if not vary_L else r'$p$'
leg_loc = 'best'

nplots = 2 if mode == "hist" else 1
if mode == "erode": # various options for offline decoding plots 
    just_te = args.just_te and vary_L
    just_f = args.just_f 
    max_te = args.max_te 
    nplots = 1 if just_te or just_f or max_te else 2

figsize = (plotsize*nplots,plotsize)
fig,ax = plt.subplots(1,nplots,figsize=figsize)
plt.subplots_adjust(wspace=.4)
if nplots > 1: 
    for i in range(nplots):
        ax[i].minorticks_off()
else:
    ax.minorticks_off()

# plot history of spin configuration 
if "hist" in mode: 
    xlab = r'$x$'; ylab = r'$t$'
    hist = 1-data[0]["hist"].T
    if "powerdw" in fins[0]: 
        charge_hist = data[0]["charge_hist"].T
        field_hist = data[0]["field_hist"].T
        ax[0].imshow(hist,origin='lower',cmap='binary')
        ax[1].imshow(charge_hist,origin='lower',cmap='coolwarm',vmin=-1,vmax=1)
        img = ax[2].imshow(data[0]["field_hist"].T,origin='lower',cmap='rainbow')
        cbar = fig.colorbar(img,fraction=0.046, pad=0.04, ax = ax[2]) 
        cbar.minorticks_on()
        cbar.set_label(r'$F$')

    if "1dloc" in fins[0]: 
        ax[0].axis('off'); ax[1].axis('off')
        ax[0].imshow(hist,origin='lower',cmap='binary')
        force_hist = data[0]["field_hist"].T[:,:,0] - data[0]["field_hist"].T[:,:,1]
        zerolocs = np.where(force_hist == 0)
        power = .5
        inverted_forces = 1/(np.abs(force_hist + 1e-9))**power 
        inverted_forces[zerolocs] = 0
        img = ax[1].imshow(inverted_forces,origin='lower',cmap='coolwarm')
        cbar = fig.colorbar(img,fraction=0.046, pad=0.04, ax = ax[1])
        ax[0].set(xlabel=r'${\sf x}$',ylabel=r'${\sf t}$')
        ax[1].set(xlabel=r'${\sf x}$',ylabel=r'${\sf t}$')

elif "erode" in mode: 
    if args.max_te: # plot the max erosion time 
        for f in range(len(fins)): 
            els = data[f]["Ls"]
            p = data[f]["p"]
            erode_stats = data[f]["erode_stats"].T 
            erode_maxs = erode_stats[:,-1]
            ax.scatter(els,erode_maxs,label=r'$%.2f$'%p,c=cmap(f/len(fins)),ec='k',lw=lw,marker='o',s=70,alpha=1)
        ylab=r'$T_{\sf dec,\, max}$'

    elif (just_te and vary_L): # just plot the erosion times vs L 
        for f in range(len(fins)): 
            els = data[f]["Ls"]
            p = data[f]["p"]
            erode_times = data[f]["erode_times"].T  
            ax.scatter(els,erode_times,label=r'$%.2f$'%p,c=cmap(f/len(fins)),ec='k',lw=lw,marker='o',s=70,alpha=1)
            
            powerfit = not True # fit erode_times to power law?
            if powerfit: # T_e as a power in L 
                mp = int(len(els)/2)
                fits = np.polyfit(np.log(els)[:mp],np.log(erode_times)[:mp],1)
                ax.plot(els,np.exp(fits[1]) * els**(fits[0]),c=cmap(f/len(fins)),ls='-',lw=3*lw,alpha=1,label=r'$%.2f$'%(-fits[0]))
                ax.set_yscale("log")
            else: # fit T_e to log(L) 
                mp = -1
                fits = np.polyfit(np.log(els[:mp]),erode_times[:mp],1)
                ax.plot(els,np.log(els) * fits[0] + fits[1],c=cmap(f/len(fins)),ls='-',lw=3*lw,alpha=.6,zorder=-100)
        ax.set_xscale("log")

    elif just_f:
        for f in range(len(fins)): 
            els = data[f]["Ls"]
            ps = data[f]["ps"]
            p = data[f]["p"]
            L = data[f]["L"]
            erode_times = data[f]["erode_times"].T  
            fs = data[f]["erode_frac"].T 
            if vary_L:
                ax.scatter(els,1-fs,label=r'$%.3f$'%p,c=cmap(f/len(fins)),ec='k',lw=lw,marker='o',s=70,alpha=1)
                mp = -1 
                fits = np.polyfit(els[:mp],np.log((1-fs)[:mp]),1)
                pc = .055
                pc = .5 
                print(p,-fits[0]/np.log(pc/p))
                ax.plot(els,np.exp(els * fits[0] + fits[1]),c=cmap(f/len(fins)),ls='-',lw=3*lw,alpha=.6,zorder=-100)
                ax.set_yscale("log")
            else:
                ax.plot(ps,1-fs,label=r'$%d$'%L,c=cmap(f/len(fins)),mec='k',lw=lw,marker='o',ms=ms)

        if not vary_L:
            pass 
    else: 
        if vary_L:
            for f in range(len(fins)): 
                els = data[f]["Ls"]
                p = data[f]["p"]
                erode_fracs = data[f]["erode_frac"].T
                erode_times = data[f]["erode_times"].T 
                infs= 1-erode_fracs + 1e-9
                ax[0].plot(els,infs,label=r'$%.2f$'%p,c=cmap(f/len(fins)),mec='k',lw=lw,marker='o',ms=ms)
                ax[1].plot(els,erode_times,label=r'$%.2f$'%p,c=cmap(f/len(fins)),mec='k',lw=lw,marker='o',ms=ms)
                # fit erode_times to power law: 
                mp = int(len(els)/2)
                fits = np.polyfit(np.log(els)[:mp],np.log(erode_times)[:mp],1)
                # power = .5
                ax[1].plot(els,np.exp(fits[1]) * els**(fits[0]),c=cmap(f/len(fins)),ls='--',lw=3*lw,alpha=.7,label=r'$%.2f$'%(-fits[0]))
            
            ax[0].set(xlabel=r'$L$',ylabel=r'$1-f$')
            ax[1].set(xlabel=r'$L$',ylabel=r'$\langle T_e \rangle$')
            ax[0].set_xscale("log")
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")

        else:
            for f in range(len(fins)): 
                L = data[f]["L"]
                ps = data[f]["ps"]
                erode_fracs = data[f]["erode_frac"].T
                erode_times = data[f]["erode_times"].T / np.log(L)
                ax[0].plot(ps,erode_fracs,label=L,c=cmap(f/len(fins)),mec='k',lw=lw,marker='o',ms=ms)
                ax[1].plot(ps,erode_times,label=L,c=cmap(f/len(fins)),mec='k',lw=lw,marker='o',ms=ms)


elif "Ft" in mode: 
    plot_inf = False # if true, plots infidelity 

    if vary_L: 
        xlab = r"$L$"
        ylab = r"$\varepsilon$" if plot_inf else r"$F_T$" 
        leg_loc = 'lower left'
        leg_title = r'$p$'
        power_fit = not False
        for f in range(len(fins)): 
            if plot_inf:
                power = .5
                try: 
                    times = data[f]["Ts"] 
                    print("times = ",times)
                except: 
                    times = 100 * np.sqrt(data[f]["Ls"])
                ys = (1-data[f]["Ft"])/ times # this is essentially the logical error rate / cycle ("essentially" since we're not getting it in the steady state)
            else: 
                ys = data[f]["Ft"] / data[f]["Ls"]
            if args.error_rate: 
                ys = data[f]["logical_error_rate"]

            ax.plot(data[f]["Ls"],ys,label=data[f]["ps"][0],c=cmap(f/len(fins)),mec='k',lw=lw,marker='o',ms=ms)
            if power_fit: 
                mp = int(len(data[f]["Ls"])/2)
                fits = np.polyfit(np.log(data[f]["Ls"][:mp]),np.log(ys)[:mp],1)
                ax.plot(data[f]["Ls"],np.exp(fits[1]) * data[f]["Ls"]**fits[0],c=cmap(f/len(fins)),ls='--',lw=3*lw,alpha=.7,label=r'$%.2f$'%(-fits[0]))

        ax.set_yscale('log')
        ax.set_xscale('log')

        if power_fit: 
            ax.set_xscale('log')

    else: # vary p  
        xlab = r"$p$"
        leg_loc = 'lower left'
        leg_title = r'$L$'
        # leg_loc = 'lower left'
        normalize = False 
        if plot_inf:
            ylab = r"$1-F_T$" if not normalize else r"$(1-F_T)/T$" 
            ax.set_title(ylab)
        else: 
            ylab = r"$F_T$"
        if args.error_rate: 
            ylab = r'$\varepsilon$'

        for f in range(len(fins)): 
            if plot_inf:
                ys = (1-data[f]["Ft"]) / (data[f]["L"] if normalize else 1)
            else: 
                ys = data[f]["Ft"]
            if args.error_rate: 
                ys = data[f]["logical_error_rate"]
                # T = data[f]["Ts"][1]
                # print(T)
                # ys = (1-(np.abs((2 * np.sqrt(data[f]["Ft"]) -1)))**(1/T))/2
            ax.plot(data[f]["ps"],ys,label=data[f]["L"],c=cmap(f/len(fins)),mec='k',lw=lw,marker='o',ms=ms)
        # ax.set_yscale('log'); ax.set_xscale('log')

elif "trel" in mode: 
    if vary_L: 
        for f in range(len(fins)): 
            Ls = data[f]["Ls"]
            ys = data[f]["trels"]
            samps = data[f]["samps"]
            print("p = ",data[f]["p"])

            ax.plot(Ls,ys,label=data[f]["p"],c=cmap(f/len(fins)),mec='k',lw=lw,marker='o',ms=ms)

        ax.set_xscale('log')
        leg_loc = 'best'

    else: 
        uselogps = not True 
        ylab = r'$t_{\sf mem}$'
        fits = np.zeros(len(fins))
        Ls = np.zeros(len(fins))
        for f in range(len(fins)): 
            # get linear fit between ln(1/p) and trel: 
            Ls[f] = data[f]["L"]
            if uselogps:
                xs = np.log(1/data[f]["ps"])
                xlab = r'$\ln(1/p)$'
            else: 
                xs = data[f]["ps"]
                xlab = r'$p$'
                ax.set_xscale('log')
            ys = data[f]["trels"]

            ax.plot(xs,ys,label=int(data[f]["L"] / (2 if "tlv" in fins[0] else 1)),c=cmap(f/len(fins)),mec='k',lw=lw,marker='o',ms=ms)
            try: 
                print(data[f]["qs"])
            except: 
                pass 
            
            mp = int(len(xs)/2) # do the fit for the data in the second half of the curve
            fits[f] = np.polyfit(np.log(1/xs[:mp]),np.log(ys[:mp]),1)[0]


        ax.set_yscale('log')
        # ax.set_xscale('log')

        # add inset plot showing the slope of the linear fits
        add_inset = False 
        if add_inset:
            ax_inset = ax.inset_axes([0.15,0.6,0.35,0.35])
            ax_inset.plot(Ls,fits,marker='o',ms=ms,mec='k',c='r',lw=lw)
            ax_inset.set(xlabel=r'$L$',ylabel=r'$\alpha$')


# set titles and labels etc. 
try: 
    for i in range(nplots): 
        ax[i].legend(title=leg_title,fontsize=legend_fontsize,title_fontsize=legend_title_fontsize,loc=leg_loc)
except:
    ax.set(xlabel=xlab,ylabel=ylab)
    ax.set_title(title)
    ax.legend(title=leg_title,fontsize=legend_fontsize,title_fontsize=legend_title_fontsize,loc=leg_loc)

plt.show()

if args.save: 
    fig.savefig(args.save,bbox_inches='tight')
    print("saved to ",args.save)