##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  plotting.py                                                       #
#    Standard matplotlib imports for various plotting functions.     #
#    Used extensively by SimulationAnalyzer classes.                 #
#                                                                    #
#====================================================================#


from developer import unavailable
success = False
try:
    import matplotlib
    gui_envs = ['GTKAgg','TKAgg','Qt4Agg','WXAgg']
    for gui in gui_envs:
        try:
            matplotlib.use(gui,warn=False, force=True)
            from matplotlib import pyplot
            success = True
            break
        except:
            continue
        #end try
    #end for
    from matplotlib.pyplot import figure,plot,xlabel,ylabel,title,show,ylim,legend,xlim,rcParams,savefig,bar,xticks,yticks,subplot,grid,setp,errorbar,loglog,semilogx,semilogy,text

    params = {'legend.fontsize':14,'figure.facecolor':'white','figure.subplot.hspace':0.,
          'axes.labelsize':16,'xtick.labelsize':14,'ytick.labelsize':14}
    rcParams.update(params)
except (ImportError,RuntimeError):
    success = False
#end try
if not success:
    figure,plot,xlabel,ylabel,title,show,ylim,legend,xlim,rcParams,savefig,bar,xticks,yticks,subplot,grid,setp,errorbar,loglog,semilogx,semilogy,text = unavailable('matplotlib.pyplot','figure','plot','xlabel','ylabel','title','show','ylim','legend','xlim','rcParams','savefig','bar','xticks','yticks','subplot','grid','setp','errorbar','loglog','semilogx','semilogy','text')
    pyplot     = unavailable('matplotlib','pyplot')
    matplotlib = unavailable('matplotlib')
#end if


# savefig(savefile,format='png',bbox_inches ='tight',pad_inches=.3)
