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
try:
    from matplotlib.pyplot import figure,plot,xlabel,ylabel,title,show,ylim,legend,xlim,rcParams,savefig,bar,xticks,yticks,subplot,grid,setp,errorbar,loglog,semilogx,semilogy,text

    params = {'legend.fontsize':14,'figure.facecolor':'white','figure.subplot.hspace':0.,
          'axes.labelsize':16,'xtick.labelsize':14,'ytick.labelsize':14}
    rcParams.update(params)
except (ImportError,RuntimeError):
   figure,plot,xlabel,ylabel,title,show,ylim,legend,xlim,rcParams,savefig,bar,xticks,yticks,subplot,grid,setp,errorbar,loglog,semilogx,semilogy,text = unavailable('matplotlib.pyplot','figure','plot','xlabel','ylabel','title','show','ylim','legend','xlim','rcParams','savefig','bar','xticks','yticks','subplot','grid','setp','errorbar','loglog','semilogx','semilogy','text')
#end try


# savefig(savefile,format='png',bbox_inches ='tight',pad_inches=.3)
