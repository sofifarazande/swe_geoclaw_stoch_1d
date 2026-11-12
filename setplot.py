"""
setplot.py -- VisClaw plots: free surface eta = h + z with the bed z,
velocity u, and particle activity b.  q[0]=h, q[1]=hu, q[2]=z, q[3]=b.
"""

def eta(current_data):
    q = current_data.q
    return q[0,:] + q[2,:]

def velocity(current_data):
    q = current_data.q
    h = q[0,:]
    hu = q[1,:]
    u = 0.0*h
    wet = h > 1.e-6
    u[wet] = hu[wet] / h[wet]
    return u

def bed(current_data):
    return current_data.q[2,:]

def activity(current_data):
    return current_data.q[3,:]


def setplot(plotdata):

    plotdata.clearfigures()

    # ---- Figure: free surface eta and bed z together ----
    plotfigure = plotdata.new_plotfigure(name='surface_and_bed', figno=0)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Free surface (eta=h+z) and bed (z)'
    plotaxes.xlabel = 'x (m)'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = eta
    plotitem.plotstyle = '-'
    plotitem.color = 'b'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = bed
    plotitem.plotstyle = '-'
    plotitem.color = 'k'

    # ---- Figure: velocity ----
    plotfigure = plotdata.new_plotfigure(name='velocity', figno=1)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Depth-averaged velocity u'
    plotaxes.xlabel = 'x (m)'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = velocity
    plotitem.plotstyle = '-'
    plotitem.color = 'b'

    # ---- Figure: particle activity b ----
    plotfigure = plotdata.new_plotfigure(name='activity', figno=2)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Particle activity b'
    plotaxes.xlabel = 'x (m)'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = activity
    plotitem.plotstyle = '-'
    plotitem.color = 'r'

    plotdata.printfigs = True
    plotdata.print_format = 'png'
    plotdata.print_framenos = 'all'
    plotdata.print_fignos = 'all'
    plotdata.html = True
    plotdata.latex = False

    return plotdata
