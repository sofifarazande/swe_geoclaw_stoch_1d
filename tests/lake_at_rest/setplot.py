"""
setplot.py for the lake-at-rest test: plot eta=h+z (should stay flat),
the bed z, and velocity u (should stay at zero) at each frame.
"""

def eta(current_data):
    q = current_data.q
    return q[0,:] + q[2,:]

def bed(current_data):
    return current_data.q[2,:]

def velocity(current_data):
    q = current_data.q
    h = q[0,:]
    hu = q[1,:]
    u = 0.0*h
    wet = h > 1.e-6
    u[wet] = hu[wet] / h[wet]
    return u


def setplot(plotdata):

    plotdata.clearfigures()

    plotfigure = plotdata.new_plotfigure(name='surface_and_bed', figno=0)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [0.0, 0.45]
    plotaxes.title = 'Free surface (eta=h+z) and bed (z) - lake at rest'
    plotaxes.xlabel = 'x (m)'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = eta
    plotitem.plotstyle = '-'
    plotitem.color = 'b'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = bed
    plotitem.plotstyle = '-'
    plotitem.color = 'k'

    plotfigure = plotdata.new_plotfigure(name='velocity', figno=1)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Velocity u (should stay at 0)'
    plotaxes.xlabel = 'x (m)'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = velocity
    plotitem.plotstyle = '-'
    plotitem.color = 'r'

    plotdata.printfigs = True
    plotdata.print_format = 'png'
    plotdata.print_framenos = 'all'
    plotdata.print_fignos = 'all'
    plotdata.html = True
    plotdata.latex = False

    return plotdata
