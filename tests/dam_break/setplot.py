"""setplot.py for the dam-break test: h and velocity u at each frame."""

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

    plotfigure = plotdata.new_plotfigure(name='depth', figno=0)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [0.0, 0.6]
    plotaxes.title = 'Water depth h'
    plotaxes.xlabel = 'x (m)'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-'
    plotitem.color = 'b'

    plotfigure = plotdata.new_plotfigure(name='velocity', figno=1)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Velocity u'
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
