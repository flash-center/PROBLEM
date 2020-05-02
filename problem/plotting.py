import matplotlib as mpl
mpl.use('qt5agg')
import matplotlib.pyplot as plt
plt.ion()


class PlotPROBLEM:
    N = 2

    def __init__(self, x, y, F):
        f, ax = plt.subplots(1,1)
        ax.pcolormesh(x, y, F)
        self.points = (
            ax.plot(x[::self.N,::self.N], y[::self.N,::self.N], '-', lw=0.1, color='w'),
            ax.plot(x[::self.N,::self.N].T, y[::self.N,::self.N].T, '-', lw=0.1, color='w')
        )
        plt.draw()
        plt.pause(0.01)
        plt.show()

    def update(self, x, y, title):
        for (p, px, py) in zip(self.points[0], x[::self.N, ::self.N].T, y[::self.N, ::self.N].T):
            p.set_data(px, py) 
        for (p, px, py) in zip(self.points[1], x[::self.N, ::self.N], y[::self.N, ::self.N]):
            p.set_data(px, py)
        plt.title(title)
        plt.draw()
        plt.pause(0.01)
