# The Normalize class is largely based on code provided by Sarah Graves.

import numpy as np
import numpy.ma as ma

import matplotlib.cbook as cbook
from matplotlib.colors import Normalize

import pylab as plt
import numpy as np
from scipy.ndimage import gaussian_filter


class MyNormalize(Normalize):
    '''
    A Normalize class for imshow that allows different stretching functions
    for astronomical images.
    '''

    def __init__(self, stretch='linear', exponent=5, vmid=None, vmin=None,
                 vmax=None, clip=False):
        '''
        Initalize an APLpyNormalize instance.

        Optional Keyword Arguments:

            *vmin*: [ None | float ]
                Minimum pixel value to use for the scaling.

            *vmax*: [ None | float ]
                Maximum pixel value to use for the scaling.

            *stretch*: [ 'linear' | 'log' | 'sqrt' | 'arcsinh' | 'power' ]
                The stretch function to use (default is 'linear').

            *vmid*: [ None | float ]
                Mid-pixel value used for the log and arcsinh stretches. If
                set to None, a default value is picked.

            *exponent*: [ float ]
                if self.stretch is set to 'power', this is the exponent to use.

            *clip*: [ True | False ]
                If clip is True and the given value falls outside the range,
                the returned value will be 0 or 1, whichever is closer.
        '''

        if vmax < vmin:
            raise Exception("vmax should be larger than vmin")

        # Call original initalization routine
        Normalize.__init__(self, vmin=vmin, vmax=vmax, clip=clip)

        # Save parameters
        self.stretch = stretch
        self.exponent = exponent

        if stretch == 'power' and np.equal(self.exponent, None):
            raise Exception(
                "For stretch=='power', an exponent should be specified")

        if np.equal(vmid, None):
            if stretch == 'log':
                if vmin > 0:
                    self.midpoint = vmax / vmin
                else:
                    raise Exception(
                        "When using a log stretch, if vmin < 0, then vmid has to be specified")
            elif stretch == 'arcsinh':
                self.midpoint = -1. / 30.
            else:
                self.midpoint = None
        else:
            if stretch == 'log':
                if vmin < vmid:
                    raise Exception(
                        "When using a log stretch, vmin should be larger than vmid")
                self.midpoint = (vmax - vmid) / (vmin - vmid)
            elif stretch == 'arcsinh':
                self.midpoint = (vmid - vmin) / (vmax - vmin)
            else:
                self.midpoint = None

    def __call__(self, value, clip=None):

        #read in parameters
        method = self.stretch
        exponent = self.exponent
        midpoint = self.midpoint

        # ORIGINAL MATPLOTLIB CODE

        if clip is None:
            clip = self.clip

        if cbook.iterable(value):
            vtype = 'array'
            val = ma.asarray(value).astype(np.float)
        else:
            vtype = 'scalar'
            val = ma.array([value]).astype(np.float)

        self.autoscale_None(val)
        vmin, vmax = self.vmin, self.vmax
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin == vmax:
            return 0.0 * val
        else:
            if clip:
                mask = ma.getmask(val)
                val = ma.array(np.clip(val.filled(vmax), vmin, vmax),
                               mask=mask)
            result = (val - vmin) * (1.0 / (vmax - vmin))

            # CUSTOM APLPY CODE

            # Keep track of negative values
            negative = result < 0.

            if self.stretch == 'linear':

                pass

            elif self.stretch == 'log':

                result = ma.log10(result * (self.midpoint - 1.) + 1.) \
                    / ma.log10(self.midpoint)

            elif self.stretch == 'sqrt':

                result = ma.sqrt(result)

            elif self.stretch == 'arcsinh':

                result = ma.arcsinh(result / self.midpoint) \
                    / ma.arcsinh(1. / self.midpoint)

            elif self.stretch == 'power':

                result = ma.power(result, exponent)

            else:

                raise Exception("Unknown stretch in APLpyNormalize: %s" %
                                self.stretch)

            # Now set previously negative values to 0, as these are
            # different from true NaN values in the FITS image
            result[negative] = -np.inf

        if vtype == 'scalar':
            result = result[0]

        return result

    def inverse(self, value):

        # ORIGINAL MATPLOTLIB CODE

        if not self.scaled():
            raise ValueError("Not invertible until scaled")

        vmin, vmax = self.vmin, self.vmax

        # CUSTOM APLPY CODE

        if cbook.iterable(value):
            val = ma.asarray(value)
        else:
            val = value

        if self.stretch == 'linear':

            pass

        elif self.stretch == 'log':

            val = (ma.power(10., val * ma.log10(self.midpoint)) - 1.) / \
                (self.midpoint - 1.)

        elif self.stretch == 'sqrt':

            val = val * val

        elif self.stretch == 'arcsinh':

            val = self.midpoint * \
                ma.sinh(val * ma.arcsinh(1. / self.midpoint))

        elif self.stretch == 'power':

            val = ma.power(val, (1. / self.exponent))

        else:

            raise Exception("Unknown stretch in APLpyNormalize: %s" %
                            self.stretch)

        return vmin + val * (vmax - vmin)


class DraggableColorbar(object):
    def __init__(self, cbar, mappable):
        self.cbar = cbar
        self.mappable = mappable
        self.press = None
        self.cycle = sorted(
            [i for i in dir(plt.cm) if hasattr(getattr(plt.cm, i), 'N')])
        self.index = self.cycle.index(cbar.get_cmap().name)

    def connect(self):
        """connect to all the events we need"""
        self.cidpress = self.cbar.patch.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.cbar.patch.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.cbar.patch.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)
        self.keypress = self.cbar.patch.figure.canvas.mpl_connect(
            'key_press_event', self.key_press)

    def on_press(self, event):
        """on button press we will see if the mouse is over us and store some data"""
        if event.inaxes != self.cbar.ax:
            return
        self.press = event.x, event.y

    def key_press(self, event):
        if event.key == 'down':
            self.index += 1
        elif event.key == 'up':
            self.index -= 1
        if self.index < 0:
            self.index = len(self.cycle)
        elif self.index >= len(self.cycle):
            self.index = 0
        cmap = self.cycle[self.index]
        self.cbar.set_cmap(cmap)
        self.cbar.draw_all()
        self.mappable.set_cmap(cmap)
        self.mappable.get_axes().set_title(cmap)
        self.cbar.patch.figure.canvas.draw()

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.press is None:
            return
        if event.inaxes != self.cbar.ax:
            return
        xprev, yprev = self.press
        dx = event.x - xprev
        dy = event.y - yprev
        self.press = event.x, event.y
        # print 'x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f'%(x0, xpress, event.xdata, dx, x0+dx)
        scale = self.cbar.norm.vmax - self.cbar.norm.vmin
        perc = 0.03
        if event.button == 1:
            self.cbar.norm.vmin -= (perc * scale) * np.sign(dy)
            self.cbar.norm.vmax -= (perc * scale) * np.sign(dy)
        elif event.button == 3:
            self.cbar.norm.vmin -= (perc * scale) * np.sign(dy)
            self.cbar.norm.vmax += (perc * scale) * np.sign(dy)
        self.cbar.draw_all()
        self.mappable.set_norm(self.cbar.norm)
        self.cbar.patch.figure.canvas.draw()

    def on_release(self, event):
        """on release we reset the press data"""
        self.press = None
        self.mappable.set_norm(self.cbar.norm)
        self.cbar.patch.figure.canvas.draw()

    def disconnect(self):
        """disconnect all the stored connection ids"""
        self.cbar.patch.figure.canvas.mpl_disconnect(self.cidpress)
        self.cbar.patch.figure.canvas.mpl_disconnect(self.cidrelease)
        self.cbar.patch.figure.canvas.mpl_disconnect(self.cidmotion)


np.random.seed(1111)

# Create empty image
nx, ny = 256, 256
image = np.zeros((ny, nx))

# Set number of stars
n = 10000

# Generate random positions
r = np.random.random(n) * nx
theta = np.random.uniform(0., 2. * np.pi, n)

# Generate random fluxes
f = np.random.random(n) ** 2

# Compute position
x = nx / 2 + r * np.cos(theta)
y = ny / 2 + r * np.sin(theta)

# Add stars to image
# ==> First for loop and if statement <==
for i in range(n):
    if x[i] >= 0 and x[i] < nx and y[i] >= 0 and y[i] < ny:
        image[int(y[i]), int(x[i])] += f[i]

# Convolve with a gaussian
image = gaussian_filter(image, 1)
# Add noise
image += np.random.normal(3., 0.01, image.shape)
img = plt.imshow(image, cmap="jet")
cbar = plt.colorbar(format='%05.2f')
cbar.set_norm(MyNormalize(vmin=image.min(),
                          vmax=image.max(), stretch='linear'))
cbar = DraggableColorbar(cbar, img)
cbar.connect()

plt.show()
