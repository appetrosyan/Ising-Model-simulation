import multiprocessing
from os import mkdir
from subprocess import run, PIPE
from numpy import loadtxt, log, sqrt, mean, std, zeros, sum as np_sum, nonzero, ravel, logical_and
from os.path import exists
import matplotlib.pyplot as plt

USE_DISK = False
breadth = 1
base_duration = .75
t_c = 2/log(1 + sqrt(2))
silent=True



try:
    # This will print a nice progressbar, if tqdm is installed, else,
    # give some useful feedback
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, desc='', **kwargs):
        """Popular progressbar framework"""
        try:
            total = len(iterable)
        except TypeError:
            total = None
        n = 0
        for it in iterable:
            yield it
            n += 1
            print('\r', desc, n, '/', total)

def fun(f, q_in, q_out):
    """Helper for parmap"""
    while True:
        i, x = q_in.get()
        if i is None:
            break
        q_out.put((i, f(x)))


def parmap(f: callable, X, nprocs=multiprocessing.cpu_count()):
    """
    A multiprocessing map, that can work with member functions and
    lambda expressions.
    """
    q_in = multiprocessing.Queue(1)
    q_out = multiprocessing.Queue()
    proc = [multiprocessing.Process(target=fun, args=(f, q_in, q_out)) for
            _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()
    sent = [q_in.put((i, x)) for i, x in enumerate(X)]
    [q_in.put((None, None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]
    [p.join() for p in proc]
    return [x for i, x in sorted(res)]

class Simulation:
    def __init__(self, magnetic_field=0.0, exchange_energy=1.0,
                 lattice_size: int = 64, temperature=1.0,
                 duration: int = 50, grain_size: int = 0,
                 dry_run: bool = False):
        self.magnetic_field = magnetic_field
        self.exchange_energy = exchange_energy
        self.lattice_size = lattice_size
        self.temperature = temperature
        self.grain_size = grain_size
        self.duration = duration
        if not dry_run:
            file_path = 'data/' + str(self) + '.csv'
            if USE_DISK:
                self.data = self.load_from_disk(file_path)
            else:
                self.data = self.run()

    def times(self):
        return self.data[:, 0]

    def magnetizations(self):
        return self.data[:, 1]

    def mean_energies(self):
        return self.data[:, 2]

    def mean_variances(self):
        return self.data[:, 3]

    def load_from_disk(self, file_path):
        if not exists(path=file_path):
            self.run(file_path)
        data = loadtxt(file_path)
        data_duration = data[-1, 0] + self.duration/256 + 10
        if self.duration > data_duration:
            print('duration mismatch, ' + str(self.duration) + '!=' + str(
                data_duration), end='. re-', flush=True)
            self.run(file_path)
            data = loadtxt(file_path)
        return data

    def __str__(self):
        return '(H=' + str(self.magnetic_field) + ')(J=' + str(
            self.exchange_energy) + ')(T=' + str(
            self.temperature) + ')(N=' + str(
            self.lattice_size) + ')' + str(self.grain_size)

    def run(self, file_path=None):
        if file_path is None:
            file_path = 'data/' + str(self) + '.csv'
        command = ['./main', '-d', str(self.duration), '-t',
                   str(self.temperature), '-n', str(self.lattice_size),
                   '-j', str(self.exchange_energy), '-H',
                   str(self.magnetic_field), '-c', str(self.grain_size)]
        if self.duration > 50000:
            command.append('-p')
            command.append(str(int(self.duration/50000)))
        if USE_DISK:
            command.append('-f')
            command.append(file_path)
            run(command)
        else:
            r = run(command, stdout=PIPE)
            raw_data = loadtxt(r.stdout.decode().split(sep='\n'),
                               delimiter='\t')
            return raw_data


def save_plot(title, legend_loc='best'):
    plt.legend(loc=legend_loc)
    # plt.suptitle(title)
    file_name = title.replace(' ', '_')
    if not exists('figures/'):
        mkdir('figures/')
    plt.savefig('figures/%s.eps'%file_name)
    if not silent:
        plt.show()
    plt.figure()


def smart_duration(temperature, multiplier=1.):
    """
    A function to scale the simulation runtime based on distance to
    critical temperature.
    """
    return int(((10**5)*base_duration*multiplier)/(
            (temperature - 2.4)**2 + breadth))

def pretty_pm(data, i):
    """Helper function to format latex. """
    args, cov = data
    return '(' + '{:01.2f}'.format(args[i]) + r'\pm' + '{:01.2f}'.format(
        sqrt(cov[i, i])) + ')'


def multi_run(sim: Simulation, re_runs: int = 2, take_last: int = 300,
              multiplier: float = 1.0):
    """
    Run a given simulation multiple times.

    Parameters:
    -----------

    sim: Simulation
    the simulation to be run.
    re_runs: int
    number of repetitions.
    take_last: int
    Number of samples to take from the end of simulations to form stats

    Returns:
    List:
    """
    sim.duration = smart_duration(sim.temperature, multiplier=multiplier)
    magnetizations = []
    sigma_magnetizations = []
    energies = []
    sigma_energies = []
    for _ in range(re_runs):
        sim.data = sim.run()
        sim.duration -= take_last
        last_magnetizations = sim.magnetizations()[-take_last:]
        magnetizations.append(abs(mean(last_magnetizations)))
        sigma_magnetizations.append(std(last_magnetizations))
        last_energies = sim.mean_energies()[-take_last:]
        energies.append(mean(last_energies))
        sigma_energies.append(std(last_energies))
    return [mean(magnetizations), mean(sigma_magnetizations),
            mean(energies), mean(sigma_energies),
            std(sigma_energies)/mean(sigma_energies)]


def auto_cov(x_i, t):
    """
    Compute empirical autocovariance. 
    """
    length = len(x_i)
    # use sample mean estimate from whole series
    x_m = mean(x_i)
    # construct copies of series shifted relative to each other,
    # subtracting the mean
    end_padded_series = zeros(length + t)
    end_padded_series[:length] = x_i - x_m
    start_padded_series = zeros(length + t)
    start_padded_series[t:] = x_i - x_m
    return 1./(length - 1)*np_sum(start_padded_series*end_padded_series)

def find_root(xdata, ydata, offset):
    """
    Parameters:
    -----------
    xdata: iterable
    ydata: iterable
  
    Returns: 
    --------
    x_0 value of the most likely root. 
    """
    # 
    crossing_indices, = nonzero(ravel(logical_and(ydata[1:] <= offset,
                                           ydata[:-1] > offset)))
    index = crossing_indices[0]
    # Linear interpolation 
    x_0, y_0 = xdata[index], ydata[index]
    x_f, y_f = xdata[index+1], ydata[index+1]
    d_x, d_y = x_0 - x_f, y_0 - y_f 
    return x_0 - (y_0 - offset)/d_y * d_x
