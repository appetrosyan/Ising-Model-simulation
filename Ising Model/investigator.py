#! /usr/local/bin/python3
from os import mkdir
from os.path import exists
from subprocess import run, PIPE
from numpy import loadtxt, linspace, log, sqrt, append, exp,  \
    array, shape, inf, std, mean, argmax
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


class Simulation:
    def __init__(self, magnetic_field=0.0, exchange_energy=1.0,
                 lattice_size=64, temperature=1.0, duration=50,
                 grain_size=0):
        self.magnetic_field = magnetic_field
        self.exchange_energy = exchange_energy
        self.lattice_size = lattice_size
        self.temperature = temperature
        self.grain_size = grain_size
        self.duration = duration
        file_path = 'data/' + str(self) + '.csv'
        if use_disk:
            self.data = self.load_from_disk(file_path)
        else:
            self.data = self.run()

    def times(self):
        return self.data[:, 0]

    def mean_magnetizations(self):
        return self.data[:, 1]

    def mean_energies(self):
        return self.data[:, 2]

    def load_from_disk(self, file_path):
        if not exists(path=file_path):
            self.run(file_path)
        data = loadtxt(file_path)
        data_duration = data[-1, 0] + self.duration/200 + 10
        if self.duration > data_duration:
            print('duration mismatch, ' + str(self.duration) + '!=' + str(
                data_duration), end='. re-', flush=True)
            self.run(file_path)
            data = loadtxt(file_path)
        return data

    def run(self, file_path=None):
        if file_path is None:
            file_path = 'data/' + str(self) + '.csv'
        # print('running simulation... ', self.duration, end=' ', flush=True)
        command = ['./main', '-d', str(self.duration), '-t',
                   str(self.temperature), '-n', str(self.lattice_size),
                   '-j', str(self.exchange_energy), '-H',
                   str(self.magnetic_field), '-c', str(self.grain_size)]
        if self.duration > 500:
            command.append('-p')
            command.append(str(int(self.duration/200)))
        if use_disk:
            command.append('-f')
            command.append(file_path)
            # print('Done! ')
            run(command)
        else:
            r = run(command, stdout=PIPE)
            raw_data = loadtxt(r.stdout.decode().split(sep='\n'),
                               delimiter='\t')
            # print('Done! ')
            return raw_data

    def __str__(self):
        return '(H=' + str(self.magnetic_field) + ')(J=' + str(
            self.exchange_energy) + ')(T=' + str(
            self.temperature) + ')(N=' + str(
            self.lattice_size) + ')' + str(self.grain_size)

# ------------------------------------------------------------------------

def plot_time_evolution(simulations):
    fig, ax = plt.subplots(nrows=2, sharex='col')
    for s in simulations:
        ax[0].plot(s.times(), s.mean_magnetizations(),
                   label='T=' + str(s.temperature))
    ax[0].set_ylabel('mean magnetisation / arb. u.')
    for s in simulations:
        ax[1].plot(s.times(), s.mean_energies(),
                   label='T=' + str(s.temperature))
    ax[1].set_ylabel('temperature / arb. u.')
    plt.xlabel('time steps / arb. u.')


def save_plot(title):
    plt.legend(loc='best')
    plt.suptitle(title)
    file_name = title.replace(' ', '_')
    if not exists('figures/'):
        mkdir('figures/')
    plt.savefig('figures/%s.png'%file_name)
    plt.show()


# -------------------------------------------------------------------------

def investigate_time_evolution(temps=None, grain_sizes=None):
    if temps is None:
        temps = [0.12, 0.56, 2.1, 2.3, 3.5]
    if grain_sizes is None:
        grain_sizes = [0, 3, 10]
    for g in grain_sizes:
        sims = [
            Simulation(temperature=t, grain_size=g, duration=(300 + 50*g))
            for t in temps]
        plot_time_evolution(sims)
        save_plot('Evolution from chequerboard pattern of grain size '
                  '' + str(g) + 'x' + str(g))


# -------------------------------------------------------------------------


def smart_duration(temperature, multiplier=1.):
    return int(((10**3)*base_duration*multiplier)/(
                (temperature - t_c)**2 + breadth))


def investigate_temperature_dependence(temps=None, lattice_sizes=None,
                                       **kwargs):
    if temps is None:
        temps = append(linspace(1.5, 2.2, 8), linspace(2.2, 2.3, 8))
        temps = append(temps, linspace(2.3, 4, 8))
    if lattice_sizes is None:
        lattice_sizes = [50, 250, 300]
    fig, ax = plt.subplots(nrows=2, sharex='col')
    ax[0].set_ylabel('mean magnetisation / arb. u.')
    ax[1].set_ylabel('mean energy / arb. u.')
    ax[0].axvline(x=t_c, color='k', ls='--')
    ax[1].axvline(x=t_c, color='k', ls='--')
    for l in lattice_sizes:
        print(l)
        simulations = [Simulation(lattice_size=l, temperature=t,
                                  duration=smart_duration(l, t)) for t in
                       temps]
        data = [multi_run(s, **kwargs) for s in simulations]
        terminal_magnetisations = array(data).T[0]
        sigma_magnetisations = array(data).T[1]
        terminal_energies = array(data).T[2]
        ax[0].plot(temps, terminal_magnetisations, '-',
                   label='N = ' + str(l))
        ax[1].plot(temps, terminal_energies, '-', label='N =' + str(l))
    plt.xlabel('temperature / arb. u.')
    plt.xticks(
        append(linspace(min(temps), t_c, 3), linspace(t_c, max(temps), 3)))
    save_plot('Critical temperature')


# -------------------------------------------------------------------------


def theoretical_capacity(x, a, b, c, d):
    return c + (b*x)/((x - a)**2 + d)


def investigate_heat_capacity(lattice_sizes=None, temps=None, **kwargs):
    if temps is None:
        temps = linspace(1.5, 3, 60)
    if lattice_sizes is None:
        lattice_sizes = [16,32, 34, 36]

    fig, axes = plt.subplots(nrows=len(lattice_sizes), sharex='col')
    crit_temps = []
    for l, ax in zip(lattice_sizes, axes):
        crit_temps.append(fit_and_plot_capacity(ax, l, temps, **kwargs))

    fig.set_size_inches(10.5, 10.5)
    plt.xlabel('temperature / arb. u.')
    save_plot('Heat capacity')
    return crit_temps


def fit_and_plot_capacity(ax, l, temps, **kwargs):
    """
    Plot the heat capacity of simulations at given temperature and lattice size.
    Afterwards fit a lorentzian and plot.

    Parameters:
    -----------
    ax : pyplot.axis
    what to plot to.

    l : int
    lattice size

    temps: numpy.array
    Array of temperatures where to evaluate heat capacity.

    Returns:
    -------
    popt[0]: float
    most likely critical temperature.
    """
    global use_disk
    use_disk = False
    simulations = [Simulation(lattice_size=l, temperature=t, duration=2)
                   for t in temps]
    # sigmas = [stdev(s.mean_energies[:]) for s in simulations]
    sigmas = array([multi_run(s, **kwargs) for s in simulations]).T[3]
    meta_sigmas = array([multi_run(s, **kwargs) for s in simulations]).T[4]
    # print(meta_sigmas)
    Cs = [sigma**2/(temp**2)*10**3 for temp, sigma in zip(temps, sigmas)]
    C_errs = Cs[:]*meta_sigmas[:]
    try:
        popt, pcov = curve_fit(theoretical_capacity, temps*10, Cs,
                               sigma=meta_sigmas, bounds=(
                [min(temps) - .2, 0, 0, 0],
                [max(temps) + .2, inf, inf, .7]))
    except RuntimeError:
        popt = [temps[argmax(Cs)], (max(Cs) - min(Cs))/4, min(Cs),
                (max(temps) - min(temps))/4]
        print('I\'m too dumb to fit')

    ax.axvline(x=popt[0], ls='--', color='g')
    ax.plot(temps, theoretical_capacity(temps*10, *popt), 'g-',
            label='N = ' + str(l) + 'd = ' + str(popt[3]) + ' fit')
    ax.errorbar(temps*10, Cs, fmt='b.', yerr=C_errs, label='N = ' + str(l))
    ax.set_ylabel(r'C $\cdot 10^3$/ arb. u.')
    ax.axvline(x=t_c, ls='-.', color='k')
    ax.set_xticks(
        append(linspace(min(temps), t_c, 6), linspace(t_c, max(temps), 6)))
    ax.legend(loc='best')
    return popt[0]


def multi_run(sim, re_runs:int=2, take_last:int=300):
    """
    Re run the Ising model simulation multiple times, and gather statistics.

    Parameters:
    ----------
    sim: simulation
    re_runs: int
    number of times to repeat simulation
    take_last: int
    How many of the final points to take statistics over.

    Returns:
    list:

    """
    global use_disk
    use_disk = False
    sim.duration = smart_duration(sim.temperature)
    print(sim.duration)
    magnetizations = []
    sigma_magnetizations = []
    energies = []
    sigma_energies = []
    for i in range(re_runs):
        sim.data = sim.run()
        # Make each run take 50% longer, so that we
        # can see if a system is still settling
        sim.duration *= 2

        last_magnetizations = sim.mean_magnetizations()[-take_last:]
        magnetizations.append(abs(mean(last_magnetizations)))
        sigma_magnetizations.append(std(last_magnetizations))

        last_energies = sim.mean_energies()[-take_last:]
        energies.append(mean(last_energies))
        sigma_energies.append(std(last_energies))

    return [mean(magnetizations), mean(sigma_magnetizations),
            mean(energies), mean(sigma_energies),
            std(sigma_energies)/mean(sigma_energies)]


# -------------------------------------------------------------------------

def finite_size_scale(N, t_inf, a, v):
    return t_inf + a*(N**(-1/v))


def investigate_finite_size_scaling(critical_temperatures, lattice_sizes,
                                    **kwargs):
    if critical_temperatures is None:
        critical_temperatures = investigate_heat_capacity(lattice_sizes,
                                                          **kwargs)
    args, cov = curve_fit(finite_size_scale, lattice_sizes,
                                  critical_temperatures)
    plt.plot(lattice_sizes, critical_temperatures, 'b+', label='data')
    plt.plot(lattice_sizes, finite_size_scale(lattice_sizes, *args), 'r-',
             label='fit')
    plt.ylabel('critical temperature / arb. u.')
    plt.xlabel('Lattice size')
    save_plot('Finite size scaling')
    return args[0], sqrt(cov[0, 0])


# -------------------------------------------------------------------------

t_c = 2/log(1 + sqrt(2))
use_disk = False
breadth = 1
base_duration = 50
sizes = [16, 32, 40, 50, 70]
# investigate_time_evolution()
# investigate_temperature_dependence(lattice_sizes=sizes, re_runs=4)

critical_temps = investigate_heat_capacity(lattice_sizes=sizes,
                                           take_last=300)

# temp_inf = investigate_finite_size_scaling(critical_temps, sizes)

#  print(temp_inf)
# print((t_c - temp_inf[0])/ temp_inf[1], ' Standard errors away')
