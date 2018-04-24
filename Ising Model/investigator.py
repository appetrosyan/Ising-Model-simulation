#! /usr/local/bin/python3
from AnalysisTools import Simulation, parmap, tqdm, t_c, use_disk, \
    smart_duration, save_plot, pretty_pm, multi_run, auto_cov
from numpy import linspace, sqrt, append, exp, array, argmax, arange
from scipy.interpolate import splrep, splev
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


# ------------------------------------------------------------------------


def plot_time_evolution(simulations):
    fig, ax = plt.subplots(nrows=2, sharex='col')
    for s in simulations:
        ax[0].plot(s.times(), s.mean_magnetizations(),
                   label='T=' + str(s.temperature))
    ax[0].set_ylabel('mean magnetisation / spin')
    for s in simulations:
        ax[1].plot(s.times(), s.mean_energies(),
                   label='T=' + str(s.temperature))
    ax[1].set_ylabel('mean energy / (J/spin)')
    plt.xlabel('time steps')


def investigate_time_evolution(temps=None, grain_sizes=None, h=0):
    if temps is None:
        temps = [0.12, 0.56, 2.1, 2.3, 3.5]
    if grain_sizes is None:
        grain_sizes = [0, 3, 10]
    for g in grain_sizes:
        sims = [
            Simulation(temperature=t, grain_size=g, duration=(300 +
                                                              50*g),
                       magnetic_field=h)
            for t in temps]
        plot_time_evolution(sims)
        save_plot('Evolution from chequerboard pattern of grain size '
                  '' + str(g) + 'x' + str(g))


# -------------------------------------------------------------------------

def magnetisation_fit(x, centre, breadth, ):
    return 1/(exp(-(x - centre)/breadth) - 1)


def investigate_magnetisation_et_al(lattice_sizes, temps=None,
                                    h=0, **kwargs):
    if temps is None:
        temps = linspace(1.2, 3, 32)
    fig, ax = plt.subplots(nrows=2, sharex='col')
    ax[0].set_ylabel('mean magnetisation / spin')
    ax[1].set_ylabel('mean energy / (J/spin)')
    ax[0].axvline(x=t_c, color='k', ls='--')
    ax[1].axvline(x=t_c, color='k', ls='--')
    res = []
    for l in tqdm(lattice_sizes):
        simulations = [Simulation(lattice_size=l, temperature=t,
                                  duration=smart_duration(l, t),
                                  dry_run=True) for t in temps]
        data = parmap(lambda s: multi_run(s, **kwargs), tqdm(simulations))
        terminal_magnetisations = array(data).T[0]
        sigma_magnetisations = array(data).T[1]
        terminal_energies = array(data).T[2]
        ax[0].errorbar(temps, terminal_magnetisations,
                       yerr=sigma_magnetisations, fmt='-',
                       label='N = ' + str(l))
        ax[1].plot(temps, terminal_energies, '-', label='N =' + str(l))
        popt, _ = curve_fit(magnetisation_fit, temps,
                            terminal_magnetisations)
        res.append(popt[0])
    plt.xlabel('temperature / J')
    plt.xticks(
        append(linspace(min(temps), t_c, 3), linspace(t_c, max(temps), 3)))
    save_plot('Critical temperature from magnetisation')
    return res


# -------------------------------------------------------------------------

def investigate_autocorrelation(lattice_sizes, max_t: int = 10**6,
                                temps=None, take_every=3):
    if temps is None:
        temps = [2.2, t_c, 2.3]
    time = linspace(0, max_t, 50, dtype=int)
    fig, axes = plt.subplots(nrows=len(lattice_sizes[0::take_every]),
                             sharex='col', sharey='col')

    for ax, l in zip(axes, lattice_sizes[0::take_every]):
        for temp in tqdm(temps):
            s = Simulation(duration=max_t, lattice_size=l,
                           temperature=temp)
            autocov = parmap(
                lambda t: auto_cov(s.mean_magnetizations(), t), time)
            norm = auto_cov(s.mean_magnetizations(), 0)
            ax.plot(time, autocov/norm, '+', label=s.temperature, )
            ax.axhline(y=1/exp(1), ls='--')
            ax.set_yticks([0, 1/exp(1), 1], )
            ax.set_yticklabels(['0', r'$\frac{1}{e}$', '1'])
            # ax.set_yscale('log')
            ax.set_ylabel(r'$a(\tau)$' + r'N = ' + str(l))
    plt.xlabel(r'$\tau$')
    save_plot('Autocorrelation', legend_loc='upper right')


# -------------------------------------------------------------------------


def theoretical_capacity(x, a, b, c, d):
    return c + (b*10*x)/((x - a)**2 + d/10)


def investigate_heat_capacity(lattice_sizes, temps=None, **kwargs):
    if temps is None:
        temps = append(linspace(1.6, 4.25, 12), linspace(1.8, 2.8, 24))
        temps.sort()

    fig, axes = plt.subplots(nrows=len(lattice_sizes), sharex='col',
                             sharey='col')
    crit_temps = []
    for l, ax in zip(tqdm(lattice_sizes, desc='Lattice sizes'), axes):
        # get upscaled and smoothed capacity vs temperature data.
        up_temps, caps = plot_capacity(ax, l, temps, **kwargs)
        main_peak = argmax(caps)
        ax.plot(up_temps[main_peak], caps[main_peak], 'bo', )
        crit_temps.append(up_temps[main_peak])
    fig.set_size_inches(10.5, 10.5)
    plt.xlabel('temperature / arb. u.')
    save_plot('Heat capacity', )
    return crit_temps


def plot_capacity(ax, l, temps, **kwargs):
    """
    Plot the heat capacity of simulations at given temperature and lattice
    size.

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
    cap_errs, caps = generate_capacity_data(l, temps, **kwargs, )
    ax.errorbar(temps, caps, fmt='b+', yerr=cap_errs,
                label='N = ' + str(l))
    # smooth data
    smoothed = splrep(temps, caps, s=1/256)
    up_temps = linspace(min(temps), max(temps), 1000)
    ax.plot(up_temps, splev(up_temps, smoothed),
            label='N = ' + str(l) + ' smoothed')
    ax.set_ylabel(r'$C_V$/ arb. u.')
    ax.axvline(x=t_c, ls='-.', color='k')
    ax.set_xticks(
        append(linspace(min(temps), t_c, 6), linspace(t_c, max(temps), 6)))
    ax.legend(loc='best')
    return up_temps, splev(up_temps, smoothed)


def generate_capacity_data(l, temps, **kwargs, ):
    simulations = [Simulation(lattice_size=l, temperature=t, dry_run=True)
                   for t in temps]
    data = parmap(lambda s: multi_run(s, **kwargs),
                  tqdm(simulations, desc='Temperature samples'))
    sigmas = array(data).T[3]
    meta_sigmas = array(data).T[4]
    capacities = array(
        [sigma**2/(temp**2)*l**2 for temp, sigma in zip(temps, sigmas)])
    sigma_capacities = capacities[:]*meta_sigmas[:]
    return sigma_capacities, capacities


# -------------------------------------------------------------------------

def finite_size_scale(N, t_inf, a, v):
    return t_inf + a*(N**(-1/v))


def investigate_finite_size_scaling(critical_temperatures, lattice_sizes,
                                    source='capacity', **kwargs):
    """
    Test the hypothesis, that the critical temperature scales with the
    lattice size as $T_c (N) = T_C(\inf) + a N^{-\frac{1}{\nu}}$
    """
    if critical_temperatures is None:
        critical_temperatures = investigate_heat_capacity(lattice_sizes,
                                                          **kwargs)
        source = 'capacity'
    args, cov = curve_fit(finite_size_scale, lattice_sizes,
                          critical_temperatures)
    data = (args, cov)
    lbl = r'$ ' + pretty_pm(data,
                            0) + '$' + r' $' + r'+ \frac{' + pretty_pm(
        data, 1) + r'}{N^ {-1/' + pretty_pm(data, 2) + r'}}' + r'$'
    plt.plot(lattice_sizes, critical_temperatures, 'b+', label='data')
    plt.plot(lattice_sizes, finite_size_scale(lattice_sizes, *args), 'r-',
             label=lbl)
    plt.ylabel('critical temperature / arb. u.')
    plt.xlabel('Lattice size')
    save_plot('Finite size scaling, from ' + source)
    return args[0], sqrt(cov[0, 0])


# -------------------------------------------------------------------------



# -------------------------------------------------------------------------


sizes = arange(16, 120, 10)
# investigate_time_evolution()
# investigate_autocorrelation(sizes, max_t=1*10**2)
critical_temps_M = investigate_magnetisation_et_al(lattice_sizes=sizes[::3], re_runs=4)
# WARNING: this can take upwards of half an hour to run on Ryzen 5
# 1600.
# critical_temps = investigate_heat_capacity(lattice_sizes=sizes,
#                                            take_last=100000, re_runs=1)
# temp_inf_m = investigate_finite_size_scaling(critical_temps_M,
#                                              sizes[::3],
#                                              source='Magnetisation')
# temp_inf = investigate_finite_size_scaling(critical_temps, sizes)
# print(temp_inf)
# print(temp_inf_m)
# print((t_c - temp_inf[0])/temp_inf[1], ' Standard errors away.  '
#                                        'capacity')
# print((t_c - temp_inf_m[1]/temp_inf[1]), ' Standard errors away. '
#                                          'Magnetisation')