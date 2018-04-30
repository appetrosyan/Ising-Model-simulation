#! /usr/local/bin/python3
from AnalysisTools import Simulation, parmap, tqdm, t_c, USE_DISK, \
    smart_duration, save_plot, pretty_pm, multi_run, auto_cov
from numpy import linspace, sqrt, append, exp, array, argmax, arange, mean, \
    polyfit, std
from scipy.interpolate import splrep, splev
from scipy.optimize import curve_fit, brentq
import matplotlib.pyplot as plt
from matplotlib import rc

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font', **{'family': 'serif', 'serif': ['Times']})
rc('text', usetex=True)
plt.rcParams["font.size"] = 14


# ------------------------------------------------------------------------


def time_plot(simulations):
    """ Plot evolution of magnetisation and mean energy."""
    fig, ax = plt.subplots(nrows=2, sharex='col')
    for s in simulations:
        ax[0].plot(s.times(), s.magnetizations(),
                   label='T=' + str(s.temperature))
    ax[0].set_ylabel(r'magnetisation / spin')
    for s in simulations:
        ax[1].plot(s.times(), s.mean_energies(),
                   label='T=' + str(s.temperature))
    ax[1].set_ylabel(r'energy per spin / J')
    plt.xlabel('time steps')


def investigate_time_evolution(temps=None, grain_sizes=None, h=0):
    """
    Plot evolution of system set to chequerboard pattern of different
    grain_sizes.

    Parameters: 
    -----------
    temps: iterable 
    
    grain_sizes: iterable
    list of sizes of the chequerboard pattern. 

    h: float
    external magnetic field applied. 
    """
    if temps is None:
        temps = [0.12, 0.56, 2.1, 2.3, 3.5]
    if grain_sizes is None:
        grain_sizes = [0, 3, 10]
    for g in grain_sizes:
        sims = [
            Simulation(temperature=t, grain_size=g, duration=(300 + 50*g),
                       magnetic_field=h) for t in temps]
        time_plot(sims)
        save_plot('Evolution from chequerboard pattern of grain size '
                  '' + str(g) + 'x' + str(g))


# -------------------------------------------------------------------------

def magnetisation_fit(x, centre, breadth, ):
    """Helper function to produce fit."""
    return 1/(exp(-(x - centre)/breadth) - 1)


def investigate_magnetisation_critical_transitions(lattice_sizes,
                                                   temps=None, **kwargs):
    """
    Plot the magnetisation and mean energy as a function of temperature. 
    Determine critical temperature by approximating the transition to a
    Fermi-Dirac distribution.
    
    Returns:
    --------
    List of critical temperatures. 
    """
    if temps is None:
        temps = linspace(1.2, 3, 32)
    fig, ax = plt.subplots(nrows=2, sharex='col')
    ax[0].set_ylabel('magnetisation / spin')
    ax[1].set_ylabel('energy per spin/ J')
    ax[0].axvline(x=t_c, color='k', ls='--')
    ax[1].axvline(x=t_c, color='k', ls='--')
    res = []
    for l in tqdm(lattice_sizes):
        simulations = [Simulation(lattice_size=l, temperature=t,
                                  duration=smart_duration(l, t),
                                  dry_run=True) for t in temps]
        data = parmap(lambda s: multi_run(s, **kwargs), tqdm(simulations))
        terminal_magnetisations = array(data).T[0]
        # The data.T[1] holds the relative errors, not absolute ones.
        sigma_magnetisations = array(data).T[1]
        sigma_magnetisations = terminal_magnetisations*sigma_magnetisations
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
                                temps=None, ):
    """
    Plots autocorrelation for different temperatures and sizes as a 
    function of time. 

    Parameters:
    -----------
    lattice_sizes: iterable
    
    max_t: int
    Maximum time to run simulations for. 

    temps: iterable

    take_every: int
    Only run for 'take_every'-th lattice size. 

    Returns:
    --------
    List $\tau_e$ , labelled by the temperature, and size. 
    
    """
    if temps is None:
        temps = linspace(1.8, 3.0, 50)
    time = linspace(0, max_t, 50, dtype=int)
    fig, axes = plt.subplots(nrows=len(lattice_sizes), sharex='col',
                             sharey='col')
    tau_e, print_interval = [], int(len(temps)/5)
    for ax, l in zip(axes, tqdm(lattice_sizes)):
        t_e, counter = [], 0
        for temp in tqdm(temps):
            s = Simulation(duration=max_t, lattice_size=l,
                           temperature=temp)
            autocorr, sigma_autocorr = generate_autocorr_data(s, time)
            t_e.append(find_tau_e(autocorr, time, max_t))
            if counter % print_interval == 0:
                plot_autocorrelation(autocorr, sigma_autocorr, ax, s, time)
            counter += 1
        tau_e.append(t_e)
    tau_e = array(tau_e)
    plt.xlabel(r'$\tau$')
    fig.set_size_inches(10.5, 10.5)
    save_plot('Autocorrelation', legend_loc='upper right')
    plt.errorbar(temps, mean(tau_e, axis=0), yerr=std(tau_e, axis=0),
                 fmt='-')
    plt.axvline(x=t_c)
    plt.xlabel('temperature / (J / $k_B$)')
    plt.ylabel(r'$\tau_e$ / time steps')
    save_plot('Coherence lifetime vs. temperature')
    plt.errorbar(lattice_sizes, mean(tau_e, axis=1),
                 yerr=std(tau_e, axis=1), fmt='.')
    plt.xlabel('Lattice size / arb. units')
    plt.ylabel(r'$\tau_e$ / time steps')
    save_plot('Coherence lifetime vs. lattice size')
    return tau_e


def find_tau_e(autocorr, time, max_t):
    smoothed = splrep(time, autocorr)
    root = brentq(lambda x: splev(x, smoothed) - 1/exp(1), 1.0, max_t)
    return root


def generate_autocorr_data(s, time):
    data = []
    for i in range(3):
        autocov = parmap(lambda t: auto_cov(s.magnetizations(), t), time)
        norm = auto_cov(s.magnetizations(), 0)
        data.append(autocov/norm)
    arr = array(data)
    return mean(arr, axis=0), std(arr, axis=0)


def plot_autocorrelation(autocorr, sigma_autocorr, ax, s, time):
    ax.errorbar(time, autocorr, yerr=sigma_autocorr, fmt='-',
                label='T = {:01.3}'.format(s.temperature))
    ax.axhline(y=1/exp(1), ls='--')
    ax.set_yticks([0, 1/exp(1), 1], )
    ax.set_yticklabels(['0', r'1/e', '1'])
    ax.set_ylabel('N = ' + str(s.lattice_size))


# -------------------------------------------------------------------------


def investigate_heat_capacity(lattice_sizes, temps=None, skip=3, **kwargs):
    """
    Plot the heat capacity versus temperature for different lattice sizes.
    Find the critical temperature for each size by doing a spline smoothing
    and finding the peak. 

    Parameters: 
    -----------
    lattice_sizes: iterable

    temps: iterable 

    **kwargs: keyword_arguments
    Arguments to be passed to multi_run()

    Returns: 
    List of critical temperatures. 
    """
    if temps is None:
        temps = append(linspace(1.6, 4.25, 12), linspace(1.8, 2.8, 24))
        temps.sort()

    fig, _ = plt.subplots()
    crit_temps, counts = [], 0
    for l in tqdm(lattice_sizes, desc='Lattice sizes'):
        for t in temps:
            pass
        caps, cap_errs = generate_data(temps, l, **kwargs, )
        up_temps, up_caps = upscale(caps, temps)
        main_peak = argmax(up_caps)
        crit_temps.append(up_temps[main_peak])
        counts += 1
        if counts % skip == 0:
            plt.errorbar(temps, caps, fmt='+', yerr=cap_errs,
                         label='N = ' + str(l))
            plt.plot(up_temps[main_peak], up_caps[main_peak], 'bo', )
            plt.plot(up_temps, up_caps,
                     label='N = ' + str(l) + ' smoothed')

            plt.xticks(append(linspace(min(temps), t_c, 3),
                              linspace(t_c, max(temps), 4)))
    plt.ylabel(r'$C_V/ k_B$')
    plt.axvline(x=t_c, ls='-.', color='k')
    fig.set_size_inches([10.5, 10.5])
    plt.xlabel('temperature / $J/ k_B$')
    save_plot('Heat capacity', )
    return crit_temps


def generate_data(temps: iter, size=32, hs=0, r1=3, r2=4, **kwargs, ):
    global use_disk
    use_disk = False
    simulations = [
        Simulation(magnetic_field=hs, lattice_size=size, temperature=t,
                   dry_run=True) for t in temps]
    data = parmap(lambda s: multi_run(s, **kwargs),
                  tqdm(simulations, desc='Temperature samples'))
    sigmas = array(data).T[r1]
    meta_sigmas = array(data).T[r2]
    if r1 == 3:
        data = array(
            [sigma**2/(temp**2)*size**2 for temp, sigma in zip(temps, sigmas)])
    elif r1 == 1:
        data = array(
            [sigma**2/(temp**2)*size**2 for temp, sigma in zip(temps, sigmas)])
    else:
        data = sigmas
    sigma_data = data[:]*meta_sigmas[:]
    return data, sigma_data


def upscale(data, temps):
    smoothed = splrep(temps, data, s=1/256)
    up_temps = linspace(min(temps), max(temps), 1000)
    return up_temps, splev(up_temps, smoothed)


# -------------------------------------------------------------------------


def finite_size_scale(x, t_inf, a, v):
    """Helper function for producing a fit. """
    return t_inf + a*(x**(-1/v))


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
    plt.ylabel('$T_C / (J/k_B)$')
    plt.xlabel('lattice size')
    save_plot('Finite size scaling, from ' + source)
    return args[0], sqrt(cov[0][0])


# -------------------------------------------------------------------------


def investigate_chi(temps, sizes):
    """
    Produce a plot of magnetic susceptibility vs temperature.
    Returns list of critical transition temperatures.
    """
    plot_magnetic_response(temps)
    t_cs = []
    for size in tqdm(sizes):
        chi_data = array([chi(t, size) for t in tqdm(temps)]).T
        chis, chierrs = chi_data[0], chi_data[1]
        plt.errorbar(temps, chis, yerr=chierrs, label=str(size), fmt='+')
        t_cs.append(temps[argmax(chis)])
    plt.axvline(x=t_c)
    plt.ylabel('susceptibility / spin')
    plt.xlabel('temperature / (J/$k_B$)')
    save_plot('Susceptibility vs Temperature')
    return t_cs


def plot_magnetic_response(temps):
    """Plot magnetisation vs H, for different temperatures."""
    counts = 0
    p = int(len(temps)/6)
    for t in tqdm(temps):
        chi(t, sz=64, plot_q=(counts%p == 0))
        counts += 1
    plt.ylabel('Magnetisation / spin')
    plt.xlabel('Applied field / (J/$\mu$)')
    save_plot('Magnetisation vs. external field')


def chi(temp, sz=64, hs=None, plot_q=False):
    """
    Compute magnetic susceptibility chi, for given size and temperature.
    """
    if hs is None:
        hs = linspace(0.01, 0.1, 20)
    # We don't want the system to equilibrate, because then susceptibility
    # would be infinite.
    tau = smart_duration(temp, multiplier=.01)
    simulations = parmap(
        lambda h: Simulation(duration=tau, temperature=temp,
                             lattice_size=sz, magnetic_field=h, ), hs)
    # grain_size=1 avoids hysteresis.
    ms = parmap(lambda s: mean(s.magnetizations()), simulations)
    popt, pcov = polyfit(hs, ms, 1, cov=True)
    if plot_q:
        plt.plot(hs, ms, '.', label='T = {:0.3f}'.format(temp))
        plt.plot(hs[:], float(popt[0])*hs[:] + float(popt[1]), '-',
                 label=r'$\chi = $' '{:0.3f}'.format(popt[0]))
    return [popt[0], pcov[0][0]]


# -------------------------------------------------------------------------


def investigate_magnetic_fluctuations(hs, temps=None, skip=3, **kwargs):
    """
    Plot the magnetic  susceptibility versus temperature for different
    lattice sizes. Find the critical temperature for each H by spline
    interpolation and finding the maximum.

    Parameters:
    -----------
    hs: iterable

    temps: iterable

    **kwargs: keyword_arguments
    Arguments to be passed to multi_run()

    Returns:
    List of critical temperatures.
    """
    if temps is None:
        temps = append(linspace(1.6, 3.9, 12), linspace(1.8, 2.8, 24))
        temps.sort()

    fig, _ = plt.subplots()
    crit_temps, counts = [], 0
    for h in tqdm(hs, desc='External fields.'):
        # get upscaled and smoothed capacity vs temperature data.
        mags, _ = generate_data(temps, hs=h, r1=1, r2=1, **kwargs)
        up_temps, up_mags = upscale(mags, temps)
        main_peak = argmax(up_mags)
        plt.plot(up_temps[main_peak], up_mags[main_peak], 'bo', )
        crit_temps.append(up_temps[main_peak])
        counts += 1
        if counts%skip == 0:
            plt.plot(temps, mags, 'b+', label=r'$\vec{H}$ = ' + str(h))
            plt.plot(up_temps, up_mags,
                     label='N = ' + str(h) + ' smoothed')
            plt.ylabel(r'$\chi/ (\mu/J)')
            plt.axvline(x=t_c, ls='-.', color='k')
            plt.xticks(append(linspace(min(temps), t_c, 3),
                              linspace(t_c, max(temps), 6)))
    fig.set_size_inches(8.5, 8.5)
    fig.tight_layout()
    plt.xlabel('temperature / (J/$k_B$')
    save_plot('Magnetic susceptibility')
    return crit_temps

# -------------------------------------------------------------------------

sizes = arange(16, 120, 10)
temperatures = linspace(1.0, 4, 40)
ext_fields = linspace(0.5, 1.2, 5)


def sanity_checks():
    investigate_time_evolution()
    investigate_autocorrelation(sizes, max_t=1*10**2)


def fss_magnetisation():
    critical_temps_M = investigate_magnetisation_critical_transitions(
        lattice_sizes=sizes, re_runs=3)
    temp_inf_m = investigate_finite_size_scaling(critical_temps_M, sizes,
                                                 source='Magnetisation')
    print(temp_inf_m)
    print((t_c - temp_inf_m[1]/temp_inf_m[1]), ' Standard errors away. '
                                             'Magnetisation')


def fss_specific_heat():
    critical_temps = investigate_heat_capacity(lattice_sizes=sizes,
                                               take_last=100000, re_runs=1)


    temp_inf = investigate_finite_size_scaling(critical_temps, sizes)
    print(temp_inf)

    print((t_c - temp_inf[0])/temp_inf[1], ' Standard errors away.  '
                                           'capacity')


# -------------------------------------------------------------------------
def external_field():
    print(investigate_chi(temperatures, sizes[::3]))
    investigate_magnetic_fluctuations(hs=ext_fields, take_last=100000,
                                  re_runs=1)


if __name__ == '__main':
    sanity_checks()
    fss_magnetisation()
    fss_specific_heat()
    external_field()
