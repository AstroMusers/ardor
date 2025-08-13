"""
Created on Wed Dec 18 15:29:39 2024

@author: whitsett.n
"""

import ardor.SPI_Forward_Models.SPI_Simulation as SPI
import ardor.SPI_Forward_Models.Orbit_Model_Library as OML
import ardor.Flares.Flare as Flare
import ardor.Utils.Utils as U
import ardor.Statistical_Tests.K_Tests as K
import ardor.Plotting.Orbital_Flare_Hist as Plot
from scipy.interpolate import interp1d
from scipy.optimize import basinhopping, minimize
from scipy.stats import vonmises, uniform, norm
from scipy.integrate import simpson
from matplotlib import pyplot as plt
from astropy.stats import kuiper
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import pandas as pd
import math
def round_to_sf(x, sf):
    """
    Rounds a number to a specified number of significant figures.

    Args:
        x: The number to round.
        sf: The number of significant figures.

    Returns:
        The rounded number.
    """
    if x == 0:
        return 0.0
    
    sign = 1 if x > 0 else -1
    x = abs(x)
    
    rounding_position = sf - int(math.floor(math.log10(x))) - 1
    
    rounded_x = round(x, rounding_position)
    
    return rounded_x * sign
def Inverse_Cubic_Function(ratio, loc, e, a, period):
    model, phase = SPI.SPI_Cubic(ratio, loc, e, a, period, length=10)
    f = interp1d(np.linspace(0, 2*np.pi, num=len(model)) ,model, kind='linear')
    return f
def VM_pdf(x, loc, kappa):
    return vonmises.pdf(x, kappa= kappa, loc=loc)

def VM_pdf_norm(x, loc, kappa):
    x = U.range_shift(x, 0, 1, -np.pi, np.pi)
    loc_norm = U.range_shift(loc, 0, 1, -np.pi, np.pi)
    pdf = vonmises.pdf(x, kappa= kappa, loc=loc_norm)
    pdf_norm = U.range_shift(pdf, -np.pi, np.pi, 0, 1)
    return pdf_norm

def cubic_pdf(x, ratio, loc, e, a, period):
    return Inverse_Cubic_Function(ratio, loc, e, a, period)(x)
def VM_likelihood(params, data):
    loc, kappa = params

    return -np.sum(np.log(VM_pdf(data, loc, kappa)))
def cubic_likelihood(params,data, e, a,period):
    loc, ratio = params
    return -np.sum(np.log(cubic_pdf(data, ratio, loc, e, a, period)))
def null_likelihood(data):
    return -np.sum(np.log(uniform.pdf(data, loc = -np.pi, scale=2*np.pi)))
def null_likelihood_cubic(data):
    return -np.sum(np.log(uniform.pdf(data, loc = 0, scale=2*np.pi)))
def VM_Unbinned_likelihood(flares, rad = False):
    if np.isnan(np.mean(flares)) == False and len(flares) >= 3:
        if rad == False:
            flares = U.range_shift(flares, 0, 1, -np.pi, np.pi)
            # results = basinhopping(VM_likelihood, x0 = [0, 1], minimizer_kwargs={'args': flares,'bounds': [(-np.pi, np.pi), (0, None)]})
            results = minimize(VM_likelihood, ([0, 2]), args = (flares), bounds = [(-np.pi, np.pi), (0, None)])
            null = null_likelihood(flares)
            location = round_to_sf(U.range_shift(results.x[0], -np.pi, np.pi, 0, 1), 3)
            kappa = round_to_sf(results.x[1], 3)
            TS = 2*(null - results.fun)
            return location, kappa, np.sqrt(TS)
        elif rad == True:
            results = minimize(VM_likelihood, ([0, 2]), args = (flares), bounds = [(-np.pi, np.pi), (0, None)])
            null = null_likelihood(flares)
            location = round_to_sf(results.x[0],3)
            kappa = round_to_sf(results.x[1], 3)
            TS = 2*(null - results.fun)
            return location, kappa, np.sqrt(TS)
    else:
        return None
def Cubic_Unbinned_likelihood(flares,e,a,period):
    flares = U.range_shift(flares, 0, 1, 0, 2*np.pi)
    results = minimize(cubic_likelihood, ([0, 0]), args = (flares,e, a, period), bounds = [(0, 2*np.pi), (0,1e4)], tol=1e-8)
    null = null_likelihood_cubic(flares)
    loc = results.x[0]
    ratio = results.x[1]
    TS =2*(null - results.fun)
    return loc, ratio, TS
def VM_Unbinned_likelihood_List(DataFrame, output_dir, phase_header="Transit_Phase"):
    hosts = list(set(DataFrame["Host_ID"]))
    host_list = []
    kappa = []
    TS = []
    loc = []
    new_ID = []
    for host in hosts:
        flares = DataFrame.loc[DataFrame["Host_ID"] == host, phase_header]
        flares = U.range_shift(flares, 0, 1, -np.pi, np.pi)
        if np.isnan(np.mean(flares)) == False and len(flares) > 3:
            results = minimize(VM_likelihood, ([0, 2]), args = (flares), bounds = [(-np.pi, np.pi), (0, None)])
            null = null_likelihood(flares)
            if null == results.fun:
                print(1)
            host_list.append(host)
            loc.append(round_to_sf(U.range_shift(results.x[0], -np.pi, np.pi, 0, 1), 3))
            kappa.append(round_to_sf(results.x[1], 3))
            TS.append(round_to_sf(2*(null - results.fun), 2))
            new_ID.append(host)
    new_data = pd.DataFrame({"Host_ID": new_ID, "Kappa": kappa, "Center": loc, "TS_{VM}":TS})
    new_data.to_csv(output_dir, index=False)
def Cubic_Unbinned_likelihood_List(DataFrame, output_dir, e, a,period, TOI = False):
    hosts = list(set(DataFrame["Host_ID"]))
    host_list = []
    ratio = []
    loc = []
    new_ID = []
    TS= []
    for host in hosts:
        flares = DataFrame.loc[DataFrame["Host_ID"] == host, "Periastron_Phase"]
        results = minimize(cubic_likelihood, ([10, 5]), args = (flares,e, a, period), bounds = [(1e-10, 1e8), (0,2*np.pi)], tol=1e-8)
        null = null_likelihood_cubic(flares)
        host_list.append(host)
        loc.append(results.x[1])
        ratio.append(results.x[0])
        TS.append(round_to_sf(2*(null - results.fun), 2))
        new_ID.append(host[-3:])
    new_data = pd.DataFrame({"Host_ID": new_ID, "ratio": ratio, "Center": loc, "TS_{VM}":TS})
    new_data.to_csv(output_dir, index=False)

def UBL_Sim_Kappa(DataFrame, output_dir, star=''):
    # interest_hosts = [['p_KU', 'p_KS', 'p_AD', 'p_KS_Samp', 'p_AD_Samp', 'N']]
    data_frame = pd.DataFrame()
    set_hosts = set(DataFrame['Host_ID'])
    kappa0_kappa = []
    kappa25_kappa = []
    kappa5_kappa = []
    kappa1_kappa = []
    kappa2_kappa = []
    kappa4_kappa = []
    kappa8_kappa = []
    loc0_kappa = []
    loc25_kappa = []
    loc5_kappa = []
    loc1_kappa = []
    loc2_kappa = []
    loc4_kappa = []
    loc8_kappa = []
    TS0_kappa = []
    TS25_kappa = []
    TS5_kappa = []
    TS1_kappa = []
    TS2_kappa = []
    TS4_kappa = []
    TS8_kappa = []
    kappa_list = [kappa0_kappa,
        kappa25_kappa,
        kappa5_kappa,
        kappa1_kappa,
        kappa2_kappa,
        kappa4_kappa,
        kappa8_kappa]
    loc_list = [loc0_kappa,
        loc25_kappa,
        loc5_kappa,
        loc1_kappa,
        loc2_kappa,
        loc4_kappa,
        loc8_kappa]
    TS_list = [TS0_kappa,
        TS25_kappa,
        TS5_kappa,
        TS1_kappa,
        TS2_kappa,
        TS4_kappa,
        TS8_kappa]
    for index, hosts in enumerate(set_hosts):

        if float(hosts[-1]) == 8 or float(hosts[-1]) == 4 or float(hosts[-1]) == 0 or float(hosts[-1]) == 2 or float(hosts[-1]) == 1:
            kappa = float(hosts[-1])
        elif float(hosts[-3:]) == 0.5:
            kappa = 0.5
        elif float(hosts[-4:]) == 0.25:
            kappa = 0.25
        peri_phases = np.array(DataFrame.loc[DataFrame['Host_ID'] == hosts, 'Phase'])
        peri_phases = U.range_shift(peri_phases, 0, 1, -np.pi, np.pi)
        results = VM_Unbinned_likelihood(peri_phases)
        loc = U.range_shift(results.x[0],-np.pi, np.pi, 0, 1)
        kappas= results.x[1]
        null = null_likelihood(peri_phases)
        TS = 2*(null - results.fun)
        if kappa == 0:
            kappa0_kappa.append(kappas)
            loc0_kappa.append(loc)
            TS0_kappa.append(TS)
        elif kappa == 0.25:
            kappa25_kappa.append(kappas)
            loc25_kappa.append(loc)
            TS25_kappa.append(TS)
        elif kappa == 0.5:
            kappa5_kappa.append(kappas)
            loc5_kappa.append(loc)
            TS5_kappa.append(TS)
        elif kappa == 1:
            kappa1_kappa.append(kappas)
            loc1_kappa.append(loc)
            TS1_kappa.append(TS)
        elif kappa == 2:
            kappa2_kappa.append(kappas)
            loc2_kappa.append(loc)
            TS2_kappa.append(TS)
        elif kappa == 4:
            kappa4_kappa.append(kappas)
            loc4_kappa.append(loc)
            TS4_kappa.append(TS)
        elif kappa == 8:
            kappa8_kappa.append(kappas)
            loc8_kappa.append(loc)
            TS8_kappa.append(TS)
    for index, kappas in enumerate([0,0.25,0.5,1,2,4,8]):
        data_frame[str(star) + "Kappa_" + str(kappas)] = kappa_list[index]
        data_frame[str(star) + "Center_" + str(kappas)] = loc_list[index]
        data_frame[str(star) + "TS_{VM}_" + str(kappas)] = TS_list[index]
    data_frame.to_csv(output_dir, index = False)

def SPI_Metric(period, e, a, phases, arg_periastron, Rst, t_tran = 0, transit = False, los_factor = 1):
    '''
    Returns an SPI metric for the provided orbit and PDF of the SPI model for
    induced flares.

    Parameters
    ----------
    period: float
        Period of the orbit, in daysa : float
        Semi-major axis of the orbit, in AU
    e : float
        Eccentricity of the orbit
    a : float
        Semi-major axis of orbit, in AU
    PDF : array-like
        Probability density function of the flare model, in units of normalized
        phase [0,1], with periastron at 0
    arg_periastron : float
        Argument of periastron, in degrees
    Rst : float
        Radius of the star, in solar radii
    los_factor : float
        Factor in which you want to blind the metric to flares occuring when
        the interaction is presumably out of line-of-sight
    Returns
    -------
    log(SPI_Metric): large positive values indicate evidence for SPI

    '''
    Rst = Rst*0.00465047
    lower,upper = OML.secondary_eclipse_true_anomaly_range_exact(a, e, 90, arg_periastron, Rst, period, t_tran)
    rad_phases= U.range_shift(phases, 0,1,-np.pi,np.pi)
    rad_loc, kappa, TS = VM_Unbinned_likelihood(rad_phases, rad=True)

    PDF = vonmises.pdf(np.linspace(-np.pi,np.pi,num=2000), kappa = kappa, loc = rad_loc)
    
    PDF_phase = U.range_shift(np.linspace(-np.pi,np.pi,num=2000), -np.pi,np.pi,0,1)
    phase, dist, rot = OML.orbit_pos_v_time(period, e, a, len(PDF), phase = True, 
                                 arg_periastron=0)
    # print("Results of KU Test on Fitted Dist: ", KS[1])
    # print("Results of the KU test against the uniform CDF: ", KU_test[1])
    # print("Results of the TS for the UBL fit: ", np.sqrt(TS))
    if len(dist) < len(PDF):
        PDF = np.delete(PDF,0)
        PDF_phase = np.delete(PDF_phase,0)
    PDF = PDF/simpson(PDF,x=PDF_phase)
    
    norm_orbit = (dist)**-3
    if transit == True:
        print(lower, upper)
        for idx, values in enumerate(phase):
            a = U.boudnary_conditions(((lower - los_factor*(upper-lower)/2)/(2*np.pi)), 0, 1)
            b = U.boudnary_conditions(((upper+los_factor*(upper-lower)/2)/(2*np.pi)), 0, 1)
            if a < b:
                if values > a and values < b:
                    norm_orbit[idx] = 0
            elif b < a:
                if values < b or values > a:
                    norm_orbit[idx] = 0
    # plt.axvline(Flare.phase_folder(t_tran, period, 0))
    normalized = norm_orbit/simpson(norm_orbit, x=phase)
    fig, ax = plt.subplots(figsize=(4,4))
    ax.plot(phase, normalized*PDF, label = r'$P(r_{{norm}}|VM)$')
    ax.plot(phase, normalized, label='$r_{norm}$')
    ax.plot(phase, PDF, label = 'VM PDF')
    ax.set_xlabel("Norm. Phase")
    ax.set_ylabel("Density")
    ax.set_ylim(0.0, 6)
    ax.vlines(U.range_shift(rad_loc, -np.pi, np.pi, 0, 1), ymin=0, ymax=10, alpha=0.6, linestyle='--', label=r'$\kappa$')
    # plt.hist(U.range_shift(rad_phases, -np.pi,np.pi,0,1), density=True, bins=15, label = 'Flare Density')
    ax.legend(loc='upper right')
    ones = np.ones(len(PDF))
    SPI_Metric = simpson((PDF * (normalized)), x = phase)
    b = simpson((PDF * (ones)), x = phase)
    beta = SPI_Metric - b
    ax.text(0.03, 5.5, fr"$\beta_{{SPI}}=${beta:.3f}", fontsize=14)
    plt.show()
    return SPI_Metric - b

def create_animated_plot(functions, variables, period, e, a, phases, omega, Rst, gif_name='animation.gif'):
    """
    Create an animated GIF of multiple functions updating over a variable.

    Args:
        functions (list of callables): Each function must accept (x, variable) and return y.
        variables (list or array): Values to animate over.
        gif_name (str): Output file name.
    """
    x = np.linspace(0, 2 * np.pi, 200)

    fig, ax = plt.subplots()
    lines = []

    # Create one Line2D for each function
    for func in functions:
        line, = ax.plot([], [], lw=2)
        lines.append(line)

    ax.set_xlim(0, 2 * np.pi)
    ax.set_ylim(-2, 2)
    ax.set_title('Multi-Function Animation')

    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    def update(variable):
        for func, line in zip(functions, lines):
            y = func(x, variable)
            line.set_data(x, y)
        ax.set_title(f'Variable = {variable:.2f}')
        return lines

    anim = FuncAnimation(
        fig, update, frames=variables, init_func=init,
        blit=True, repeat=False
    )

    anim.save(gif_name, writer=PillowWriter(fps=10))
    print(f"Saved animation to {gif_name}")

    plt.close(fig)