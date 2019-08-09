'''Make the "lobster plot."

That is, the possible values of the Majorana mass given the mixing parameters
and assuming standard three-neutrino mixing, for the normal and inverted
mass hierarchy.

By the way, the current global best fit parameters from NuFit are:

    sin(t12)^2 = 0.304 +/- 0.012
    sin(t23)^2 = 0.451 +/- 0.002 (first quadrant)
    sin(t23)^2 = 0.576 +0.024 -0.037 (second quadrant)
    sin(t13)^2 = 0.021 +0.593 -0.560
        dm21sq = 7.53e-5 +/- 0.18e-5
        dm32sq = 2.44e-3 +/- 0.06e-3

.. todo:: Overlay bands to show the uncertainties on the mixing parameters.

.. moduleauthor:: Andy Mastbaum <mastbaum@hep.upenn.edu>, May 2015

'''

import numpy as np
#import ROOT
import matplotlib.pyplot as plt

def mbb(m1, s12sq, s13sq, s23sq, dm21sq, dm32sq, phi2, phi3, hierarchy):
    '''Compute the effective Majorana mass for given parameters.

    For light neutrino exchange-mediated neutrinoless double-beta decay, the
    effective Majorana mass is

        mbb = sum [ (U_ek)^2 m_k ], k = 1, 2, 3

    :param m1: Mass of the m1 eigenstate
    :param s12sq: sin^2 theta12
    :param s13sq: sin^2 theta13
    :param s23sq: sin^2 theta23
    :param dm21sq: Delta m21 squared
    :param dm32sq: Delta m32 squared
    :param phi2: First Majorana phase
    :param phi3: Second Majorana phase
    :param hierarchy: Positive for normal, negative for inverted
    :returns: mbb in the same units as m1
    '''
    t12 = np.arcsin(np.sqrt(s12sq))
    t13 = np.arcsin(np.sqrt(s13sq))
    t23 = np.arcsin(np.sqrt(s23sq))

    if hierarchy > 0:
        m2 = np.sqrt(np.square(m1) + dm21sq)
        m3 = np.sqrt(np.square(m2) + dm32sq)
    elif hierarchy < 0:
        # m3 is the lighest!
        m2 = np.sqrt(np.square(m1) + dm32sq)
        m3 = np.sqrt(np.square(m2) + dm21sq)
        m1, m3 = m3, m1

    mbb = np.abs(np.square(np.cos(t13)) * np.square(np.cos(t12)) * m1 +
                 np.square(np.cos(t13)) * s12sq * m2 * np.exp(1j * phi2) +
                 s13sq * m3 * np.exp(1j * phi3))

    return mbb


def mbb_range(m1, s12sq, s13sq, s23sq, dm21sq, dm32sq,
              hierarchy, step=0.001):
    '''Find a range of mbbs for given parameters by varying the unknown
    Majorana phases over their allowed range from 0 to 2pi.

    :param m1: Mass of the m1 eigenstate
    :param s12sq: sin^2 theta12
    :param s13sq: sin^2 theta13
    :param s23sq: sin^2 theta23
    :param dm21sq: Delta m21 squared
    :param dm32sq: Delta m32 squared
    :param hierarchy: Positive for normal, negative for inverted
    :returns: A tuple of (min, max) mbb in the same units as m1
    '''
    phis = np.arange(0, 2 * np.pi, step)

    # Create a grid of mbbs with phi2 and phi3 taking on all combinations
    mbbs = mbb(m1, s12sq, s13sq, s23sq, dm21sq, dm32sq,
               phis, phis[:,np.newaxis], hierarchy=hierarchy)

    return (np.min(mbbs), np.max(mbbs))

def mb(m1, s12sq, s13sq, s23sq, dm21sq, dm32sq, phi2, phi3, hierarchy):
    '''Compute the effective Majorana mass for given parameters.

    For light neutrino exchange-mediated neutrinoless double-beta decay, the
    effective Majorana mass is

        mbb = sum [ (U_ek)^2 m_k ], k = 1, 2, 3

    :param m1: Mass of the m1 eigenstate
    :param s12sq: sin^2 theta12
    :param s13sq: sin^2 theta13
    :param s23sq: sin^2 theta23
    :param dm21sq: Delta m21 squared
    :param dm32sq: Delta m32 squared
    :param phi2: First Majorana phase
    :param phi3: Second Majorana phase
    :param hierarchy: Positive for normal, negative for inverted
    :returns: mbb in the same units as m1
    '''
    t12 = np.arcsin(np.sqrt(s12sq))
    t13 = np.arcsin(np.sqrt(s13sq))
    t23 = np.arcsin(np.sqrt(s23sq))

    if hierarchy > 0:
        m2 = np.sqrt(np.square(m1) + dm21sq)
        m3 = np.sqrt(np.square(m2) + dm32sq)
    elif hierarchy < 0:
        # m3 is the lighest!
        m2 = np.sqrt(np.square(m1) + dm32sq)
        m3 = np.sqrt(np.square(m2) + dm21sq)
        m1, m3 = m3, m1

    mb = np.abs(np.square(np.cos(t13)) * np.square(np.cos(t12)) * m1 +
                 np.square(np.cos(t13)) * s12sq * m2 * np.exp(1j * phi2) +
                 s13sq * m3 * np.exp(1j * phi3))

    return mb
def plot(s12sq, s13sq, s23sq, dm21sq, dm32sq,
         hierarchy, step=0.005, name='graph'):
    '''Make a plot for the given parameters.

    :param s12sq: sin^2 theta12
    :param s13sq: sin^2 theta13
    :param s23sq: sin^2 theta23
    :param dm21sq: Delta m21 squared
    :param dm32sq: Delta m32 squared
    :param hierarchy: Positive for normal, negative for inverted
    :returns: A TGraph with the allowed area filled in
    '''
    def plot(s12sq, s13sq, s23sq, dm21sq, dm32sq,
         hierarchy, step=0.005, name='graph'):
    '''Make a plot for the given parameters.

    :param s12sq: sin^2 theta12
    :param s13sq: sin^2 theta13
    :param s23sq: sin^2 theta23
    :param dm21sq: Delta m21 squared
    :param dm32sq: Delta m32 squared
    :param hierarchy: Positive for normal, negative for inverted

    '''
    # Values of m_lightest to evaluate [eV]
    x = np.logspace(-5, 0, 100)

    # Loop over m_lightests to find the mbb range for each
    y = []
    for m in x:
        r = mbb_range(m, s12sq, s13sq, s23sq, dm21sq, dm32sq,
                      hierarchy=hierarchy, step=step)
        y.append(r)

    # Invert to get a tuple of (lower bounds, upper bounds)
    b = [list(t) for t in y]

    up = [y[i][0] for i in range(len(y))]
    dw = [y[i][1] for i in range(len(y))]

    #g = ROOT.TGraph(2 * len(x))
    #g.SetName(name)

    ind = np.arange(len(x))

    if hierarchy == -1: color = 'r'
    else: color = 'b'

    plt.yscale('log')
    plt.xscale('log')
    plt.fill_between(x, up, dw, color=color, alpha=0.5)

    plt.ylabel('|<$m_{\\beta \\beta}$>| (eV)', fontsize=14)
    plt.xlabel('$m_{lightest}$ (eV)', fontsize=14)

def plot_root(s12sq, s13sq, s23sq, dm21sq, dm32sq,
         hierarchy=1, step=0.005, name='graph'):
    '''Make a plot for the given parameters.

    :param s12sq: sin^2 theta12
    :param s13sq: sin^2 theta13
    :param s23sq: sin^2 theta23
    :param dm21sq: Delta m21 squared
    :param dm32sq: Delta m32 squared
    :param hierarchy: Positive for normal, negative for inverted
    :returns: A TGraph with the allowed area filled in
    '''
    # Values of m_lightest to evaluate [eV]
    x = np.logspace(-5, 0, 100)

    # Loop over m_lightests to find the mbb range for each
    y = []
    for m in x:
        r = mbb_range(m, s12sq, s13sq, s23sq, dm21sq, dm32sq,
                      hierarchy=hierarchy, step=step)
        y.append(r)

    # Invert to get a tuple of (lower bounds, upper bounds)
    y = np.array(zip(*y))

    #g = ROOT.TGraph(2 * len(x))
    #g.SetName(name)

    ind = np.arange(len(x))
    ind2 = np.arange(len(x), 2*len(x))
    x2 = np.logspace(0, -5, 100)
    #y2 = [i for i in y[1][]]

    plt.plot(ind, x, y[0])
    #plt.plot(ind2, x2, y2[1])

    '''for i in range(len(x)):
        g.SetPoint(i, x[i], y[0][i])
        g.SetPoint(len(x)+i, x[-i-1], y[1][-i-1])

    g.Draw('goff')

    return g'''


def lobsterROOT(s12sq, s13sq, s23sq, dm21sq, dm32sq):
    '''Make the whole lobster plot for the given parameters.

    :param s12sq: sin^2 theta12
    :param s13sq: sin^2 theta13
    :param s23sq: sin^2 theta23
    :param dm21sq: Delta m21 squared
    :param dm32sq: Delta m32 squared
    :returns: A tuple of (NH graph, IH graph, canvas)
    '''
    c = ROOT.TCanvas('c', 'c', 500, 500)
    c.SetLogx()
    c.SetLogy()

    # Normal hierarchy
    gn = plot(s12sq, s13sq, s23sq, dm21sq, dm32sq, 1, name='normal')
    gn.SetFillColor(ROOT.kRed - 7)

    # Inverted hierarchy
    gi = plot(s12sq, s13sq, s23sq, dm21sq, dm32sq, -1, name='inverted')
    gi.SetFillColor(ROOT.kAzure - 3)

    # Draw
    gi.Draw('af')
    gi.SetTitle('')
    gi.GetXaxis().SetTitle('Lightest #nu mass (eV)')
    gi.GetYaxis().SetTitle('|<m_{#beta#beta}>| (eV)')
    gi.GetXaxis().SetLimits(1e-4, 1)  # Oh, ROOT...
    gi.GetYaxis().SetRangeUser(1e-4, 1)

    gn.Draw('f same')

    ROOT.gPad.RedrawAxis()
    c.Update()
    return gn, gi, c

def lobster(s12sq, s13sq, s23sq, dm21sq, dm32sq):
    fig = plt.figure(figsize(10,10))
    plt.yscale('log')
    plt.xscale('log')

    #Normal hierarchy
    plot(s12sq, s13sq, s23sq, dm21sq, dm32sq, 1, name='normal')

    # Inverted hierarchy
    gi = plot(s12sq, s13sq, s23sq, dm21sq, dm32sq, -1, name='inverted')


if __name__ == '__main__':
    # Use the current best-fit parameters (and first-quadrant t23)
    s12sq = 0.304
    s23sq = 0.451
    s13sq = 0.021
    dm21sq = 7.53e-5
    dm32sq = 2.44e-3

    gn, gi, c = lobster(s12sq, s13sq, s23sq, dm21sq, dm32sq)

    raw_input()
