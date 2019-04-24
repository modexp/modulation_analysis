# these functions are used in the all_slows_plot

def rebin(xlist, incr):
    '''Rebin list: xlist with incr number of items per bin'''
    xs = []
    for i in range(0, len(xlist), incr):
        xmax = min(len(xlist), i + incr)
        xs.append(np.average(xlist[i:xmax]))
    return np.array(xs)

def rebin_errors(dxlist, incr):
    '''Rebin list of errors: dxlist with incr number of items per bin'''
    dx = []
    for i in range(0, len(xlist), incr):
        plot_dysum = 0
        for sig in dxlist[i:xmax]: plot_dysum += (sig ** 2)
        dx.append(np.sqrt(plot_dysum) / len(plot_dy[i:xmax]))
    return np.array(dx)

def load_ana(anasource, incr = 1, val_range = [-1e99, 1e99], sigmacut = 3, ch = 0, pk = 0,
             save = False, startdate = False, limits = False):
    '''Returns plots of a property of anasource of the channels and all peaks'''
    if not anasource in in_ana: raise NameError('Not a right ana property')

    if startdate:
        t0 = startdate[0].timestamp() + toffset
        t1 = startdate[1].timestamp() + toffset
        cuts = [cutter(ch_0 = ch, pk_0 = pk, t_start = t0, t_end = t1)]
    else: cuts = [cutter(ch_0 = ch, pk_0 = pk)]

    this_time = time[cuts]
    ys = globals()[anasource][cuts]

    std_ys  = np.std(ys[np.all([ys<val_range[1], ys>val_range[0]],axis = 0)])
    mean_ys = np.mean(ys[np.all([ys<val_range[1], ys>val_range[0]],axis = 0)])
    numcut       = np.any([ys < mean_ys - sigmacut * std_ys, ys > mean_ys+ sigmacut * std_ys], axis = 0)
    print(anasource, 'cutted', sum(numcut),'outliers',np.mean(ys[~numcut]), np.std(ys[~numcut]))

    time_datetime = np.array([datetime.datetime.fromtimestamp(np.int(tim-toffset)
        ) for tim in rebin(this_time[~numcut], incr)])
    return time_datetime, rebin(ys[~numcut],incr), np.mean(ys[~numcut]), np.std(ys[~numcut])

def do_ax(ax, x, y, mu, sigma, name, dec):
    '''Make an ax for the plot of slow variables over time'''
    ax.axhline(mu+sigma, c= 'r', lw = 2)
    ax.axhline(mu-sigma, c= 'r', lw = 2)
    ax.axhline(mu,  c = 'green', lw = 2)
    ax.set_ylabel(name)

    min_max = max(y)-min(y)
    min_max= 0.5 * round(min_max, dec)
    mu = round(mu,dec)
    ax.set_ylim(min(y) - 0.2*sigma, max(y)+ .2*sigma)

    ax.scatter(x, y, s= 10, edgecolor='None')
    ax.set_ylim(ax.get_yticks()[0], ax.get_yticks()[-1])
    ax.set_yticks(ax.get_yticks()[1:-1])


    # some functions and some fits, some not used anymore

    def exp_function(x, tau, a):
        '''Returns a*2^x/tau with x as a timestamp (in seconds) and tau in years'''
        f = np.power(2, -((x-t0[0])*s2y)/tau)*a
        return f

    def fitter(chan, pk =0, chi_cut = 0):
        '''Fits an exponential for the given channel and the given peak, returns the optimized values of the exponent.
        If chi_cut > 0 it cuts away the rates that come from a spectrum that is fitted with a chi2/ndf larger than
        chi_cut times the sdt from the mean'''

        # If there is a cut on the chi2 of the fits exclude the fits with high chi2s.
        if chi_cut > 0 :
            chimax = np.mean(chi2ndf[np.all([channel == chan, peak ==pk],axis=0)]) + \
            chi_cut * np.std(chi2ndf[np.all([channel == chan, peak ==pk],axis=0)])
            cuts = np.all([channel == chan, chi2ndf < chimax, peak == pk], axis=0)

        else:
            cuts = np.all([channel == chan, peak == pk], axis=0)

        xs = time[cuts]
        ys = rate[cuts]
        yerror = drate[cuts]

        # Best guess for a and tau in the exp_funtion are a = ratemean and tau = halflife
        popt, pcov = curve_fit(exp_function, xs, ys, (halflife[chan],
                np.mean(rate[np.all([channel == chan, peak == pk],axis=0)])), sigma = yerror)

        return popt, pcov


    def xy_exp_fitter_no_tau(xs, ys, yerror, chan):
        '''Fits an exponential for the given the rate (ys and yerror) over time (xs) and the channel the optized
        values of the exponent.'''
        # Best guess for a and tau in the exp_funtion are a = ratemean and tau = halflife
        tau = halflife[chan]
        def exp_function_chan(x, a):
            '''Returns a*2^x/tau with x as a timestamp (in seconds) and tau in years'''
            f = np.power(2, -((x-t0[0])*s2y)/tau)*a
            return f
        popt, pcov = curve_fit(exp_function_chan, xs, ys, (np.mean(ys[0:10])) ,
            sigma = yerror)
        return popt, pcov

    def xy_exp_fitter(xs, ys, yerror, chan):
        '''Fits an exponential for the given the rate (ys and yerror) over time (xs) and the channel the optized
        values of the exponent.'''
        # Best guess for a and tau in the exp_funtion are a = ratemean and tau = halflife
        popt, pcov = curve_fit(exp_function, xs, ys, (halflife[chan], np.mean(ys[0:10])) , sigma = yerror)
        return popt, pcov

    def sin_function(x, A, phi):
        '''returns A sin(2pi (x-phi)) with a period of 1 in the same units as x'''
        f = A * np.sin(2 * np.pi * (x) - phi / 1.)
        return f

    def sinfitter(xs, ys,yerr):
        '''Fits a sin to the xs, ys and yerror of the data provided. Assumes that the amplitude is ~0.1 and
        the phase ~0.'''
        popt, pcov = curve_fit(sin_function, xs, ys, sigma = yerr, p0=(0.10, 0), maxfev = 8000)
        return popt, pcov

    def gaussian(x, mu, sigma):
        return np.exp((-(x-mu)**2)/(2 * sigma ** 2)) / (np.sqrt(2*np.pi) * sigma)

    def poisson(x, mu):
        return np.power(mu, x) * np.exp(-mu) / np.math.factorial(x)

    def lnpoisson(x, mu):
        prod = 0
        for j in x:
            prod += np.log(np.power(mu, int(j)) * np.exp(-mu) / np.math.factorial(int(j)))
        return prod


    print("Rateplotter::\tdefined functions exp_function, fitter, xy_exp_fitter_no_tau, \
          xy_exp_fitter, sin_function, sinfitter, gaussian, poisson, lnpoisson")

          # ln prior used in lnprob, lnlike used in u
          def lnprior(theta, tau12, dtau12, A_max, fixed_a = False):
              '''This is a prior of the parameters. Theta are the parameters which can be:
              length == 1: The overall rate
              length == 2: The overall rate + the halflife
              length == 3: The overall rate + the halflife + the phase of the modulation
              length == 4: The overall rate + the halflife + the phase & amplitude of the modulation
              you can also fit for a specified value of the modulation amplitude by giving a fixed_a'''
              if len(theta) == 1:
                  A = theta
                  tau, a, phi = tau12, 0, 0
              if len(theta) == 2:
                  [A, tau] = theta
                  a, phi = 0, 0
              if len(theta) == 3:
                  [A, tau, phi] = theta
                  if fixed_a: a = fixed_a
                  else: a = 0
              if len(theta) == 4: [A, tau, phi, a] = theta
              if  0 < A < A_max and -1 < a < 1 and - np.pi / 2 < phi < np.pi / 2:
                  return - 0.5 * ((tau12 - tau) ** 2 * (1.0 / (dtau12 ** 2)) - np.log(1 / (dtau12 ** 2)))
              return -np.inf


          def lnlike(theta, x0, x1, y, yerr, tau12, dtau12, T, fixed_a = False):
              '''This is the likelihood of the data given the parameters. Theta are the parameters which can be:
              length == 1: The overall rate
              length == 2: The overall rate + the halflife
              length == 3: The overall rate + the halflife + the phase of the modulation
              length == 4: The overall rate + the halflife + the phase & amplitude of the modulation
              you can also fit for a specified value of the modulation amplitude by giving a fixed_a.
              The data consists of start times (x0), end times (x1), rate (y), uncertainty of the rate (yerr),
              the literature value of the half-life time (tau12 +/- dtau12) and the modulation period (usually 1 yr)'''
              if len(theta) == 1:
                  A = theta[0]
                  tau, a, phi = tau12, 0, 0
              if len(theta) == 2:
                  [A, tau] = theta
                  a, phi = 0, 0
              if len(theta) == 3:
                  [A, tau, phi] = theta
                  if fixed_a: a = fixed_a
                  else: a = 0
              if len(theta) == 4: [A, tau, phi, a] = theta

              # Terms without modulation
              model = A * tau * (np.power(2, - x0 / tau) - np.power(2, - x1 / tau)) / np.log(2)

              # Terms with modulation integrated over time between x0 and x1
              if len(theta) > 2:
                  argx0 = 2 * np.pi * x0 / T + phi
                  argx1 = 2 * np.pi * x1 / T + phi
                  model+= (np.power(2, - (x0-x1) / tau) * a * A * T * tau * (
                      np.power(2, x1 / tau) * ( T * np.cos(argx0) * np.log(2) -
                      2 * np.pi * tau * np.sin(argx0)) +
                      np.power(2, x0 / tau) * (- T * np.cos(argx1) * np.log(2) +
                      2 * np.pi * tau * np.sin(argx1)))/((2 * np.pi * tau) ** 2 + (T * np.log(2)) ** 2))

              inv_sigma2 = 1.0 / ((yerr * (x1-x0)) ** 2)

              return np.sum(-0.5 * (( y * (x1 - x0) - model) ** 2 * inv_sigma2) - np.log(np.sqrt(2 * np.pi * inv_sigma2)))

          # Define the probability function as likelihood * prior.
          def lnprob(theta, x0, x1, y, yerr, tau12, dtau12, A_max, T, fixed_a = False):
              lp = lnprior(theta,tau12, dtau12, A_max, fixed_a)
              if not np.isfinite(lp):
                  return -np.inf
              return lp + lnlike(theta, x0, x1, y, yerr, tau12, dtau12, T, fixed_a)
