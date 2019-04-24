class data(object):
    "A class to store the data"
    def __init__( self ) :
        self.channel = [-1]
        self.peak    = -1
        self.time    = [toffset]
        self.t0      = [toffset]
        self.t1      = [toffset]
        self.t0_yr   = -1
        self.t1_yr   = -1
        self.rate    = -1
        self.drate   = -1
        self.r_fit   = -1
        self.r_mod   = -1
        self.r_sub   = -1
        self.dr_sub  = -1
        self.temp    = -1
        self.fit     = -1
        self.mod_fit = -1
        self.rate_fit= -1
        self.tstart  = 0
        self.tend    = 1e10
        self.axes    = False
        self.fig     = None



    def __str__(self):
        s = "Data"
        s+= "\n Source %s \tchannel %i - %i"%(sourceName[self.channel[0]],self.channel[0], self.channel[-1])
        s+= "\n peak en. %.f keV \tnumber %i"%(sourceE[self.channel[0]][self.peak],self.peak)
        s+= "\n between times: %s to %s" %(
            str(datetime.datetime.fromtimestamp(self.t0[0] - toffset)),
            str(datetime.datetime.fromtimestamp(self.t1[-1]  - toffset)))
        return s

    def __repr__(self):
        return self.__str__()

    def load(self):
        print("\tdata::\tLoad channel %i and peak %i"%(self.channel[0], self.peak))
        cuts = [cutter(ch_0 = ch, pk_0 = self.peak, t_start = self.tstart, t_end = self.tend) for ch in self.channel]

        self.rate    = np.sum(rate[cuts[i]] for i in range(len(self.channel)))
        self.drate   = np.sqrt(np.sum(drate[cuts[i]] ** 2 for i in range(len(self.channel))))

        self.bgrate  = np.sum(bgrate[cuts[i]] for i in range(len(self.channel)))
        self.bgdrate = np.sqrt(np.sum(bgdrate[cuts[i]] ** 2 for i in range(len(self.channel))))

        self.temp    = temp[cuts[0]]
        self.time    = time[cuts[0]]
        self.t0      = time_start[cuts[0]]
        self.t1      = time_end[cuts[0]]
        self.t0_yr   = (self.t0 - self.t0[0]) * s2y
        self.t1_yr   = (self.t1 - self.t0[0]) * s2y


    def fit_exp(self, nwalkers = 50, nsteps = 50, burnin = 10, plot_level = []):
        '''Fit exponential decay'''
        print("\tdata::\tFit channel %i and peak %i with an exponential decay"%(self.channel[0], self.peak))
        if type(self.rate) == int: data.load(self)
        if self.channel[0] == 4: dt12 = 10 * dhalflife[self.channel[0]]
        else: dt12 = dhalflife[self.channel[0]]
        self.fit = use_minuit(self.t0_yr, self.t1_yr, self.rate, self.drate,
            halflife[self.channel[0]], dt12, self.channel[0], self.peak,
            ndim = 2, nwalkers = nwalkers, nsteps = nsteps, burnin = burnin,
            plot_level = plot_level, save = False)
        res = np.transpose(self.fit)
        A0, t12, phi, a = res[0]
        dA0, dt12, dphi, da = (res[1] + res[2]) / 2

        t_year = (self.time - self.t0[0]) * s2y

        self.r_fit  = np.power(2, -(t_year / t12)) * A0

    # This is not recommended
    def fit_mod(self, nwalkers = 100, nsteps = 150, burnin = 50, plot_level = []):
        '''Fit exponential decay with a modulation using MCMC'''
        print("\tdata::\tFit channel %i and peak %i with an exponential decay with modulation"%(
                self.channel[0], self.peak))
        if type(self.rate) == int: data.load(self)
        if 4 in self.channel: dt12 = 10 *  dhalflife[self.channel[0]]
        else: dt12 = dhalflife[self.channel[0]]

        self.mod_fit = fit_A_tau_a_phi(self.t0_yr, self.t1_yr, self.rate, self.drate,
                halflife[self.channel[0]], dt12, self.channel[0], self.peak,
                ndim = 4, nwalkers = nwalkers, nsteps = nsteps, burnin = burnin,
                plot_level = plot_level, save = False)
        res = np.transpose(self.mod_fit)
        A0, t12, phi, a = res[0]
        dA0, dt12, dphi, da = (res[1] + res[2]) / 2

        t_year = (self.time - self.t0[0]) * s2y

        self.r_fit = np.power(2, -(t_year / t12)) * A0
        self.r_mod = a * np.sin((2 * np.pi * (t_year) / 1 + phi))

    # This is recommended
    def fit_mod_minuit(self, plot_level = []):
        '''Fit exponential decay with a modulation using Minuit'''
        print("\tdata::\tFit channel %i and peak %i with an exponential decay with modulation"%(
                self.channel[0], self.peak))
        if type(self.rate) == int: data.load(self)
        if 4 in self.channel: dt12 = dhalflife[self.channel[0]]
        else: dt12 = dhalflife[self.channel[0]]

        self.mod_fit = use_minuit(self.t0_yr, self.t1_yr, self.rate, self.drate,
                halflife[self.channel[0]], dt12, self.channel[0], self.peak,
                ndim = 4, plot_level = plot_level, save = False)
        res = np.transpose(self.mod_fit)
        A0, t12, phi, a = res[0]
        dA0, dt12, dphi, da = (res[1] + res[2]) / 2

        t_year = (self.time - self.t0[0]) * s2y

        self.r_fit = np.power(2, -(t_year / t12)) * A0
        self.r_mod = a * np.sin((2 * np.pi * (t_year) / 1 + phi))

    def cal_ln(self, fit_a):
        '''Calculate the log-likelihood'''
        if type(self.rate) == int: data.load(self)

        return use_minuit(self.t0_yr, self.t1_yr, self.rate, self.drate,
                halflife[self.channel[0]], dhalflife[self.channel[0]], self.channel[0], self.peak,
                ndim = 3, fit_a = fit_a)


    def subtract_exp(self):
        '''Subtract the fitted exponential decay from the data'''
        print("\tdata::\tSubtract channel %i and peak %i with an exponential decay with modulation"%(
                self.channel[0], self.peak))
        if self.mod_fit != -1: res = np.transpose(self.mod_fit)
        elif self.fit != -1: res = np.transpose(self.fit)
        else: self.fit_exp(); res = np.transpose(self.fit)

        A0, t12, phi, a = res[0]
        dA0, dt12, dphi, da = (res[1] + res[2]) / 2

        t_year = (self.time - self.t0[0]) * s2y

        self.r_sub = 100 * ((self.rate / self.r_fit) - 1)

        self.dr_sub = np.sqrt(
            # Error due to error on the rate
            625 * np.power(2, 4 + 2 * t_year / t12) * ((self.drate / A0) ** 2) +
            # Due to error in A0
            625 * np.power(2, 4 + 2 * t_year / t12) * (dA0 ** 2) * (self.rate ** 2) / (A0 ** 4) +
            # Due to error in t12
            625 * np.power(2, 4 + 2 * t_year / t12) * (dt12 ** 2) * (t_year ** 2) * (self.rate ** 2)  *
            (np.log(2) ** 2) / ((A0 ** 2)*( t12 ** 4)) )

    def correct_temp(self):
        '''Do a temperature correction'''
        print("Correct the data for fluctuations in the data")
        if type(self.rate) == int: data.load(self)
        Tset = 30 #np.mean(self.temp)
        print("Before correction rmean",np.mean(self.rate),
              "tmean, tmax, tmin, lenT", np.mean(self.temp), np.max(self.temp), np.min(self.temp), len(self.temp) )
        self.temp = np.clip(self.temp, 10, 40)

        self.rate, self.drate = T_correct(self.rate, self.drate, self.temp, Tset, self.channel, pk = self.peak)
        self.temp = np.clip(self.temp, Tset, Tset)
        print("After correction rmean",np.mean(self.rate),
              "tmean, tmax, tmin, lenT", np.mean(self.temp), np.max(self.temp), np.min(self.temp), len(self.temp) )

    def rebin(self, incr):
        '''Rebin the data for plotting purposes'''
        xs, ys, err= [], [], []
        for i in range(0, len(self.time), incr):
            xmax = min(len(self.time), i + incr)
            xs.append(np.average(self.time[i:xmax]))
            ys.append(np.average(self.rate[i:xmax]))


            # For the error we sum the errors quadratically and devide by incr
            errsum = 0
            for sig in self.drate[i:xmax]: errsum += (sig ** 2)
            errsum  = np.sqrt(errsum) / len(self.drate[i:xmax])
            err.append(errsum)
        self.time, self.rate, self.drate = np.array(xs), np.array(ys), np.array(err)

    def show(self, sub = False, save= False, incr = 1, savename = 'test', axnum = False):
        '''Show the data'''
        def option(value, num):
            if type(num)== bool and not num: return value
            elif num == 0: return value
            else: return False
        if type(axnum) == list:
            try: ax = (self.axes)[axnum[0]]
            except TypeError:
                self.fig, (self.axes) = plt.subplots(axnum[1], sharex=True, sharey=False)
                ax = (self.axes)[axnum[0]]
                self.fig.subplots_adjust(hspace=0.0)
                self.fig.autofmt_xdate()
                plt.setp([a.get_xticklabels() for a in self.fig.axes[:-1]], visible=False)

        else: self.fig, ax = plt.subplots()
        if sub == False:
            rate, drate, r_fit, fit = self.rate, self.drate, self.r_fit, self.fit
        else:
            rate, drate, r_fit, fit = self.r_sub, self.dr_sub, self.r_mod, self.fit
        if sub == False and self.mod_fit != -1:
            rate, drate, r_fit, fit = self.rate, self.drate, (1 + self.r_mod) * self.r_fit, self.mod_fit
        if sub and self.mod_fit != -1: rate, drate, r_fit, fit = self.r_sub, self.dr_sub, self.r_mod, self.mod_fit;
        if self.mod_fit == -1 and not sub: rate, drate, r_fit, fit = self.rate, self.drate, self.r_fit, self.fit;
#        print(rate, drate, r_fit, fit )
        binned_rate_xy(self.fig,
                       ax,
                       self.channel,
                       self.peak,
                       self.time,
                       rate,
                       drate,
                       r_fit,
                       fit,
                       coloritem = self.temp,
                       lit_val = False,
                       exp_text = option(True, 0),
                       legend = False, #True if axnum[0] == 2 else False,
                       fitted_pars = 4 if self.mod_fit !=-1 else 2,
                       sub = sub,
                       save = False,
                       savename = savename,
                       jenkins = True,
                       doplot = False,
                       incr = incr)
        if save and ((axnum == False and type(axnum)==bool) or ax):
            self.fig.savefig(droppath + savename+'.pdf', dpi = 300,  bbox_inches='tight')
            self.fig.savefig(droppath + savename+'.png', dpi = 300,  bbox_inches='tight')

        plt.rcParams['figure.figsize'] = (12.0, 12.0)
