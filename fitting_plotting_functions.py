# I often make use of cuts therefore use:
def cutter(t_start = -1, t_end = -1, ch_0 = -1, ch_1 = -1, pk_0 = -1,
    pk_1 = -1, pk_2 = -1, T_low = -1, T_high = -1, P_min = -1, P_max = -1,
    chi_max = -1):
    '''Feed cuts to cut a ANA file and will return the right list. Note that we don't want -1's'''
    try: cuts = [t0>0]
    except: cuts = [t0[0]>0]
    if t_start   != -1:
        cuts = np.all([cuts, time > t_start], axis = 0)
    if t_end     != -1:
        cuts = np.all([cuts, time < t_end  ], axis = 0)
    if ch_1 < ch_0:
        ch_1, ch_0 = ch_0, ch_1
    if ch_1 != -1:
        if ch_0 != -1:
            cuts = np.all([cuts, np.any([channel == ch_0, channel == ch_1],axis = 0)], axis = 0)
        else:
            cuts = np.all([cuts, channel == ch_1], axis = 0)
    if pk_0    != -1:
        cuts = np.all([cuts, peak == pk_0], axis = 0)
    if pk_1    != -1:
        cuts = np.all([cuts, peak == pk_1], axis = 0)
    if pk_2    != -1:
        cuts = np.all([cuts, peak == pk_2], axis = 0)
    if T_low  != -1:
        cuts = np.all([cuts, T_low < temp], axis = 0)
    if T_high  != -1:
        cuts = np.all([cuts, T_high > temp],axis = 0)
    if P_min  != -1:
        cuts = np.all([cuts, P_min < pres], axis = 0)
    if P_max  != -1:
        cuts = np.all([cuts, P_max > pres], axis = 0)
    if chi_max!= -1:
        cuts = np.all([cuts, chi2ndf < chi_max], axis = 0)
    return cuts

    # Show slow variables over time

    def all_slows_plot(incr, startdate = False, zrh = False):
        plt.rcParams['figure.figsize'] = (12.0, 16.0)    # A big plot
        f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=True, sharey=False)

        # Temperature
        x, y, mu, sigma = load_ana('temp', val_range = [1, 1e6], sigmacut = 3, incr = incr, startdate = startdate)
        do_ax(ax1, x, y, mu, sigma, 'T [$^\circ$C]', 0)

        # humidity is left out for now
        x, y, mu, sigma = load_ana('humid', val_range = [-100, 1e2], sigmacut = 5, incr = incr, startdate = startdate)
        do_ax(ax2, x, y, mu, sigma, 'H [$\%$]',0)
        #     ax.set_

        # bfield
        x, y, mu, sigma = load_ana('btot', val_range = [-1e9, 1e9], sigmacut = 1, incr = incr, startdate = startdate)
        do_ax(ax3, x, y/1e3, mu/1e3, sigma/1e3, 'B [G]',2)

        # pres
        x, y, mu, sigma = load_ana('pres', val_range = [1, 1e6], sigmacut = 10, incr = incr, startdate = startdate)
        do_ax(ax4, x, y/1.0e5, mu/1.0e5, sigma/1.0e5, 'P [bar]', 1)

        if zrh==False:
            # High voltage of channel 2
            x, y, mu, sigma = load_ana('hv2', val_range = [1, 1e6], sigmacut = 10, incr = incr, startdate = startdate)
            do_ax(ax5, x, y, mu, sigma, 'HV [V]', 2)
        y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax5.yaxis.set_major_formatter(y_formatter)

        # Fine-tune figure; make subplots close to each other and hide x ticks for
        # all but bottom plot.
        f.subplots_adjust(hspace=0.)
        f.autofmt_xdate()
        # Remove some ticks to make the plot nicer
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.xlim(x[0], x[-1])

        # Save figure
        #plt.savefig(droppath +'slows.png', dpi = 100, bbox_inches='tight') # bbox_extra_artists=(lgnd,)
        #plt.savefig(droppath +'slows.pdf', dpi = 300, bbox_inches='tight') # bbox_extra_artists=(lgnd,)
        plt.show()
        plt.rcParams['figure.figsize'] = (12.0, 10.0)    # ... and big plots

        # this is the likelihood model of the data

        def fit_A_tau_a_phi(x0, x1, y, yerr, tau12, dtau12, chan, pk, ndim = 4, T = 1, fit_a = False,
                nwalkers = 200, nsteps = 300, burnin = 50, plot_level = ['guess', 'walkers', 'corner', 'result'],
                save = False):
            '''Given a data with x_start, x_end (units years) and a rate (with error), an expected halflife tau
            (and error) optimizes Rate0, and possible modulations (with amplitude a, phase phi and period T,
            with a default period of 1 year. Specify the number of parameters with ndim (2 = exp, 3 = exp &
            phi, 4 = exp, phi & a))'''

            # Estimate the parameters
            A_est   = np.mean(y[0:10])
            A_max   = 10 * A_est
            tau_est = halflife[chan]
            phi_est = 0.0

            if fit_a: a_est = fit_a
            if ndim <= 2: a_est = 0
            if not fit_a and ndim > 2: a_est = 0.001

            est_list= [A_est, tau_est, phi_est, a_est]
            label_list = ['A', 'tau', 'phi', 'a']

            x = (x0 + x1)/2

            if ndim == 2:
                print('Rateplotter::fit_A_tau_a_phi::\tFitting for tau and A')
                a, phi= 0., 0.
            elif ndim == 3:
                assert(fit_a), 'If you want to fit for given a specify fit_a'
                print('Rateplotter::fit_A_tau_a_phi::\tFitting for tau, A and phase for a = %f'%fit_a)

            if 'guess' in plot_level:
                # Plot the dataset and the estimated model.
                xl = np.linspace(0, x[-1], 50)
                fig0 = plt.errorbar(x, y, yerr = yerr, linestyle = 'None', marker = 'o', c = 'b', label = 'Data')
                plt.plot(xl, A_est * np.power(2, -xl / tau_est) * (1 + a_est * np.sin(2 * np.pi * xl / T - phi_est)),
                    "k", lw = 3, alpha = 0.5, label = 'Estimation', linestyle = '--')
                plt.xlabel("$time$")
                plt.ylabel("$rate$")



            # Set up the sampler.
            #     pos = [result["x"] * (1 + 1e-3*np.random.randn(ndim)) for i in range(nwalkers)]
            pos = [[A_est, tau_est, a_est, phi_est][0:ndim] + 1e-1 * np.random.randn(ndim) for i in range(nwalkers)]
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = (x0, x1, y, yerr, tau12, dtau12, A_max, T), threads = 20)

                # Clear and run the production chain.
            print("Rateplotter::fit_A_tau_a_phi::\tRunning MCMC...")
            sampler.run_mcmc(pos, nsteps, rstate0 = np.random.get_state())
            print("Rateplotter::fit_A_tau_a_phi::\tDone.")

            if 'walkers' in plot_level:
                fig, axes = plt.subplots(max(ndim,2), 1, sharex=True, figsize=(8, 3 * ndim))
                for i in range(ndim):
                    axes[i].plot(sampler.chain[:, :, i].T, color="k", alpha=0.1)
                    axes[i].yaxis.set_major_locator(MaxNLocator(5))
                    axes[i].axhline(est_list[i], color="g", lw=2)
                    axes[i].set_ylabel(label_list[i])
                    axes[i].axvline(burnin, color="black", lw=2)

                axes[ndim - 1].set_xlabel("step number")
                if save: plt.savefig(droppath + 'walkers_ch%i_pk%i_%i.png'%(ch,pk,ndim), dpi = 300, bbox_inches='tight')
                plt.show()


            print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))



            # Make the triangle plot.
            samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))

            if 'corner' in plot_level:
                fig = corner.corner(samples, labels=label_list[0:ndim],
                    truths=[A_est, tau_est, phi_est, a_est * 0 ][0:ndim])
                if save: plt.savefig(droppath + 'corner_ch%i_pk%i_%i.png'%(ch,pk,ndim), dpi = 300, bbox_inches='tight')
                plt.show()

            # Compute the best estimates and erros
            mcmc_list = list(map(lambda v: [v[1], v[2]-v[1], v[1]-v[0]],zip(*np.percentile(samples, [16, 50, 84], axis=0))))

            print("Rateplotter::fit_A_tau_a_phi::\tMCMC result:")
            for i in range(ndim):
                print("\t\t\t%s\t= %.3g \t+ %.3g \t- %.3g \t(estimate = %.3g)"%(
                    label_list[i], mcmc_list[i][0], mcmc_list[i][1], mcmc_list[i][2], est_list[i]))
            for i in range(ndim, len(est_list)):
                mcmc_list.append([est_list[i],0,0])
            if 'result' in plot_level:
                # Plot the dataset and the fitted model.
                plt.errorbar(x, y, yerr = yerr, linestyle = 'None', marker = 'o', c = 'b', label = 'Data')
                xl = np.linspace(x[0], x[-1], 100)
                plt.plot(xl, mcmc_list[0][0] * np.power(2, -xl / mcmc_list[1][0]) * (1 + mcmc_list[3][0] * np.sin(2 * np.pi * xl / T - mcmc_list[2][0])),
                    c = 'r', lw = 3, label = 'MCMC result', linestyle = '--')
                plt.xlabel("$time$")
                plt.ylabel("$rate$")
                lgnd1 = plt.legend( loc = 3, borderaxespad = 0., markerscale = 1, numpoints = 1)
                if save: plt.savefig(droppath + 'res_ch%i_pk%i_%i.png'%(ch,pk,ndim), dpi = 300, bbox_extra_artists=(lgnd1,), bbox_inches='tight')
                plt.show()

                y_mcmc =  mcmc_list[0][0] * np.power(2, - x / mcmc_list[1][0])
                print('Delta rate should be %.3f'%(100 * np.mean(y/y_mcmc-1)))
                plt.errorbar(x, 100*(y/y_mcmc - 1) , yerr = 100*yerr/y_mcmc, linestyle = 'None', marker = 'o', c = 'b', label = 'sub Data')
                plt.axhline(0, color = 'red', lw = 4)
                plt.xlabel("$time$")
                plt.ylabel("$\Delta rate[\%]$")
                lgnd1 = plt.legend( loc = 3, borderaxespad = 0., markerscale = 1, numpoints = 1)
                plt.show()

            return mcmc_list


        print("Rateplotter::\tdefined fit function fit_A_tau_a_phi")

        def use_minuit(x0, x1, y, yerr, tau12, dtau12, chan, pk, ndim = 4, T = 1, fit_a = False,
                nwalkers = 200, nsteps = 300, burnin = 50, plot_level = ['guess', 'walkers', 'corner', 'result'],
                save = False):
            '''Given a data with x_start, x_end (units years) and a rate (with error), an expected halflife tau
            (and error) optimizes Rate0, and possible modulations (with amplitude a, phase phi and period T,
            with a default period of 1 year. Specify the number of parameters with ndim (2 = exp, 3 = exp &
            phi, 4 = exp, phi & a))'''

            # Estimate the parameters
            A_est   = np.mean(y[0:10])
            A_max   = 10 * A_est
            tau_est = halflife[chan]
            phi_est = 0.0

            if fit_a: a_est = fit_a
            if ndim <= 2: a_est = 0
            if not fit_a and ndim > 2: a_est = 0.001

            est_list= [A_est, tau_est, phi_est, a_est]
            label_list = ['A', 'tau', 'phi', 'a']

            x = (x0 + x1)/2
            fixed_a = fit_a

            if ndim == 1:
                def f(amp): return -lnlike([amp], x0, x1, y, yerr, tau12, dtau12, T, fixed_a)
                m = Minuit(f,
                        amp       = A_est,      limit_amp = (0, A_max),          error_amp = 1,
                        errordef  = 1.0, print_level = 0)
            if ndim == 2:
                def f(amp, t12): return -lnlike([amp,t12], x0, x1, y, yerr, tau12, dtau12, T, fixed_a)
                m = Minuit(f,
                    amp       = A_est,      limit_amp = (0, A_max),          error_amp = 1,
                    t12       = tau_est,    limit_t12 = (0, 1000),          error_t12 = dtau12,
                    errordef  = 1.0, print_level = 0)
            if ndim == 3:
                def f(amp, t12, phi): return -lnlike([amp,t12, phi], x0, x1, y, yerr, tau12, dtau12, T, fixed_a)
                m = Minuit(f,
                        amp       = A_est,      limit_amp = (0, A_max),          error_amp = 5,
                        t12       = tau_est,    limit_t12 = (5, 80),            error_t12 = dtau12,
                        phi       = phi_est,    limit_phi = (-2 * np.pi, 2 * np.pi),    error_phi = 0.1,
                        errordef  = 1.0, print_level = 0)
            if ndim == 4:
                def f(amp, t12, phi, a): return -lnlike([amp,t12, phi, a], x0, x1, y, yerr, tau12, dtau12, T, fixed_a)
                m = Minuit(f,
                        amp       = A_est,      limit_amp = (0, A_max),          error_amp = 5,
                        t12       = tau_est,    limit_t12 = (5, 80),            error_t12 = dtau12,
                        phi       = phi_est,    limit_phi = (-2 * np.pi, 2 * np.pi),    error_phi = 0.1,
                        a         = a_est,      limit_a   = (-0.1, 0.1),        error_a   = 0.0001,
                        errordef  = 1.0, print_level = 0)


            m.migrad()
            m.minos()
            m.migrad()
            res = [m.values[label] for label in ['amp','t12', 'phi', 'a'][0:ndim]]


            res_list = [[m.values[label], m.errors[label], m.errors[label]] for label in ['amp','t12', 'phi', 'a'][0:ndim]]
            for i in range(ndim + 1, 4 + 1):
                res_list.append([-1, 0, 0])

            print("\tMINUIT RESULTS")
            [print("\t\t%s \t= \t %.3g, \t+/- %.3g"%(label, m.values[label], m.errors[label])) for label in ['amp','t12', 'phi', 'a'][0:ndim]]
            if ndim != 3: return res_list
            else: return m.fval


        print("Rateplotter::\tdefined fit function use_minuit")

        def binned_rate_xy(fig, ax1, chanlist,
                pk,    # Specify source (or channel) and pk
                plot_x,
                plot_y,
                plot_dy,
                fit_y,
                fit_res,               # format: [[A,da], [tau, dtau], [phi, dphi], [a, da]]
                fitted_pars = 2,       # 2-> A & t12 3 +> phi 4+>a
                sub = False,           # Plot the residuals of the rate - exponential fit
                save = False,          # Save the figure
                savename = False,      # Specify the name of the figure
                doplot = True,         # Show the plot
                incr = 1,              # The number of fits per bin
                jenkins = False,       # Only works if sub = True, plots the jenkins hypothesis in the residual
                bgrateplot = False,    # Very specific option, used to see the fitted BG rate in this source
                exp_text = False,      # Show the results of the exponential fit
                lit_val = False,       # Show the literature value of the halflife
                legend = True,         # If true, plot the legend
                coloritem = temp,      # An ana_file propperty can be used as a colorbar
                colorbar_name = 'Temperature [$^\circ$C]',# The label that goes with the colorbar
                extra_text = False     # If some extra note needs to be placed in the figure, put it here.
                ):
            '''Plots the rate of a given source and peak number using the offered rate over time. It will return a figure of the rate of this
            source over time. '''
            plt.rcParams['figure.figsize'] = (12.0, 14.0)    # big plots
            # Check the source name, if a channel is given, interpreted the corresponding source
            chan = chanlist[0]

            if type(coloritem) == np.ndarray:
                show_col = True
                std_col_bar  = np.std(coloritem)
                mean_col_bar = np.mean(coloritem)
                col_bar      = np.clip(coloritem, mean_col_bar - 2 * std_col_bar, mean_col_bar + 2 * std_col_bar)
                numcut       = sum(np.any([coloritem < mean_col_bar - 2 * std_col_bar, coloritem > mean_col_bar + 2 * std_col_bar], axis = 0))
                print('Rateplotter::binned_rate::\tClipped %i items for the colorbar = %.0e of total'%(numcut, numcut/len(col_bar)))
            else: show_col = False

            xs, ys, dys, xlist, cdata, fitted_y = [], [], [], plot_x, [], []

            # Effectively re-bin with a number of 'incr' measured rates per bin
            for i in range(0, len(xlist), incr):
                xmax = min(len(xlist), i + incr)
                xs.append(np.average(plot_x[i:xmax]))
                ys.append(np.average(plot_y[i:xmax]))
                if type(coloritem) == np.ndarray: cdata.append(np.average(col_bar[i:xmax]))
                # now the errors
                plot_dysum = 0
                for sig in plot_dy[i:xmax]: plot_dysum += (sig ** 2)
                plot_dysum  = np.sqrt(plot_dysum) / len(plot_dy[i:xmax])
                dys.append(plot_dysum)
                fitted_y.append(np.average(fit_y[i:xmax]))

            plot_x, plot_y, plot_dy, fitted_y = np.array(xs), np.array(ys), np.array(dys), np.array(fitted_y)
            xsyears = (plot_x - toffset) * s2y


            fit_res = np.transpose(fit_res)

            amp_fit, tau_fit, phi_fit, a_fit = fit_res[0]
            damp_fit, dtau_fit, dphi_fit, da_fit = (fit_res[1] + fit_res[2]) / 2
            fit_res = np.transpose(fit_res)

            # The plotting part
            plot_x_date = [datetime.datetime.fromtimestamp(np.int(x - toffset)) for x in plot_x]

            try: ax1
            except NameError:
                fig, ax1 = plt.subplots()
                print("Binned_rate_xy::\tnoaxes to plot on now making a new figure")

            # Plot the rate subtracted from the fit and express in %'s
            if sub == True:
                ax1.set_ylabel('$\Delta$ Rate [%]')
                #         ax1.errorbar(plot_x_date, plot_y, yerr = plot_dy, linestyle = 'None', markersize = 2.5, color = 'blue')
                ax1.errorbar(plot_x_date, plot_y, yerr = plot_dy, linestyle = 'None',  label = 'Data - exponential decay', color = 'blue', capsize = 5, elinewidth = 2, zorder = 0, marker = 'None' if show_col else "o")

            # Plot the rates and fit it with
            elif sub == False:
                ax1.errorbar(plot_x_date, plot_y, yerr = plot_dy, linestyle = 'None', marker = "o", markersize = 5, label = 'Data', color = 'b',  capsize = 5, elinewidth = 2, zorder = 0)
                ax1.set_ylabel('Rate [Hz]')
            if show_col: tocol = ax1.scatter(plot_x_date, plot_y, marker = "o", s = 75 if sub else 50, c = cdata, cmap = plt.cm.rainbow, zorder = 3)
            if not bgrateplot and not sub: ax1.plot(plot_x_date, fitted_y , color = 'red', linewidth = 2.0, linestyle = 'dashed', label = 'Fitted exponential decay')



            props = dict(boxstyle = 'round', facecolor = 'white', alpha = 1)
            if jenkins and sub:
                if fitted_pars > 2 and sub: ax1.plot(plot_x_date, 100 * a_fit * np.cos(2 * np.pi * (xsyears-xsyears[0])/1 + phi_fit), label = 'Cosine for T = 1 y', color='red', linewidth = 2, linestyle = '--')
                # Now compare this to jenkins findings (an amplitude of 0.034-0.015 % and a maximum at feb 1
                # +/- 14 day <- plot_dyor not taken into account here)
                jen_xs = np.linspace(int(xsyears[0]) - 1, xsyears[-1], num = 100)
                jen_xs -= int(xsyears[0]) - 1 # Shift so xs start at 0 i.e. at 1st of january
                jen_plot_y_min = sin_function(jen_xs, 0.034, (- 0.25 + 29 / 365.25) * 2 * np.pi ) # Small sin (note last arg. is phi)
                jen_plot_y_max = sin_function(jen_xs, 0.15 , (- 0.25 + 29 / 365.25) * 2 * np.pi ) # Large sin
                jen_xs += int(xsyears[0]) - 1 # Shift back
                # Cut anything that is calculated for a date smaller that the first (real) date
                jen_plot_y_min = jen_plot_y_min[jen_xs > xsyears[0]]
                jen_plot_y_max = jen_plot_y_max[jen_xs > xsyears[0]]
                jen_xs     = jen_xs[jen_xs > xsyears[0]]
                jen_exp_fit = np.power(2, -((jen_xs - jen_xs[0]) / tau_fit)) * amp_fit
                jen_xs_yr = jen_xs


                # Convert to nice dates
                jen_xs = [datetime.datetime.fromtimestamp(np.int(x * y2s)) for x in jen_xs]
                if not sub: jen_plot_y_min, jen_plot_y_max = jen_exp_fit * (1 + 0.01 * jen_plot_y_min), jen_exp_fit * (1 + 0.01 * jen_plot_y_max )
                if fitted_pars > 2 and not sub:
                    ax1.plot(plot_x_date,
                             np.power(2, -((xsyears - xsyears[0]) / tau_fit)) * amp_fit *
                             (1 + a_fit * np.cos(2 * np.pi * xsyears/1 + phi_fit)), label = 'MCMC sinus for T = 1 y', color='red', linewidth = 2, linestyle = '--')
                # And plot
                ax1.plot(jen_xs, jen_plot_y_min, label = "O'Keefe et al. 2013, amplitude \n$0.034 - 0.15\%$ max. at 29 jan", color = 'g', linestyle = '-', linewidth = 2)
                ax1.plot(jen_xs, jen_plot_y_max, color = 'g', linestyle = '-', linewidth = 2)
                ax1.fill_between(jen_xs, jen_plot_y_min, jen_plot_y_max, color= 'g', alpha = 0.2) # nice region
                if sub: ax1.plot(jen_xs, jen_plot_y_min * 0, c = 'r', label = 'No modulation', linewidth = 2)

                # Add a little textbox about the fitresults $\propto 1 + a \cdot\cos(2\pi t/T-\phi)$
                if fitted_pars > 2:
                    textstr = ("fitted model:\n$a = %.2G \pm %.2G$"%(100 * a_fit, 100 * da_fit) + "%"+"\n$\phi= %.2G \pm %.2G$\n$T=1$ yr" %(phi_fit, dphi_fit))
                    ax1.text(1.05 +.25*show_col, 1-0.5*legend, textstr, color = 'r', transform = ax1.transAxes, fontsize = 24, verticalalignment='top', bbox=props)

                if legend: lgnd1 = ax1.legend(bbox_to_anchor = (1.05 +.25*show_col, 1), loc = 2, borderaxespad = 0., markerscale = 1, numpoints = 1)
            elif legend:   lgnd1 = ax1.legend(bbox_to_anchor = (1.05 +.25*show_col, 1), loc = 2, borderaxespad = 0., markerscale = 1, numpoints = 1)

            if show_col:
                cb1 = fig.colorbar(tocol, ax=[ax1] , orientation = 'vertical')#, fig.add_axes([0.85, 0.15, 0.05, 0.7])

                cb1.set_label(colorbar_name)

            plt.xticks(rotation=45)
            if np.max(abs(plot_y) + abs(plot_dy)) < 0.15: ax1.set_ylim(-0.15, 0.15)

            textstr = 'Source: %s at %.1f keV \n$t_{1/2}=%.5g\pm %.2g^{stat} $ year' % (sourceLatex[chan], sourceE[chan][pk], tau_fit, dtau_fit)
            if lit_val: textstr += '\nLiterature: $%.5g\pm%.2g$ year'%(halflife[chan], dhalflife[chan])
            if exp_text and chan >= -2: ax1.text(1.05 +.25*show_col, 0.4 if sub else 1 -0.5*legend, textstr,  transform = ax1.transAxes, fontsize = 24, verticalalignment='top', bbox=props)
            if extra_text: ax1.text(0.05, 0.95, extra_text, transform = ax1.transAxes, fontsize = 24, verticalalignment='top', bbox=props)
            if save: fig.savefig(droppath + savename, dpi = 300, bbox_extra_artists=(lgnd1 if legend else "",), bbox_inches='tight')

            if doplot: plt.show()
