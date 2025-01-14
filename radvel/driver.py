"""
Driver functions for the radvel pipeline.\
These functions are meant to be used only with\
the `cli.py` command line interface.
"""
from __future__ import print_function

import radvel
from radvel.likelihood import GPLikelihood
from radvel.plot import orbit_plots, mcmc_plots
from radvel.mcmc import statevars

import os
import sys
import copy
import collections
from collections import OrderedDict

import pandas as pd
import numpy as np
from numpy import inf
from astropy import constants as c

if sys.version_info[0] < 3:
    import ConfigParser as configparser
else:
    import configparser
    from configparser import NoSectionError


def plots(args):
    """
    Generate plots

    Args:
        args (ArgumentParser): command line arguments
    """

    config_file = args.setupfn
    conf_base = os.path.basename(config_file).split('.')[0]
    statfile = os.path.join(
        args.outputdir, "{}_radvel.stat".format(conf_base)
    )

    status = load_status(statfile)

    P, post = radvel.utils.initialize_posterior(config_file,
                                                decorr=args.decorr)

    assert status.getboolean('fit', 'run'), \
        "Must perform max-liklihood fit before plotting"
    post = radvel.posterior.load(status.get('fit', 'postfile'))

    for ptype in args.type:
        print("Creating {} plot for {}".format(ptype, conf_base))

        if ptype == 'rv':
            args.plotkw['uparams'] = post.uparams
            args.plotkw['status'] = status
            if 'saveplot' not in args.plotkw:
                saveto = os.path.join(
                    args.outputdir, conf_base+'_rv_multipanel.pdf'
                )
            else:
                saveto = args.plotkw['saveplot']
                args.plotkw.pop('saveplot')
            P, _ = radvel.utils.initialize_posterior(config_file)
            if hasattr(P, 'bjd0'):
                args.plotkw['epoch'] = P.bjd0

            if args.gp or isinstance(post.likelihood, GPLikelihood):
                GPPlot = orbit_plots.GPMultipanelPlot(
                    post, saveplot=saveto, **args.plotkw
                )
                GPPlot.plot_multipanel()
            else:
                RVPlot = orbit_plots.MultipanelPlot(
                    post, saveplot=saveto, **args.plotkw
                )
                RVPlot.plot_multipanel()

                # check to make sure that Posterior is not GP, print warning if it is
                if isinstance(post.likelihood, radvel.likelihood.CompositeLikelihood):
                    like_list = post.likelihood.like_list
                else:
                    like_list = [post.likelihood]
                for like in like_list:
                    if isinstance(like, radvel.likelihood.GPLikelihood):
                        print("WARNING: GP Likelihood(s) detected. \
You may want to use the '--gp' flag when making these plots.")
                        break

        if ptype == 'corner' or ptype == 'auto' or ptype == 'trend':
            assert status.getboolean('mcmc', 'run'), \
                "Must run MCMC before making corner, auto, or trend plots"

            chains = pd.read_csv(status.get('mcmc', 'chainfile'))
            autocorr = pd.read_csv(status.get('mcmc', 'autocorrfile'))

        if ptype == 'auto':
            saveto = os.path.join(args.outputdir, conf_base+'_auto.pdf')
            Auto = mcmc_plots.AutoPlot(autocorr, saveplot=saveto)
            Auto.plot()

        if ptype == 'corner':
            saveto = os.path.join(args.outputdir, conf_base+'_corner.pdf')
            Corner = mcmc_plots.CornerPlot(post, chains, saveplot=saveto)
            Corner.plot()

        if ptype == 'trend':
            nwalkers = status.getint('mcmc', 'nwalkers')
            nensembles = status.getint('mcmc', 'nensembles')

            saveto = os.path.join(args.outputdir, conf_base+'_trends.pdf')
            Trend = mcmc_plots.TrendPlot(post, chains, nwalkers, nensembles, saveto)
            Trend.plot()

        if ptype == 'derived':
            assert status.has_section('derive'), \
                "Must run `radvel derive` before plotting derived parameters"

            P, _ = radvel.utils.initialize_posterior(config_file)
            chains = pd.read_csv(status.get('derive', 'chainfile'))
            saveto = os.path.join(
                args.outputdir, conf_base+'_corner_derived_pars.pdf'
            )

            Derived = mcmc_plots.DerivedPlot(chains, P, saveplot=saveto)
            Derived.plot()

        savestate = {'{}_plot'.format(ptype): os.path.relpath(saveto)}
        save_status(statfile, 'plot', savestate)


def fit(args):
    """Perform maximum a posteriori fit

    Args:
        args (ArgumentParser): command line arguments
    """

    config_file = args.setupfn
    conf_base = os.path.basename(config_file).split('.')[0]
    print("Performing maximum a posteriori fitting for {}".format(conf_base))

    P, post = radvel.utils.initialize_posterior(config_file, decorr=args.decorr)

    post = radvel.fitting.maxlike_fitting(post, verbose=True)

    postfile = os.path.join(args.outputdir,
                            '{}_post_obj.pkl'.format(conf_base))
    post.writeto(postfile)

    savestate = {'run': True,
                 'postfile': os.path.relpath(postfile)}
    save_status(os.path.join(args.outputdir,
                '{}_radvel.stat'.format(conf_base)),
                'fit', savestate)


def mcmc(args):
    """Perform MCMC error analysis

    Args:
        args (ArgumentParser): command line arguments
    """

    config_file = args.setupfn
    conf_base = os.path.basename(config_file).split('.')[0]
    statfile = os.path.join(args.outputdir,
                            "{}_radvel.stat".format(conf_base))

    if args.save or args.proceed:
        backend_loc = os.path.join(args.outputdir, conf_base+'_rawchain.h5')
    else:
        backend_loc = None

    status = load_status(statfile)
    P, post = radvel.utils.initialize_posterior(config_file,
                                                decorr=args.decorr)

    if status.getboolean('fit', 'run'):
        print("Loading starting positions from previous MAP fit")

        post = radvel.posterior.load(status.get('fit', 'postfile'))

    msg1 = (
            "Running MCMC for {}, N_walkers = {}, N_steps = {}, N_ensembles = {}, Min Auto Factor = {}, "
            ).format(conf_base, args.nwalkers, args.nsteps, args.ensembles, args.minAfactor)

    msg2 = (
            "Max Auto Relative-Change = {}, Max G-R = {}, Min Tz = {} ..."
            ).format(args.maxArchange, args.maxGR, args.minTz)

    print(msg1 + '\n' + msg2)

    chains = radvel.mcmc(post, nwalkers=args.nwalkers, nrun=args.nsteps, ensembles=args.ensembles,
                         minAfactor=args.minAfactor, maxArchange=args.maxArchange, burnAfactor=args.burnAfactor,
                         burnGR=args.burnGR, maxGR=args.maxGR, minTz=args.minTz, minsteps=args.minsteps,
                         minpercent=args.minpercent, thin=args.thin, serial=args.serial, save=args.save,
                         savename=backend_loc, proceed=args.proceed, proceedname=backend_loc, headless=args.headless)

    mintz = statevars.mintz
    maxgr = statevars.maxgr
    minafactor = statevars.minafactor
    maxarchange = statevars.maxarchange

    # Convert chains into synth basis
    synthchains = chains.copy()
    for par in post.params.keys():
        if not post.vector.vector[post.vector.indices[par]][1]:
            synthchains[par] = post.vector.vector[post.vector.indices[par]][0]

    synthchains = post.params.basis.to_synth(synthchains)
    synth_quantile = synthchains.quantile([0.159, 0.5, 0.841])

    # Get quantiles and update posterior object to median
    # values returned by MCMC chains
    post_summary = chains.quantile([0.159, 0.5, 0.841])

    for k in chains.columns:
        if k in post.params.keys():
            post.vector.vector[post.vector.indices[k]][0] = post_summary[k][0.5]

    post.vector.vector_to_dict()

    print("Performing post-MCMC maximum likelihood fit...")
    post = radvel.fitting.maxlike_fitting(post, verbose=False)

    final_logprob = post.logprob()
    final_residuals = post.likelihood.residuals().std()
    final_chisq = np.sum(post.likelihood.residuals()**2 / (post.likelihood.errorbars()**2))
    deg_of_freedom = len(post.likelihood.y) - len(post.likelihood.get_vary_params())
    final_chisq_reduced = final_chisq / deg_of_freedom
    post.vector.vector_to_dict()
    synthparams = post.params.basis.to_synth(post.params)

    print("Calculating uncertainties...")
    post.uparams = {}
    post.medparams = {}
    post.maxparams = {}
    for par in synthparams.keys():
        maxlike = synthparams[par].value
        med = synth_quantile[par][0.5]
        high = synth_quantile[par][0.841] - med
        low = med - synth_quantile[par][0.159]
        err = np.mean([high, low])
        if maxlike == -np.inf and med == -np.inf and np.isnan(low) and np.isnan(high):
            err = 0.0
        else:
            err = radvel.utils.round_sig(err)
        if err > 0.0:
            med, err, errhigh = radvel.utils.sigfig(med, err)
            maxlike, err, errhigh = radvel.utils.sigfig(maxlike, err)
        post.uparams[par] = err
        post.medparams[par] = med
        post.maxparams[par] = maxlike

    print("Final loglikelihood = %f" % final_logprob)
    print("Final RMS = %f" % final_residuals)
    print("Final reduced chi-square = {}".format(final_chisq_reduced))
    print("Best-fit parameters:")
    print(post)

    print("Saving output files...")
    saveto = os.path.join(args.outputdir, conf_base+'_post_summary.csv')
    post_summary.to_csv(saveto, sep=',')

    postfile = os.path.join(args.outputdir,
                            '{}_post_obj.pkl'.format(conf_base))
    post.writeto(postfile)

    csvfn = os.path.join(args.outputdir, conf_base+'_chains.csv.bz2')
    chains.to_csv(csvfn, compression='bz2')

    auto = pd.DataFrame()
    auto['autosamples'] = statevars.autosamples
    auto['automin'] = statevars.automin
    auto['automean'] = statevars.automean
    auto['automax'] = statevars.automax
    auto['factor'] = statevars.factor
    autocorr = os.path.join(args.outputdir, conf_base+'_autocorr.csv')
    auto.to_csv(autocorr, sep=',')

    savestate = {'run': True,
                 'postfile': os.path.relpath(postfile),
                 'chainfile': os.path.relpath(csvfn),
                 'autocorrfile': os.path.relpath(autocorr),
                 'summaryfile': os.path.relpath(saveto),
                 'nwalkers': statevars.nwalkers,
                 'nensembles': args.ensembles,
                 'maxsteps': args.nsteps*statevars.nwalkers*args.ensembles,
                 'nsteps': statevars.ncomplete,
                 'nburn': statevars.nburn,
                 'minafactor': minafactor,
                 'maxarchange': maxarchange,
                 'minTz': mintz,
                 'maxGR': maxgr}
    save_status(statfile, 'mcmc', savestate)

    statevars.reset()


def ic_compare(args):
    """Compare different models and comparative statistics including
          AIC and BIC statistics.

    Args:
        args (ArgumentParser): command line arguments
    """

    config_file = args.setupfn
    conf_base = os.path.basename(config_file).split('.')[0]
    statfile = os.path.join(args.outputdir,
                            "{}_radvel.stat".format(conf_base))

    status = load_status(statfile)

    P, post = radvel.utils.initialize_posterior(config_file,
                                                decorr=args.decorr)

    assert status.getboolean('fit', 'run'), \
        "Must perform max-liklihood fit before running Information Criteria comparisons"
    post = radvel.posterior.load(status.get('fit', 'postfile'))

    choices = ['nplanets', 'e', 'trend', 'jit', 'gp']
    statsdictlist = []
    paramlist = []
    compareparams = args.type

    ipost = copy.deepcopy(post)

    if args.simple:
        statsdictlist += radvel.fitting.model_comp(ipost, params=[], verbose=args.verbose)
    else:
        if hasattr(args, 'fixjitter') and args.fixjitter:
            for param in ipost.params:
                if len(param) >= 3 and param[0:3] == 'jit':
                    ipost.params[param].vary = False

        for compareparam in compareparams:
            assert compareparam in choices, \
                "Valid parameter choices for 'ic -t' are combinations of: "\
                + " ".join(choices)
            paramlist.append(compareparam)
            if hasattr(args, 'mixed') and not args.mixed:
                statsdictlist += radvel.fitting.model_comp(ipost, params=[compareparam], verbose=args.verbose)
        if hasattr(args, 'mixed') and not args.mixed:
            new_statsdictlist = []
            for dicti in statsdictlist:
                anymatch = False
                for seendict in new_statsdictlist:
                    if collections.Counter(dicti['Free Params'][0]) == \
                            collections.Counter(seendict['Free Params'][0]):
                        anymatch = True
                        continue
                if not anymatch:
                    new_statsdictlist.append(dicti)
            statsdictlist = new_statsdictlist

        if not hasattr(args, 'mixed') or (hasattr(args, 'mixed') and args.mixed):
            statsdictlist += radvel.fitting.model_comp(ipost, params=paramlist, verbose=args.verbose)

    savestate = {'ic': statsdictlist}
    save_status(statfile, 'ic_compare', savestate)


def tables(args):
    """Generate TeX code for tables in summary report

    Args:
        args (ArgumentParser): command line arguments
    """

    config_file = args.setupfn
    conf_base = os.path.basename(config_file).split('.')[0]
    statfile = os.path.join(args.outputdir,
                            "{}_radvel.stat".format(conf_base))
    status = load_status(statfile)

    assert status.getboolean('mcmc', 'run'), \
        "Must run MCMC before making tables"

    P, post = radvel.utils.initialize_posterior(config_file)
    post = radvel.posterior.load(status.get('fit', 'postfile'))
    chains = pd.read_csv(status.get('mcmc', 'chainfile'))
    minafactor = status.get('mcmc', 'minafactor')
    maxarchange = status.get('mcmc', 'maxarchange')
    maxgr = status.get('mcmc', 'maxgr')
    mintz = status.get('mcmc', 'mintz')
    if 'derive' in status.sections() and status.getboolean('derive', 'run'):
        dchains = pd.read_csv(status.get('derive', 'chainfile'))
        chains = chains.join(dchains, rsuffix='_derived')
        derived = True
    else:
        derived = False
    report = radvel.report.RadvelReport(P, post, chains, minafactor, maxarchange, maxgr, mintz, derived=derived)
    tabletex = radvel.report.TexTable(report)
    attrdict = {'priors': 'tab_prior_summary', 'rv': 'tab_rv',
                'params': 'tab_params', 'derived': 'tab_derived',
                'crit': 'tab_crit'}
    for tabtype in args.type:
        print("Generating LaTeX code for {} table".format(tabtype))

        if tabtype == 'ic_compare':
            assert status.has_option('ic_compare', 'ic'), \
                "Must run Information Criteria comparison before making comparison tables"

            compstats = eval(status.get('ic_compare', 'ic'))
            report = radvel.report.RadvelReport(
                P, post, chains, minafactor, maxarchange, maxgr, mintz, compstats=compstats
            )
            tabletex = radvel.report.TexTable(report)
            tex = tabletex.tab_comparison()
        elif tabtype == 'rv':
            tex = getattr(tabletex, attrdict[tabtype])(name_in_title=args.name_in_title, max_lines=None)
        elif tabtype == 'crit':
            tex = getattr(tabletex, attrdict[tabtype])(name_in_title=args.name_in_title)
        else:
            if tabtype == 'derived':
                assert status.has_option('derive', 'run'), \
                    "Must run `radvel derive` before making derived parameter table"
            assert tabtype in attrdict, 'Invalid Table Type %s ' % tabtype
            tex = getattr(tabletex, attrdict[tabtype])(name_in_title=args.name_in_title)

        saveto = os.path.join(
            args.outputdir, '{}_{}.tex'.format(conf_base, tabtype)
        )
        with open(saveto, 'w+') as f:
            f.write(tex)

        savestate = {'{}_tex'.format(tabtype): os.path.relpath(saveto)}
        save_status(statfile, 'table', savestate)


def derive(args):
    """Derive physical parameters from posterior samples

    Args:
        args (ArgumentParser): command line arguments
    """

    config_file = args.setupfn
    conf_base = os.path.basename(config_file).split('.')[0]
    statfile = os.path.join(args.outputdir,
                            "{}_radvel.stat".format(conf_base))
    status = load_status(statfile)

    msg = "Multiplying mcmc chains by physical parameters for {}".format(
        conf_base
    )
    print(msg)

    assert status.getboolean('mcmc', 'run'), "Must run MCMC before making tables"

    P, post = radvel.utils.initialize_posterior(config_file)
    post = radvel.posterior.load(status.get('fit', 'postfile'))
    chains = pd.read_csv(status.get('mcmc', 'chainfile'))

    try:
        mstar = np.random.normal(
            loc=P.stellar['mstar'], scale=P.stellar['mstar_err'],
            size=len(chains)
            )
    except AttributeError:
        print("Unable to calculate derived parameters, stellar parameters not defined the config file.")
        return

    if (mstar <= 0.0).any():
        num_nan = np.sum(mstar <= 0.0)
        nan_perc = float(num_nan) / len(chains)
        mstar[mstar <= 0] = np.abs(mstar[mstar <= 0])
        print("WARNING: {} ({:.2f} %) of Msini samples are NaN. The stellar mass posterior may contain negative \
values. Interpret posterior with caution.".format(num_nan, nan_perc))

    # Convert chains into synth basis
    synthchains = chains.copy()
    for par in post.params.keys():
        if not post.params[par].vary:
            synthchains[par] = post.params[par].value

    synthchains = post.params.basis.to_synth(synthchains)

    savestate = {'run': True}
    outcols = []
    for i in np.arange(1, P.nplanets + 1, 1):
        # Grab parameters from the chain
        def _has_col(key):
            cols = list(synthchains.columns)
            return cols.count('{}{}'.format(key, i)) == 1

        def _get_param(key):
            if _has_col(key):
                return synthchains['{}{}'.format(key, i)]
            else:
                return P.params['{}{}'.format(key, i)].value

        def _set_param(key, value):
            chains['{}{}'.format(key, i)] = value

        def _get_colname(key):
            return '{}{}'.format(key, i)

        per = _get_param('per')
        k = _get_param('k')
        e = _get_param('e')

        mpsini = radvel.utils.Msini(k, per, mstar, e, Msini_units='earth')
        _set_param('mpsini', mpsini)
        outcols.append(_get_colname('mpsini'))

        mtotal = mstar + (mpsini * c.M_earth.value) / c.M_sun.value      # get total star plus planet mass
        a = radvel.utils.semi_major_axis(per, mtotal)               # changed from mstar to mtotal
        
        _set_param('a', a)
        outcols.append(_get_colname('a'))

        musini = (mpsini * c.M_earth) / (mstar * c.M_sun)
        _set_param('musini', musini)
        outcols.append(_get_colname('musini'))

        try:
            rp = np.random.normal(
                loc=P.planet['rp{}'.format(i)],
                scale=P.planet['rp_err{}'.format(i)],
                size=len(chains)
            )

            _set_param('rp', rp)
            _set_param('rhop', radvel.utils.density(mpsini, rp))

            outcols.append(_get_colname('rhop'))
        except (AttributeError, KeyError):
            pass

    # Get quantiles and update posterior object to median
    # values returned by MCMC chains
    quantiles = chains.quantile([0.159, 0.5, 0.841])
    csvfn = os.path.join(args.outputdir, conf_base+'_derived_quantiles.csv')
    quantiles.to_csv(csvfn, columns=outcols)

    # saved derived paramters to posterior file
    postfile = os.path.join(args.outputdir,
                            '{}_post_obj.pkl'.format(conf_base))
    post.derived = quantiles[outcols]
    post.writeto(postfile)
    savestate['quantfile'] = os.path.relpath(csvfn)

    csvfn = os.path.join(args.outputdir, conf_base+'_derived.csv.bz2')
    chains.to_csv(csvfn, columns=outcols, compression='bz2')
    savestate['chainfile'] = os.path.relpath(csvfn)

    print("Derived parameters:", outcols)

    save_status(statfile, 'derive', savestate)


def report(args):
    """Generate summary report

    Args:
        args (ArgumentParser): command line arguments
    """

    config_file = args.setupfn
    conf_base = os.path.basename(config_file).split('.')[0]
    print("Assembling report for {}".format(conf_base))

    statfile = os.path.join(args.outputdir,
                            "{}_radvel.stat".format(conf_base))

    status = load_status(statfile)

    P, post = radvel.utils.initialize_posterior(config_file)
    post = radvel.posterior.load(status.get('fit', 'postfile'))
    chains = pd.read_csv(status.get('mcmc', 'chainfile'))
    minafactor = status.get('mcmc', 'minafactor')
    maxarchange = status.get('mcmc', 'maxarchange')
    maxgr = status.get('mcmc', 'maxgr')
    mintz = status.get('mcmc', 'mintz')
    if 'derive' in status.sections() and status.getboolean('derive', 'run'):
        dchains = pd.read_csv(status.get('derive', 'chainfile'))
        chains = chains.join(dchains, rsuffix='_derived')
        derived = True
    else:
        derived = False
    try:
        compstats = eval(status.get('ic_compare', 'ic'))
    except (KeyError, NoSectionError):
        print("WARNING: Could not find {} model comparison \
in {}.\nPlease make sure that you have run `radvel ic -t {}` (or, e.g., `radvel \
ic -t nplanets e trend jit gp`)\
\nif you would like to include the model comparison table in the \
report.".format(args.comptype,
                statfile,
                args.comptype))
        compstats = None

    report = radvel.report.RadvelReport(P, post, chains, minafactor, maxarchange, maxgr, mintz, compstats=compstats,
                                        derived=derived)
    report.runname = conf_base

    report_depfiles = []
    for ptype, pfile in status.items('plot'):
        report_depfiles.append(pfile)

    with radvel.utils.working_directory(args.outputdir):
        rfile = os.path.join(conf_base+"_results.pdf")
        report_depfiles = [os.path.basename(p) for p in report_depfiles]
        report.compile(
            rfile, depfiles=report_depfiles, latex_compiler=args.latex_compiler
        )


def save_status(statfile, section, statevars):
    """Save pipeline status

    Args:
        statfile (string): name of output file
        section (string): name of section to write
        statevars (dict): dictionary of all options to populate
           the specified section
    """

    config = configparser.RawConfigParser()

    if os.path.isfile(statfile):
        config.read(statfile)

    if not config.has_section(section):
        config.add_section(section)

    for key, val in statevars.items():
        config.set(section, key, val)

    with open(statfile, 'w') as f:
        config.write(f)


def load_status(statfile):
    """Load pipeline status

    Args:
        statfile (string): name of configparser file

    Returns:
        configparser.RawConfigParser
    """

    config = configparser.RawConfigParser()
    gl = config.read(statfile)

    return config
