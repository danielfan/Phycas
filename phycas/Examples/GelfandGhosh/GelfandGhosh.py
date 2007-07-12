# This example program simulates data from a four-taxon tree using the JC model,
# drawing edge lengths in the model tree from an exponential distribution with mean
# 0.02. These data are then analyzed by a model having a different edge length
# prior mean and the Gelfand-Ghosh measure is computed so that the model used 
# for the analysis can be compared to other models.
#
# This is a more advanced example that involves more than simply changing
# default values and saying go. If this is the first example you have looked at,
# you might want to consider looking at Simplest or Paradox first.

from phycas import *
from math import exp

# Here are some numbers we will use later on in this script. They are defined
# here to make it easy to find later should we wish to run the script again
# with different values.
rnseed = 7654321
num_sites = 2000
ncycles = 20000
true_edgelen_mean = 0.1
analysis_prior_means = [0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0]

#############################################################
#################### Simulate a data set ####################
#############################################################
def simulate(file_name):
    # Create a Phycas object named simulator. Note that simulator is local to this
    # function and will be destroyed before the function returns. This is ok because
    # by that point the simulated data will be saved in a file.
    simulator = Phycas()

    # Create a pseudorandom number generator
    rng = ProbDist.Lot()
    rng.setSeed(rnseed)

    # Set up the substitution model that will be used for both simulation and the MCMC analysis
    jcmodel = Likelihood.JCModel() # use the Jukes-Cantor (1969) model
    jcmodel.setNGammaRates(1)      # specifying one rate means *no* rate heterogeneity

    # Create a tree
    true_tree = Phylogeny.Tree()
    true_tree.buildFromString('(1,2,(3,4))')

    # Create a true edge length probability distribution
    # If setLot were not called, true_edgelen_dist would use its own random number generator
    # and the edge lengths generated would be different each time this script is run
    true_edgelen_dist = ProbDist.ExponentialDist(1.0/true_edgelen_mean)
    true_edgelen_dist.setLot(rng)

    # Now use true_edgelen_dist to assign random edge lengths to the tree
    tm = Phylogeny.TreeManip(true_tree)
    tm.setRandomEdgeLengths(true_edgelen_dist)

    # Create a likelihood object to orchestrate both simulations and likelihood calculations
    simulator.likelihood = Likelihood.TreeLikelihood(jcmodel)

    # Prepare the tree for simulation (i.e. equip nodes with transition matrices)
    simulator.likelihood.prepareForSimulation(true_tree)

    # Create a SimData object to hold the simulated data
    sim_data = Likelihood.SimData()

    # Simulate num_sites of data and store in sim_data
    # Use the function simulateFirst (rather than just simulate) in order
    # to force calculation of transition probabilities
    simulator.likelihood.simulateFirst(sim_data, true_tree, rng, num_sites)

    # Define the names of the taxa to use when the simulated data set is saved to a file
    taxon_labels = ['A', 'B', 'C', 'D']

    # Save simulated data to a NEXUS file using taxon_labels, datatype=dna and
    # using the symbols a, c, g, and t to represent the nucleotides
    sim_data.saveToNexusFile(file_name, taxon_labels, 'dna', ('a','c','g','t'))

    # Add a paup block to the end of the simulated data file to make it easy to
    # verify that the simulation worked
    sim_file = file(file_name, 'a') # the 'a' means open for appending
    sim_file.write('\nbegin trees;\n')
    sim_file.write('  utree truth = %s;\n' % (true_tree.makeNewick()))
    sim_file.write('end;\n\n')
    sim_file.write('begin paup;\n')
    sim_file.write('  set criterion=likelihood;\n')
    sim_file.write('  lset nst=1 basefreq=equal rates=equal;\n')
    sim_file.write('  describetrees 1 / plot=phylogram brlens;\n')
    sim_file.write('end;\n')

####################################################################
#################### Analyze the simulated data ####################
####################################################################
def analyze(prior_mean):
    # Create a Phycas object named analyzer. Note that analyzer is local to this
    # function and will be destroyed before the analyze function returns.
    analyzer = Phycas()

    # Set various options before starting the MCMC run. To see all available options, you
    # can go to the source: take a look at the beginning part (the __init__ function)
    # of the Phycas.py file
    #   phycas/Examples/GelfandGhosh/GelfandGhosh.py  <-- you are here
    #   phycas/Phycas/Phycas.py                       <-- here is Phycas.py

    # Use the same random number seed in the analysis that we used for the simulation
    analyzer.random_seed = rnseed

    # Specify the model to be used for the analysis
    analyzer.default_model = 'jc' # use the Jukes-Cantor (1969) model
    analyzer.num_rates = 1        # specifying 1 rate means no rate heterogeneity

    # The JC model has only edge length parameters, so the only priors we need to worry
    # about are the edge length priors. The master_edgelen_dist specifies the prior used
    # for all edge lengths. Tell Phycas to not use a hyperprior for edge lengths, otherwise
    # it will ignore your master_edgelen_dist specification.
    analyzer.using_hyperprior = False 
    analyzer.master_edgelen_dist = ProbDist.ExponentialDist(1.0/prior_mean)

    # Don't allow polytomies
    analyzer.allow_polytomies = False

    # Specify the data file
    analyzer.data_file_name = 'simulated.nex'

    # Start with a random tree
    analyzer.starting_tree_source = 'random'

    # Tell analyzer that we want to run the MCMC analysis for 20000 cycles.
    analyzer.ncycles = ncycles
    analyzer.sample_every = 10          # save tree and parameters every 10 cycles
    analyzer.report_every = ncycles/20  # report progress only 20 times during the run

    # Tell analyzer to calculate the Gelfand-Ghosh measure as the MCMC analysis runs
    # gg_kvect is a list containing all the values of k for which you wish to calculate the
    # Gelfand-Ghosh measure. The value of k is the relative importance given to goodness-of-fit
    # as opposed to predictive variance. Here we weight them equally by supplying only the
    # value 1.0 for k.
    analyzer.gg_do = True       # if False, no posterior predictive simulations would be performed 
    analyzer.gg_kvect = [1.0]   # could be, e.g. [0.5, 1.0, 2.0] if you wanted to try three different k values
    analyzer.gg_outfile = None  # do not need to save a file of Gelfand-Ghosh results

    # Finally, call the mcmc method, which prepares and then runs the MCMC sampler
    analyzer.mcmc()

    # This function returns a tuple (a type of list, in this case consisting of two values)
    # The first element of the tuple is Pm, and the other value is Gm. The sum of these
    # two values is GGCm, the overall Gelfand-Ghosh value, but this need not be saved.
    # The [0] after gg_Gm is necessary because Gm is different for each value of k.
    # We specified only one value for k (i.e. 1.0), so only one value of Gm has been
    # computed and gg_Gm is thus a list containing only one element; nevertheless,
    # because gg_Gm is a list, we must ask for the first (zeroth) element.
    return (analyzer.gg_Pm, analyzer.gg_Gm[0])

if __name__ == '__main__':
    # Call the simulate function (defined above) to simulate a data set and save it
    # in the file named 'simulated.nex'
    simulate('simulated.nex')

    # Open a file that will save the Gelfand-Ghosh values
    outf = file('output.txt', 'w')  # the 'w' means "open the file for writing"
    outf.write('prior mean\tPm\tGm\tGGCm\n')    # '\t' means tab, '\n' means newline

    # For each value of prior mean, perform an MCMC analysis
    for m in analysis_prior_means:
        Pm, Gm = analyze(m)
        GGCm = Pm + Gm
        print '**'
        print '** Analysis with prior mean',m,'yielded Pm =',Pm,', Gm =',Gm
        print '**'

        # Write out a line comprising four values separated by tabs
        # Each %f represents one floating point value from the 4-valued tuple that follows
        # Each '\t' is a tab character, and the final '\n' starts a new line in the file
        outf.write('%f\t%f\t%f\t%f\n' % (m, Pm, Gm, GGCm))

    # Close the file and we're done!
    outf.close()
