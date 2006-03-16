29-Nov-2005
-----------

Timings and parameter estimates for the nyldna4.nex dataset:

          PAUP* 4b10      MrBayes 3.1   MCMCSimple 0.1  MCMCSimple 0.2
        +--------------+---------------+---------------+--------------+
lnL     | -7124.42582  |   -7107.889   |   -7050.815   |  -7128.651   |
        +--------------+---------------+---------------+--------------+
TL      |   0.22693    |      0.218    |      0.225    |     0.229    |
        +--------------+---------------+---------------+--------------+
kappa   |   1.86054    |      1.854    |      1.844    |     1.876    |
        +--------------+---------------+---------------+--------------+
A       |   0.284433   |      0.285    |      0.281    |     0.284    |
        +--------------+---------------+---------------+--------------+
C       |   0.177087   |      0.176    |      0.175    |     0.179    |
        +--------------+---------------+---------------+--------------+
G       |   0.207635   |      0.209    |      0.206    |     0.207    |
        +--------------+---------------+---------------+--------------+
T       |   0.330845   |      0.33     |      0.327    |     0.33     |
        +--------------+---------------+---------------+--------------+
hyper   |     ---      |      ---      |      0.183    |     0.23     |
        +--------------+---------------+---------------+--------------+
evals   |     ---      |   1,000,001   |    742,018    |   786,597    |
        +--------------+---------------+---------------+--------------+
secs    |     ---      |     51.377    |     54.375    |    32.815    |
        +--------------+---------------+---------------+--------------+
ev/sec  |     ---      |     19464     |     13646     |    23971     |
        +--------------+---------------+---------------+--------------+

The difference between MCMCSimple versions 0.1 and 0.2 are that the 
ParamManager and the MCMC param classes (KappaParam, EdgeLenParam,
HyperPriorParam and BaseFreqParam) were moved from Python to C++. The
python versions of these were slowing things down because the slice
samplers used their __call__ function many times. Now these calls are
to operator() and do not cross the Python/C++ boundary, resulting in
a speedup of approximately 75%).

Version 0.2 was committed to CVS about 2pm 29-Nov-2005 (neither version
was tagged, however, because a lot of cleaning up still needs to be done).

