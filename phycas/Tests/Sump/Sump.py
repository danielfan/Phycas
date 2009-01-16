from phycas import *

if __name__ == '__main__':
    sump.out.log            = 'logfile.txt'
    sump.out.log.mode       = REPLACE
    sump.file               = 'test.p'
    #sump.file               = 'k2.n10k.p'
    sump.burnin             = 1
    sump()
    
