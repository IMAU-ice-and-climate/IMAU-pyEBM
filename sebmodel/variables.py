from dataclasses import dataclass
import numpy as np
from numba import int32, float64    # import the types
import numpy as np

from numba.experimental import jitclass

# from globals import *

spec = [
    ('dzdeep', float64),               # a simple scalar field
    ('dz0', float64),          # an array field
]


@jitclass(spec)
class info:

    
    def __init__(self):
        """ Init INFO Class"""


        self.dzdeep = 2.0
        self.dz0 = 0.1

    def read_info(self,file):
        with open(file) as f:
            lines = f.readlines()
            
            # global lcomment
            
        self.lcomment = int(lines[0].split()[0])
        self.lerrorgap = int(lines[1].split()[0])
        self.lmc = int(lines[2].split()[0])
        self.lhourlysnowoutput = int(lines[3].split()[0])
        self.ibyear = int(lines[4].split()[0])
        self.ilyear = int(lines[5].split()[0])
        self.tstep = int(lines[6].split()[0])
        self.mbsumdy, mbwindy= float(lines[7].split()[0]), float(lines[7].split()[1])
        self.dz0 = float(lines[8].split()[0])
        self.dzdeep = float(lines[9].split()[0])
        self.zdeep = float(lines[10].split()[0])
        self.densice, self.densfirn= float(lines[11].split()[0]), float(lines[11].split()[1])
        self.rhosn, self.rhosnprec, self.trhoprec = float(lines[12].split()[0]), float(lines[12].split()[1]), int(lines[12].split()[2])
        self.tpcond = int(lines[13].split()[0])
        self.T10m = float(lines[14].split()[0])
        self.dsnow, self.dfirn  = float(lines[15].split()[0]), float(lines[15].split()[1])
        self.luseacc = int(lines[16].split()[0])
        self.tpdens, self.lrefr = int(lines[17].split()[0]), int(lines[17].split()[1])
        self.tpirre, self.cirre = int(lines[18].split()[0]), float(lines[18].split()[1])
        self.lslush = int(lines[19].split()[0])
        self.surfangle, self.tausteep, self.tauhor,self.tau1,self.slfact = float(lines[20].split()[0]), float(lines[20].split()[1]), float(lines[20].split()[2]), float(lines[20].split()[3]), float(lines[20].split()[4])
        self.accyear = float(lines[21].split()[0])
        self.lz0m, self.zll, self.zul, self.Hmax = float(lines[22].split()[0]), float(lines[22].split()[1]), float(lines[22].split()[2]), float(lines[22].split()[3])
        self.z0msn, self.z0mice = float(lines[23].split()[0]), float(lines[23].split()[1])
        self.lz0h = int(lines[24].split()[0])
        self.tcalc = int(lines[25].split()[0])
        self.extrapolation = int(lines[26].split()[0])
        self.lsnet = int(lines[27].split()[0])
        self.albmin, self.albmax = float(lines[28].split()[0]), float(lines[28].split()[1])
        self.emis = float(lines[29].split()[0])
        self.lalbedo, self.solzenyes, self.SSAfresh, self.radrefr = int(lines[30].split()[0]), int(lines[30].split()[1]), float(lines[30].split()[2]), float(lines[30].split()[3])
        self.albsnow, self.albice, self.albfirn, self.soot = float(lines[31].split()[0]), float(lines[31].split()[1]), float(lines[31].split()[2]), float(lines[31].split()[3])
        self.snowstar,self.tstarwet,self.tstardry0,self.tstardry10 = float(lines[32].split()[0]), float(lines[32].split()[1]), float(lines[32].split()[2]), float(lines[32].split()[3])
        self.penetration = int(lines[33].split()[0])
        self.dzrad,self.zradmax = float(lines[34].split()[0].replace(",","")), float(lines[34].split()[1])
        self.radiussn, self.radiusice = float(lines[35].split()[0]), float(lines[35].split()[1])
        self.lwcloud = int(lines[36].split()[0])
        self.lwmax = float(lines[37].split()[0]), float(lines[37].split()[1]), float(lines[37].split()[2])
        self.lwmin = float(lines[38].split()[0]), float(lines[38].split()[1]), float(lines[38].split()[2])
        self.depthin = float(lines[39].split()[0]), float(lines[39].split()[1]), float(lines[39].split()[2]), float(lines[39].split()[3]), float(lines[39].split()[4])
        self.lclimtemp,self.climtemp = int(lines[40].split()[0]), float(lines[40].split()[1])
        self.lclimprec,self.climprec = int(lines[41].split()[0]), float(lines[41].split()[1])
        self.lclimrad,self.climrad = int(lines[42].split()[0]), float(lines[42].split()[1])
        self.lclimws,self.climws = int(lines[43].split()[0]), float(lines[43].split()[1])

        
