# -*- coding: utf-8 -*-

from __future__ import print_function, division

import sys
import os
from os import path
import scipy as sp
import time
import subprocess

class Noise(object):

    def __init__(self):
        self.params = {'SEEING': '0.80',
                       'ELEV': '4000',
                       'ZA': '45.00',
                       'LUNARPHASE': '0',
                       'LUNARANGLE': '60.0',
                       'LUNARZA': '30.0',
                       'EBV': '0.00',
                       'FIELDANG': '0.675',
                       'DECENT': '0.0',
                       'T_EXP': '450',
                       'N_EXT': '8',
                       'REFF': '0.3',
                       }
        self.spectrograph_config = 'config/PFS.20151204.dat'
        self.output_noise_file = 'out/ref.noise.dat'
        self.overwrite = False

        self.HOME_DIR = path.dirname(path.abspath(__file__))
        self.NOISE_SRC = 'bin/gsnoise.x'
        if not os.path.exists('out'):
            os.mkdir('out')
        if not os.path.exists(self.ETC_SRC):
            exit("Unable to find ETC engine; please run make first and try again")
        return None

    def run(self):
        start = time.time()
        ''' run NOISE '''
        try:
            print('##### starting to run NOISE ... (it takes a few min.) #####')
            args = [
                os.path.join(self.HOME_DIR, self.NOISE_SRC),
                os.path.join(self.HOME_DIR, self.spectrograph_config),
                '-', # observation config is taken from stdin
                self.output_noise_file,
            ]
            proc = subprocess.Popen(args, stdin=subprocess.PIPE)
            obs = ""
            for k in self.params:
                obs += '{} {}\n'.format(k, self.params[k])
            proc.communicate(obs.encode())
        except OSError as e:
            exit('Execution error of "%s" (%s)' % self.NOISE_SRC, e)
            
        elapsed_time = time.time() - start
        print("##### finished (elapsed_time: %.1f[sec]) #####" % (elapsed_time))

        return 0
