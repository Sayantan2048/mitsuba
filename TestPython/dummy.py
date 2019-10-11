import os, sys
mitdir = 'D:/projects/mitsuba/dist'
sys.path.append('D:/projects/mitsuba/dist/python/3.6')
os.environ['PATH'] = mitdir + os.pathsep + os.environ['PATH']

import mitsuba

print("Hello")