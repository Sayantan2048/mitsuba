# NOTE: remember to specify paths using FORWARD slashes (i.e. '/' instead of
# '\' to avoid pitfalls with string escaping)
# Configure the search path for the Python extension module

import sys, os

# mtsdir = os.environ['MITSUBADIR']
mtsdir = 'D:/projects/mitsuba/dist'
#mtsdir = 'C:/Users/sherholz/Develop/guiding-jirka-merge/Mitsuba0.5.0/dist'
if sys.platform.startswith('win'):
    sys.path.append(mtsdir + '/python/3.6')
else:
    sys.path.append(mtsdir + '/python/3.6')
os.environ['PATH'] = mtsdir + os.pathsep + os.environ['PATH']

import mitsuba
import  mitsuba.core as mtsCore
import mitsuba.render as mtsRender
import mitsuba.importance as mtsImp
#from mitsuba.render import SceneHandler


class MyFormatter(mtsCore.Formatter):
    def format(self, logLevel, sourceClass, sourceThread, message, filename, line):
        #return '%s (log level: %s, thread: %s, class %s, file %s, line %i)'(message, str(logLevel), sourceThread.getName(), sourceClass, filename, line)
        return "LOG: " + message
class MyAppender(mtsCore.Appender):
    def append(self, logLevel, message):
        print(message)
    def logProgress(self, progress, name, formatted, eta):
        print('Progress message: ' + formatted)
