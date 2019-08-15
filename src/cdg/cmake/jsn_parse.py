########################################################################
# Python script to parse primary JSON file for 
# info
########################################################################

import json
import sys
from pprint import pprint

jdata = open(sys.argv[1])

data = json.load(jdata)

...

jdata.close()
