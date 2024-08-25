from sequence_utils import *
import subprocess
NETFILE="tap-b/net/SiouxFalls_net.txt"
TRIPFILE="tap-b/net/SiouxFalls_trips.txt"
net = create_network(NETFILE, TRIPFILE, 1, 1)
net.not_fixed = set([])
net.art_links = {}
net.maxruntime=str(10)
popen=subprocess.run(['tap-b/bin/tap','0.0001','1',NETFILE,TRIPFILE])