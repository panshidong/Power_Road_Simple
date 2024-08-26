import subprocess
def run_tapb(netfile,tripfile,gap='0.0001',classes='1'):
    """
    Runs tap-b and generate results
    Call rule: popen=subprocess.run(['tap-b/bin/tap','0.0001','1',NETFILE,TRIPFILE])
    It would generate a "flow.txt" under the run folder
    """
    popen=subprocess.run(['tap-b/bin/tap',gap,classes,netfile,tripfile],stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    #popen=subprocess.run(['tap-bak/bin/tap',gap,classes,netfile,tripfile])