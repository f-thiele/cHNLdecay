from subprocess import Popen, PIPE

def get_prod_BR(HNLmassinMeV, HNLcoupling, B_ID, meson_ID, lepton_ID):
	
	# 2 body (leptonic) decay
	if(meson_ID==0): 
		output = Popen(['./cHNLdecay', '--mainmode', '1', '--BmesonID', str(B_ID), '--generations', str(lepton_ID), '--mass', str(HNLmassinMeV), '--angle', str(HNLcoupling)], stdout=PIPE)
	
	# 3 body (semileptonic) decay
	else:
		output = Popen(['./cHNLdecay', '--mainmode', '1', '--BmesonID', str(B_ID), '--generations', str(lepton_ID), '--mass', str(HNLmassinMeV), '--angle', str(HNLcoupling)], stdout=PIPE)
	
	
	out = output.stdout.read()
	out_BR = float(out.split()[0]);
	
	return out_BR;
