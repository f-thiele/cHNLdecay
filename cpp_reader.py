from subprocess import Popen, PIPE

output = Popen(['./cHNLdecay', '--generations', '13', '--mass', '20000', '--angle', '3.5e-6'], stdout=PIPE)
out = output.stdout.read()

print "out: ", out
out_pw = out.split()[0];

print out_pw;
