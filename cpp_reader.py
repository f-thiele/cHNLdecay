from subprocess import Popen, PIPE

output = Popen(['./cHNLdecay', '--mainmode', '1', '--BmesonID', '521', '--generations', '13', '--mass', '3500', '--angle', '3.5e-6'], stdout=PIPE)
out = output.stdout.read()

out_pw = out.split()[0];

print out_pw;
