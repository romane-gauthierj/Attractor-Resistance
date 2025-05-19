#!/Users/romane/repos/icr/attractors_resistance/res_attr/bin/python3.12
from __future__ import print_function
import os, sys, maboss

def main(argv):

	if len(argv) < 3 or len(argv) > 6:
		print(
			("Bad arguments, MBSS_FormatTable.py <file.bnd> <file.cfg> "
			"<(optional)threshold> <(optional)-mb MaBoSS executable name>"), 
			file=sys.stderr
		)
		exit()

	bndFile = argv[1]
	cfgFile = argv[2]
	mabossCommand = None
	probCutoff = None
	
	if len(argv) > 3:
		index = 3
		while index < len(argv):
			if argv[index].startswith("-mb"):
				mabossCommand = argv[index+1]
				index += 2
			else:
				probCutoff = float(argv[index])
				index += 1

	simulationName = os.path.splitext(os.path.basename(argv[2]))[0]
	
	simulation = maboss.load(bndFile, cfgFile, command=mabossCommand)
	result = simulation.run(command=mabossCommand, prefix=simulationName)
	
	result.save(simulationName)
	
	table = result.get_states_probtraj_full()
	table.to_csv(os.path.join(simulationName, ("%s_probtraj_table.csv" % simulationName)), sep="\t")
	
	result.write_statdist_table(
		os.path.join(simulationName, "%s_statdist_table.csv" % simulationName), 
		probCutoff
	)
	
if __name__ == '__main__':
    main(sys.argv)