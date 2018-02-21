from pylab import *
from prody import *
import Bio.PDB as pdb
import numpy as np;
import matplotlib.pyplot as plt


#tested on this template
#only concerned main chain
PDB_FILE='5fs4.pdb';
def histogram(pdb_file):
	backbone = parsePDB(pdb_file, subset='bb');
	coordinates = backbone.getCoords()
	names = backbone.getNames();
	numCA = int(names.size/4);
	caDistances =[];
	psi =[];
	for i in range(0, numCA):
		for j in range(i, numCA):
			indexI= 4*i+1
			indexJ= 4*j+1;
			CAi = coordinates[indexI, 0:3];
			CAj = coordinates[indexJ, 0:3];
			currentDistance = np.linalg.norm(CAi-CAj);
			caDistances.append(currentDistance);

	plt.hist(caDistances,50);
	plt.xlabel("CA distances");
	plt.ylabel("Frequency");
	plt.title("Histogram of CA distances");
	plt.show()

histogram(PDB_FILE);
histogram("2b7h.pdb");
