from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import sys
from datetime import date

def oamCompute(delta_eta, eta_max, IC = "Born", alpha = 0.25):
	# IC = "Born" or IC = "Ones"
	n_steps = round(eta_max / delta_eta)
	logPrint("n steps ="+str(n_steps))
	if IC == "Ones":
		G3 = np.ones((n_steps + 1, n_steps + 1))
		G4 = np.ones((n_steps + 1, n_steps + 1))
		G5 = np.ones((n_steps + 1, n_steps + 1))
		G6 = np.ones((n_steps + 1, n_steps + 1))
		
		Gm3 = np.ones((n_steps + 1, n_steps + 1, n_steps + 1))
#         Gm4 = np.ones((n_steps + 1, n_steps + 1, n_steps + 1))
		Gm5 = np.ones((n_steps + 1, n_steps + 1, n_steps + 1))
		Gm6 = np.ones((n_steps + 1, n_steps + 1, n_steps + 1))
	else: 
		Nc = 3
		CF = ((Nc**2)-1)/(2*Nc)
		coef = ((alpha**2)*(CF**2) * np.pi)/(4*Nc)
		G3 = np.array([[coef for j in range(n_steps + 1)] for i in range(n_steps + 1)], dtype='longdouble')
		G4 = np.array([[0 for j in range(n_steps + 1)] for i in range(n_steps + 1)], dtype='longdouble')
		G5 = np.array([[0 for j in range(n_steps + 1)] for i in range(n_steps + 1)], dtype='longdouble')
		G6 = np.array([[0 for j in range(n_steps + 1)] for i in range(n_steps + 1)], dtype='longdouble')
		Gm3 = np.array([[[coef for j in range(n_steps + 1)] for k in range(n_steps + 1)] for i in range(n_steps + 1)], dtype='longdouble')
#         Gm4 = np.array([[[0 for j in range(n_steps + 1)] for k in range(n_steps + 1)] for i in range(n_steps + 1)])
		Gm5 = np.array([[[0 for j in range(n_steps + 1)] for k in range(n_steps + 1)] for i in range(n_steps + 1)], dtype='longdouble')
		Gm6 = np.array([[[0 for j in range(n_steps + 1)] for k in range(n_steps + 1)] for i in range(n_steps + 1)], dtype='longdouble')
		
	ic = 0 # one step difference between initial conditions
	d = (delta_eta ** 2)
	for j in range(1, n_steps + 1):
		for i in range(j + 1):
			
			G3[i, j] = G3[i, j-1] + ic
			G4[i, j] = G4[i, j-1] + ic
			G5[i, j] = G5[i, j-1] + ic
			G6[i, j] = G6[i, j-1] + ic
			
			Gm3[i, i,j] = G3[i, j]
#             Gm4[i, i,j] = G4[i, j]
			Gm5[i, i,j] = G5[i, j]
			Gm6[i, i,j] = G6[i, j]
			
			# G IR loop
			for ii in range(i):
				G3[i, j] += d * (2*G3[ii, j-i+ii] + G5[ii, j-i+ii] + 3*G6[ii,j-i+ii])
				G4[i, j] += 0.5 * d * (2*G3[ii, j-i+ii] + 7*G4[ii,j-i+ii] + 7*G5[ii,j-i+ii] - 0.5*G6[ii,j-i+ii])
				G5[i, j] += 0.5 * d * (-2*G3[ii, j-i+ii] - 4*G4[ii,j-i+ii] - 11*G5[ii,j-i+ii] - 3*G6[ii,j-i+ii])
				G6[i, j] += 0.5 * d * (4*G3[ii, j-i+ii] + 2*G6[ii,j-i+ii])
				
			# G UV loop
			for ii in range(i, j):
				G3[i, j] += d * (Gm3[i,ii,j-1] + 0.5*Gm5[i,ii,j-1] + 3*Gm6[i,ii,j-1])
				
			# Gm UV loop 
			for k in range(i+1, j+1):
				Gm3[i, k, j] = Gm3[i, k-1, j-1] + ic
				for ii in range(k-1, j):
					Gm3[i, k, j] += d * (Gm3[i,ii,j-1] + 0.5*Gm5[i,ii,j-1] + 3*Gm6[i,ii,j-1])
			
			# Gm IR loop
			for k in range(i+1, j+1):
#                 Gm4[i, k, j] = Gm4[i, k-1, j-1] 
				Gm5[i, k, j] = Gm5[i, k-1, j-1] 
				Gm6[i, k, j] = Gm6[i, k-1, j-1] 
				
		if(np.mod(j*delta_eta, 1) == 0): logPrint(np.round(j * delta_eta, 4))
	
	logG3 = np.log(np.abs(G3) + 0.0000001)
	logG4 = np.log(np.abs(G4) + 0.0000001)
	logG5 = np.log(np.abs(G5) + 0.0000001)
	logG6 = np.log(np.abs(G6) + 0.0000001)
	logGm3 = np.log(np.abs(Gm3) + 0.0000001)
#     logGm4 = np.log(np.abs(Gm4) + 0.0000001)
	logGm5 = np.log(np.abs(Gm5) + 0.0000001)
	logGm6 = np.log(np.abs(Gm6) + 0.0000001)
	
	return logG3, logG4, logG5, logG6, logGm3, logGm5, logGm6 #, logGm4


def logPrint(text):
	original = sys.stdout
	f = open('data/logs/run2.log', 'a')
	sys.stdout = f
	print(text)
	sys.stdout = original
	f.close()


if __name__ == '__main__':
	# mdeltas = [10, 10, 20, 20, 30, 40, 50, 60, 60, 70]
	# deltas = [0.0125, 0.016, 0.025, 0.032, 0.0375, 0.05, 0.0625, 0.075, 0.08, 0.1]

	mdeltas = [50, 60, 60, 70]
	deltas = [0.0625, 0.075, 0.08, 0.1]

	logPrint('Started at '+str(date.today()))

	for mdelta, delta in zip(mdeltas, deltas):
		for ieta in range(10, mdelta+1, 10):
			logPrint("Starting (del, eta)="+str(delta)+" "+str(ieta))
			G3i, G4i, G5i, G6i, Gm3i, Gm5i, Gm6i = oamCompute(delta, ieta)
			np.savetxt("data/surface/G3_"+str(ieta)+"_"+str(delta)[2:]+".dat", G3i)
			np.savetxt("data/surface/G4_"+str(ieta)+"_"+str(delta)[2:]+".dat", G4i)
			np.savetxt("data/surface/G5_"+str(ieta)+"_"+str(delta)[2:]+".dat", G5i)
			np.savetxt("data/surface/G6_"+str(ieta)+"_"+str(delta)[2:]+".dat", G6i)
			logPrint("Finished (del, eta)="+str(delta)+" "+str(ieta))

