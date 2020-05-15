import matplotlib.pyplot as plt
import numpy as np

colors = {"matrixMultiplication":'m', "matrixMultiplicationForPosition":'y', "matrixMultiplicationTranspose":'c', "matrixMultiplicationStrassen":'g', "matrixMultiplicationStrassen_Opti":'k'}

# PLOT 
def plotDataFromFolder(folder, x_axis, title):
	algorithms = ["matrixMultiplication","matrixMultiplicationForPosition", "matrixMultiplicationTranspose", "matrixMultiplicationStrassen", "matrixMultiplicationStrassen_Opti"]
	alg_time = []

	for alg_name in algorithms:
		file = open(folder+alg_name+".txt","r") 
		data = file.readlines()
		file.close() 

		alg_time.append(np.array([float(x) for x in data[0].split(";")[:-1]], dtype = np.float32))


	for i in range(len(alg_time)):
		plt.plot(x_axis, alg_time[i], label=algorithms[i], color=colors[algorithms[i]])


	plt.ylabel("Tiempo [s]")
	plt.xlabel("Tamaño de la matriz")
	plt.legend()
	plt.grid()
	plt.title(title)

	plt.show()

plotDataFromFolder("C/", range(100,2901, 100), "Matrices Cuadradas")
plotDataFromFolder("C/nonSquared/", range(100, 1000, 100), "Matrices No Cuadradas")