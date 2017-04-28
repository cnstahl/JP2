from syk_susy_2 import *
import manipulate

H = hamiltonian(9,1)
print("Hamiltonian built")
w,v = manipulate.eigen(H)
print(len(w), " total states")
print(sum(i< 1e-10 for i in w), " ground states")

mid = 0
for idx, val in enumerate(w):
	if val > 5e-0 and val < 5e-9:
		print(idx, "'th eval is", val)
		mid = 1
		break
if mid == 1:
	for idx, val in enumerate(w):
		if val > 5e-9:
			print(idx, "'th eval is", val)
			break