from syk_susy_2 import *
import matplotlib.pyplot as plt

N=10 
H= hamiltonian(N, 1)

plt.imshow(np.log(abs(H)), cmap='YlOrRd', interpolation='nearest')

plt.title("Hamiltonian for $N=10$", fontsize=20, family='serif')
#plt.savefig("../data/hamil10.pdf")
plt.show()

N=5
Ha= hamiltonian(N, 1)

plt.imshow(np.log(abs(Ha)), cmap='YlOrRd', interpolation='nearest')

plt.title("Hamiltonian for $N=%s$" %N, fontsize=20, family='serif')
#plt.savefig("../data/hamil%s.pdf" % N)
plt.show()
