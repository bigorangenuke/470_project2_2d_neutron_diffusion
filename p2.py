import rx
import numpy as np
import matplotlib.pyplot as plt
import time
dbg = False
ddbg = False

K_TOLERANCE = 1e-4

m = 60.
M_NODES = m
N_NODES = m

REACTOR_WIDTH	= 230.
REACTOR_HEIGHT	= 230.


class ReactorSolver():
	
	############################################
	#Solves 2D reactor iteratively via finite difference
	#A\phi=\frac{1}[k}B\phi = S(\phi)
	
	def __init__(self,reactor):
		if dbg: print("init ReactorSolver()")
		self.rx = reactor
		
		#Initial guesses, ones across the board
		self.phi = np.ones(((M_NODES-2)*(N_NODES-2),1))
		self.source = np.eye((M_NODES-2)*(N_NODES-2))
		self.k = 1.
		
	def A(self):
		#empty a array
		a=big_zeros()
		empty_slice = np.zeros((M_NODES-2,N_NODES-2))
		
		#iterate through indicies of this empty array
		#assign the new slice to a slice of the big array
		for (i,j),s in np.ndenumerate(empty_slice):
			a[:,i*(N_NODES-2)+j] = self.get_slice(i,j,empty_slice)
		return a
		
	def B(self):
		#empty
		b=big_zeros()
		
		#slice off the ends because they are zero
		nds = self.rx.nodes[1:-1,1:-1]
		
		#assign the production term to the diagonal of b
		for (i,j),nd in np.ndenumerate(nds):
			b[(i*(N_NODES-2))+j,(i*(N_NODES-2))+j] = nd.NuSigma_f()*self.rx.dx*self.rx.dy
		return b
		
	def get_slice(self,i,j,slice_like):
		dx = self.rx.dx
		dy = self.rx.dy
		
		tmp =np.zeros_like(slice_like)
		#tmp = np.zeros((M_NODES -2) *(N_NODES -2))
		tmp[i,j] = 2*self.rx.nodes[i,j].D()*((dx/dy)+(dy/dx))+(self.rx.nodes[i,j].Sigma_a()*dx*dy)
		
		if i:
			tmp[i-1,j] = -self.rx.nodes[i-1,j].D()*dy/dx
		if j:
			tmp[i,j-1] = -self.rx.nodes[i,j-1].D()*dx/dy
		if j<tmp.shape[1]-1:
			tmp[i,j+1] = -self.rx.nodes[i,j+1].D()*dx/dy
		if i<tmp.shape[0]-1:
			tmp[i+1,j] = -self.rx.nodes[i+1,j].D()*dy/dx
		
		#squeeze 2d array to 1d array
		return tmp.reshape((tmp.shape[0]*tmp.shape[1]))

		
	def solve(self):
		
		itr = 0
		
		#
		oldsource = 1./self.k*np.matrix(self.B())*self.phi
		#print('oldsource',oldsource.shape)
		done = False
		while not done:
			itr+=1
			#print('==Iteration %s=='%(itr))
			#===================================================================
			# print('i=%s k=%s'%(itr,self.k))
			#===================================================================
			b = np.matrix(self.B())
			ai = np.linalg.inv(np.matrix(self.A()))
			p = self.phi
			#print('a shape = %s \nb shape = %s\nphi shape = %s'%(ai.shape,b.shape,p.shape))
			#print('phi shape = %s \nsource shape = %s'%(source.))
			newphi = ai*oldsource
			newsource = 1./self.k*b*newphi
			newk = self.k*np.sum(newsource)/np.sum(oldsource)
			
			self.phi = newphi
			oldsource = newsource
			if np.abs(self.k-newk)<K_TOLERANCE:
				print('size = %s,k = %s'%(self.rx.w, newk))
				done = True
			self.k = newk
		
		return self.k,self.phi
			
			
def minimum_critical_geometry():
	#Start close to answer
	rxsz = 210
	done = False
	
	x = []
	f = []
	while not done:
		rxsz +=1
		reactor = rx.Reactor(nodes = [M_NODES,N_NODES],size = [rxsz,rxsz])
		rxsolver = ReactorSolver(reactor)
		k,phi = rxsolver.solve()
		print(k)
		x.append(rxsz)
		f.append(k)
		if k>1:
			done = True
			
	print('The minimum critical geometry is %s'%(rxsz))
	return x,f
		
def k_vs_dim(maxDim = 1000,steps = 50):
	#Get k as a function of square reactor dimensions	
	x = []
	f = []
	
	dims = np.logspace(1,np.log10(maxDim),steps)
	
	for dim in dims:
		reactor = rx.Reactor(nodes = [M_NODES,N_NODES],size = [dim,dim])
		rxsolver = ReactorSolver(reactor)
		k = rxsolver.solve()
		x.append(dim)
		f.append(k)

	return x,f
		
	
	
		
def big_zeros():
	#empty zero array 
	return np.zeros(((M_NODES-2)*(N_NODES-2),(M_NODES-2)*(N_NODES-2)))					
		

def solve():
	#Reactor dimensions
	rxsize = [REACTOR_WIDTH,REACTOR_HEIGHT]

	#Create reactor
	reactor = rx.Reactor(nodes = [M_NODES,N_NODES],size = rxsize)

	#Solve
	rxsolver = ReactorSolver(reactor)
	return rxsolver.solve()
	
def timing(start,stop,step):
	ti=[]
	t = []
	for i in range (start,stop,step):
		print(i)
		M_NODES = float(i)
		N_NODES = float(i)
		
		starttime = time.clock()
		solve()
		t.append(time.clock()-starttime)
		ti.append(i)
		print(t[-1])
	
	return ti,t
	
if __name__ == "__main__":
	#===========================================================================
	# i,t  = timing(5,60,5)
	# plt.xlabel('Matrix dimension')
	# plt.ylabel('Time for convergence (s)')
	# plt.plot(i,t)
	# plt.show()
	#===========================================================================
	
	
	
	starttime = time.clock()
	solve()
	print(M_NODES,time.clock()-starttime)

	
	
	
	
	#===========================================================================
	# newt = np.asarray(t)/np.min(t)
	# plt.plot(i,newt)
	# plt.show()
	#===========================================================================
	#===========================================================================
	# x,f = minimum_critical_geometry()
	# solve()
	#===========================================================================
	#x,f=minimum_critical_geometry()
#===============================================================================
# 	k,phi = solve()
# 	print(phi.shape)
# 	p = phi.reshape((M_NODES-2,N_NODES-2))
# 	
# 	#tabulate position using node number for plotting
# 	x = np.zeros(M_NODES-2)
# 	x *= REACTOR_WIDTH/M_NODES
# 	y = np.zeros(N_NODES-2)
# 	y *= REACTOR_HEIGHT/N_NODES
# 
# 	#normalize
# 	p /=np.max(p)
# 	
# 	ax = plt.imshow(p,extent=[1,REACTOR_WIDTH,1,REACTOR_HEIGHT])
# 	plt.xlabel('x position (cm)')
# 	plt.ylabel('y position (cm)')
# 	cb = plt.colorbar()
# 	plt.title('Flux vs. position in the core')
# 	plt.show()
#===============================================================================



	
# 	plt.plot(x,y,p)
# 	plt.xlabel('x position')
# 	plt.ylabel('y position')
# 	cb = plt.colorbar()
# 	cb.ax.set_ylabel('flux (n/cm^2/s)')
# 	plt.title('Flux vs position in the core')
# 	plt.show()




# 	x,f = k_vs_dim()
# 	plt.xlabel('Reactor Size')
# 	plt.ylabel('k')
# 	plt.title('Multiplication factor vs reactor size')
# 	plt.plot()
# 	plt.show()
	
	#solve()
	

	


	
	
	




