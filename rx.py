import numpy as np

MINIMUM_REACTOR_DIMENSION = 4

dbg = False
dbg_itr = False
class Material():
	def __init__(self,material):
		lines = material.split(',')
		lines=list((line.strip() for line in lines))
		#print(lines)
		self.name=					lines[0]
		self.Sigma_tr = 			lines[1]
		self.Sigma_a = 				lines[2]
		self.nuSigma_f=				lines[3]
		self.relative_absorption=	lines[4]
		self.atomic_mass = 			lines[5]

	def L(self):
		#Extrapolation Length
		s_tr = self.Sigma_tr
		s_a= self.Sigma_a
		return 1/np.sqrt(3*s_tr*s_a)*(1+0.4*s_a/s_tr)		

	def D(self):
		#Diffusion coefficient
		return 1./(3.*self.Sigma_tr())
	
	def __repr__(self):
		return self.name




class Node():
	def __init__(self,i,j,materials=None):
		self.i = i
		self.j = j
		self.materials = materials
		#print (materials)
		self.sigma_tr = self.get_sigma_tr()
		self.sigma_a = self.get_sigma_a()
		self.nuSigma_f = self.get_nuSigma_f()
		
		self.source = 1.
		self.phi = 1.
	def __repr__(self):
		return "(%s,%s)"%(self.i,self.j)
	
	def get_sigma_tr(self):
		#Get array of the Sigma_tr of all the materials
		return np.array(list(material.Sigma_tr for material in self.materials),dtype="float")

	def get_nuSigma_f(self):
		return np.array(list(material.nuSigma_f for material in self.materials),dtype='float')
		
	def get_sigma_a(self):
		#Get array of the Sigma_a of all the materials
		return np.array(list(material.Sigma_a for material in self.materials),dtype='float')

	def Sigma_a(self):
		#Macroscopic absorption cross section
		return np.sum(self.sigma_a)
	
	def NuSigma_f(self):
		return np.sum(self.nuSigma_f)
	
	def Sigma_tr(self):
		#Macroscopic transient cross section
		return np.sum(self.sigma_tr)
	
	def L(self):
		#Extrapolation Length
		s_tr = self.Sigma_tr
		s_a= self.Sigma_a
		return 1/np.sqrt(3*s_tr*s_a)*(1+0.4*s_a/s_tr)		

	def D(self):
		#Diffusion coefficient
		return 1./(3.*self.Sigma_tr())



class Reactor():
	def __init__(self,**kwargs):
		##############################
		#default values
		#Nodes in x
		self.m = 5
		#Nodes in y
		self.n = 5
		#Physical height
		self.h= 1.
		#Physical width
		self.w = 1.
		#Step size in x
		self.dx = 1.
		#Step size in y
		self.dy = 1.
		###############################
	
		self.set_reactorSize(size=[self.w,self.h])
		
		#Nodes in the reactor
		if "nodes" in kwargs:
			nds = kwargs["nodes"]
			if not any(nd<MINIMUM_REACTOR_DIMENSION for nd in nds):
				self.m = nds[1]
				self.n = nds[0]
			else: print ("BAD TROUBLE. One or more reactor node dimensions below MINIMUM_REACTOR_DIMENSION = %s"%(MINIMUM_REACTOR_DIMENSION))
		#Physical Dimensions of reactor
		if "size" in kwargs:
			self.set_reactorSize(kwargs["size"])

		self.materials = self.load_materials()


		thenodes = self.slice()

		for j in range(int(self.m)):
			for i in range(int(self.n)):
				thenodes[i,j] = Node(i,j,self.materials)
		self.nodes = thenodes
		
	def set_reactorSize(self,size):
		self.w = size[0]
		self.dx = self.w/(self.m-1)
		
		self.h = size[1]
		self.dy = self.h/(self.n-1)
		
	def slice(self,centralNode=None):
		if dbg: print("===Start rx.slice(%s)==="%(centralNode))
		
		#Initial instantiation of self.nodes
		if centralNode==None: return np.zeros((self.n,self.m),dtype=object)
		newSlice = np.zeros_like(self.nodes)
		
		
		
		for i in range(centralNode.i):
			for j in range(centralNode.j):
				print("Slice ij",i,j)
				newSlice[i,j]=self.nodes[i,j]
	    #print(newSlice)
		return newSlice
		
	def load_materials(self):
		#Pull out of text file and load line by line to materials
		f = open('macroscopiccrosssections.txt')
		lines = f.readlines()
		f.close()

		mats=[]
		for line in lines:
			mats.append(Material(line))
		return mats

	
	def getCoefficient(self,node=None,isCentralNode=False):
		d = 1
		s = 1
		if not node: print("Error. No node passed to rx.getCoefficient()")
		dx = self.dx
		dy = self.dy
		
		if node:
			if isCentralNode:
				d = node.D()
				s = node.Sigma_a()
				i = dx
				j = dy
			else:
				d = node.D()
		else: print('Warning. No node.')
			
		return (2*d*j+2*d*i+s*i*j)/(i*j)
	
	
	def __repr__(self):
		return "Rx(m = %s, n =  %s)"%(self.m,self.n)



if __name__=='__main__':
	
	#Node values
	m = 5
	n = 5
	nodes = [m,n]
	
	#Reactor dimensions
	w = 10. 
	h = 10.
	rxsize = [w,h]
	
	#Create reactor
	reactor = Reactor(nodes = nodes,size = rxsize)
	print(reactor)

	
	




