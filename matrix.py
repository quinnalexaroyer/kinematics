class Matrix:
	def __init__(self,*a):
		if len(a) == 0:
			self.m = []
		elif len(a) == 2 and isinstance(a[0],int) and isinstance(a[1],int):
			self.m = [[0 for i in range(a[1])] for j in range(a[0])]
		else:
			self.m = [[a[i][j] for j in range(len(a[i]))] for i in range(len(a))]
	def __add__(self,other):
		return Matrix(*[[self.m[i][j] + other.m[i][j] for i in range(len(self.m))] for j in range(len(self.m[0]))])
	def __mul__(self,other):
		return Matrix(*[[self.dot(self.m[i],other.col(j)) for j in range(len(other.m[0]))] for i in range(len(self.m))])
	def cp(self):
		return Matrix(*self.m)
	def col(self,n):
		return [self.m[i][n] for i in range(len(self.m))]
	def dot(self,x,y):
		sum = 0
		for i in range(len(x)):
			sum += x[i]*y[i]
		return sum
	def t(self):
		return Matrix(*[[self.m[i][j] for i in range(len(self.m))] for j in range(len(self.m[0]))])
	def tt(self):
		return Matrix(*[[self.dot(self.m[i],self.m[j]) for j in range(len(self.m))] for i in range(len(self.m))])
	def ttt(self):
		return self.t().tt()
	def rplus(self,row,otherRow,factor):
		for i in range(len(self.m[row])):
			self.m[row][i] += factor*self.m[otherRow][i]
			self.m[row][i] = round(self.m[row][i],7)
	def rtimes(self,row,factor):
		for i in range(len(self.m[row])):
			self.m[row][i] *= factor
			self.m[row][i] = round(self.m[row][i],7)
	def rswap(self,row1,row2):
		temp = self.m[row1]
		self.m[row1] = self.m[row2]
		self.m[row2] = temp
	def ref(self):
		pivot = [0,0]
		while pivot[0] < len(self.m):
			foundPivot = False
			while not foundPivot and pivot[1] < len(self.m[pivot[0]]):
				i = pivot[0]
				while not foundPivot and i < len(self.m):
					if self.m[i][pivot[1]] != 0:
						foundPivot = True
					else:
						i += 1
				if not foundPivot: pivot[1] += 1
			if pivot[1] < len(self.m[pivot[0]]):
				if self.m[pivot[0]][pivot[1]] == 0:
					self.rplus(pivot[0],i,1)
				for i in range(pivot[0]+1,len(self.m),1):
					self.rplus(i,pivot[0],-self.m[i][pivot[1]]/self.m[pivot[0]][pivot[1]])
			else:
				pivot[0] = len(self.m)
			pivot[0] += 1
	def rref(self):
		self.ref()
		for i in range(len(self.m)):
			j = i
			while j < len(self.m[i]) and self.m[i][j] == 0:
				j += 1
			if j < len(self.m[i]):
				for k in range(i):
					self.rplus(k,i,-self.m[k][j]/self.m[i][j])
				self.rtimes(i,1/self.m[i][j])
	def det(self):
		c = self.cp()
		c.ref()
		d = 1
		for i in range(len(c.m)):
			d *= c.m[i][i]
		return d
	def inv(self):
		v = Matrix(*[self.m[i]+[1 if i == j else 0 for j in range(len(self.m[i]))] for i in range(len(self.m))])
		v.rref()
		return Matrix(*[[v.m[i][j] for j in range(len(v.m[i])//2,len(v.m[i]),1)] for i in range(len(v.m))])
	def pseudoInvL(self):
		mmt = self.tt()
		if mmt.det() != 0:
			return self.ttt().inv()*self.t()
		else:
			print("DET = 0")
	def pseudoInvR(self):
		mmt = self.tt()
		if mmt.det() != 0:
			return self.t()*self.tt().inv()
		else:
			print("DET = 0")
