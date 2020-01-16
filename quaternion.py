import math

class Quaternion:
	round = 0
	def __init__(self,*q):
		if len(q) == 4:
			self.w = q[0]
			self.x = q[1]
			self.y = q[2]
			self.z = q[3]
		else:
			self.w = 1
			self.x = 0
			self.y = 0
			self.z = 0
	def __add__(self,other):
		if isinstance(other,(int,float,complex)):
			return Quaternion(self.w+other,self.x,self.y,self.z)
		elif isinstance(other,Quaternion):
			return Quaternion(self.w+other.w,self.x+other.x,self.y+other.y,self.z+other.z)
	def __mul__(self,other):
		if isinstance(other,(int,float,complex)):
			return Quaternion(other*self.w,other*self.x,other*self.y,other*self.z)
		if isinstance(other,Quaternion):
			w = self.w*other.w - self.x*other.x - self.y*other.y - self.z*other.z
			x = self.w*other.x + self.x*other.w + self.y*other.z - self.z*other.y
			y = self.w*other.y - self.x*other.z + self.y*other.w + self.z*other.x
			z = self.w*other.z + self.x*other.y - self.y*other.x + self.z*other.w
			return Quaternion(w,x,y,z)
	def __invert__(self):
		return self.conj()/self.norm()**2
	def __truediv__(self,other):
		if isinstance(other,(int,float,complex)):
			return Quaternion(self.w/other,self.x/other,self.y/other,self.z/other)
		else:
			return self*(~other)
	def __neg__(self):
		return Quaternion(-self.w,-self.x,-self.y,-self.z)
	def __sub__(self,other):
		return self+(-other)
	def __eq__(self,other):
		return self.w == other.w and self.x == other.x and self.y == other.y and self.z == other.z
	def __neq__(self,other):
		return not(self == other)
	def __repr__(self):
		if self.round == 0:
			return "<Quaternion %s %si %sj %sk>" % (self.w,self.x,self.y,self.z)
		else:
			s = "<Quaternion %%.%if %%.%ifi %%.%ifj %%.%ifk>" % (self.round,self.round,self.round,self.round)
			return s % (self.w,self.x,self.y,self.z)
	def copy(self):
		return Quaternion(self.w,self.x,self.y,self.z)
	def norm(self):
		return (self.w**2+self.x**2+self.y**2+self.z**2)**0.5
	def conj(self):
		return Quaternion(self.w,-self.x,-self.y,-self.z)
	def normalize(self):
		n = self.norm()
		self.w /= n
		self.x /= n
		self.y /= n
		self.z /= n
	def getAxisAngle(self):
		t = 2*math.acos(self.w)
		if math.sin(t/2) != 0:
			v = tuple([i/math.sin(t/2) for i in (self.x,self.y,self.z)])
		else:
			v = (0,0,0)
		return (t,v)
	def getAxisAngle1(self):
		t = math.acos(self.w)
		v = tuple([i/math.sin(t) for i in (self.x,self.y,self.z)])
		return (t,v)
	def qrot(self,other):
		q = other.copy()
		q.normalize()
		return q*self*(~q)
	def rot(self,angle,axis):
		q = Quaternion.fromAxisAngle(angle,axis)
		return q*self*(~q)
	def dot(self,other):
		return self.w*other.w + self.x*other.x + self.y*other.y + self.z*other.z
	def cross(self,other):
		q = self*other
		return Quaternion(0,q.x,q.y,q.z)
	def dist(self,other):
		return (other-self).norm()
	def dc(self,other):
		pr = self*other
		angle = math.acos(round(-pr.w,6))
		axis = tuple([h/math.sin(angle) for h in (pr.x,pr.y,pr.z)])
		return (angle,axis)
	def euler1(self,axis):
		l = (self.w,self.x,self.y,self.z)
		(s,t,u,v) = (l[0],l[axis+1],l[((axis+1)%3)+1],l[((axis+2)%3)+1])
		return math.atan2(2*(s*t+u*v),1-2*(t**2+u**2))
	def axisRot(self,axis):
		(x,y,z) = (Quaternion(0,1,0,0),Quaternion(0,1,0,0),Quaternion(0,1,0,0))
		x1 = self*x/self
		y1 = self*y/self
		z1 = self*z/self
	def rotCo(self):
		y = Quaternion(0,0,1,0)
		z = Quaternion(0,0,0,1)
		p = self*y/self
		pxy = Quaternion(0,p.x,p.y,0)
		pxy.normalize()
		rz = -math.atan2(pxy.x,pxy.y)
		rx = math.atan2(p.z,(p.x**2+p.y**2)**0.5)
		s = Quaternion.aa(rz,(0,0,1))*Quaternion.aa(rx,(1,0,0))
		z1 = self*z/self
		z2 = s*z/s
		zp = z1*z2
		ry = math.acos(round(-zp.w,6))
		if p.y != 0 and zp.y/p.y > 0:
			ry = -ry
		return (rx,ry,rz)
	def inLimits(self,limits):
		rc = self.rotCo()
		if limits[0] is not None and r2d(rc[0]) < limits[0]: return False
		if limits[1] is not None and r2d(rc[0]) > limits[1]: return False
		if limits[2] is not None and r2d(rc[1]) < limits[2]: return False
		if limits[3] is not None and r2d(rc[1]) > limits[3]: return False
		if limits[4] is not None and r2d(rc[2]) < limits[4]: return False
		if limits[5] is not None and r2d(rc[2]) > limits[5]: return False
		return True
	def co(self,axis):
		return (self.w,self.x,self.y,self.z)[axis]
	@staticmethod
	def setRound(n):
		Quaternion.round = n
	@staticmethod
	def fromAxisAngle(t,v):
		if isinstance(v,Quaternion):
			v = (v.x,v.y,v.z)
		w = math.cos(t/2)
		n = (v[0]**2+v[1]**2+v[2]**2)**0.5
		x = v[0]*math.sin(t/2)/n
		y = v[1]*math.sin(t/2)/n
		z = v[2]*math.sin(t/2)/n
		return Quaternion(w,x,y,z)
	@staticmethod
	def fromAxisAngle1(t,v):
		w = math.cos(t)
		n = (v[0]**2+v[1]**2+v[2]**2)**0.5
		x = v[0]*math.sin(t)/n
		y = v[1]*math.sin(t)/n
		z = v[2]*math.sin(t)/n
		return Quaternion(w,x,y,z)
	@staticmethod
	def aa(t,v):
		return Quaternion.fromAxisAngle(t,v)
	def aaget(self):
		return self.getAxisAngle()

Q = Quaternion

def degrees(t):
	return t*math.pi/180

def getDegrees(t):
	if t is not None:
		return t*180/math.pi
	else:
		return None

d = degrees
r2d = getDegrees
