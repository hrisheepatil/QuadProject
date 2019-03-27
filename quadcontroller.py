import numpy as np
import scipy
import sympy
import matplotlib as plt

class thrust:
    def __init__(self, motorcoeff, omg1,omg2,omg3,omg4):
        self.motorcoeff = motorcoeff
        self.omg1 = omg1
        self.omg2 = omg2
        self.omg3 = omg3
        self.omg4 = omg4

    def motorth1(self):
        return self.motorcoeff * (self.omg1) * (self.omg1)

    def motorth2(self):
        return self.motorcoeff * (self.omg2) * (self.omg2)

    def motorth3(self):
        return self.motorcoeff * (self.omg3) * (self.omg3)

    def motorth4(self):
        return self.motorcoeff * (self.omg4) * (self.omg4)

    def thrustmat(self):
        a =  np.array([self.motorcoeff * (self.omg1) * (self.omg1), self.motorcoeff * (self.omg2) * (self.omg2), self.motorcoeff * (self.omg3) * (self.omg3), self.motorcoeff * (self.omg4) * (self.omg4)])
        return(np.vstack(a))
    
class moment:
    def __init__(self, motorcoeff, omg1,omg2,omg3,omg4):
        self.motorcoeff = motorcoeff
        self.omg1 = omg1
        self.omg2 = omg2
        self.omg3 = omg3
        self.omg4 = omg4

    def motormo1(self):
        return self.motorcoeff * (self.omg1) * (self.omg1)

    def motormo2(self):
        return self.motorcoeff * (self.omg2) * (self.omg2)

    def motormo3(self):
        return self.motorcoeff * (self.omg3) * (self.omg3)

    def motormo4(self):
        return self.motorcoeff * (self.omg4) * (self.omg4)

    def momentmat(self):
       a = np.array([self.motorcoeff * (self.omg1) * (self.omg1), self.motorcoeff * (self.omg2) * (self.omg2), self.motorcoeff * (self.omg3) * (self.omg3), self.motorcoeff * (self.omg4) * (self.omg4)])
       return(np.vstack(a))

def bodytoworldrotmat(theta,phi,si):
    a = np.array([[1, (sin(theta)*sin(phi))/cos(phi), (cos(theta)*sin(phi))/cos(phi)],[0, cos(theta), -sin(theta)],[0, sin(theta)/cos(phi), cos(theta)/cos(phi)]])
    return a

def worldtobodyrotmat(x,y,z):
    a = np.array([[cos(y)*cos(z) - sin(x)*sin(y)*sin(z), -cos(y), cos(y)*cos(z) + cos(z)*cos(x)*cos(y)], [cos(z)*sin(y) + cos(y)*sin(x)*sin(z), cos(x)*cos(y), sin(y)*sin(z) - cos(y)*cos(z)*sin(x)], [-cos(z), sin(x), cos(x)*cos(z)]])
    return a

def eulerrate(p,q,r,theta,phi,si):
    a = np.vstack([p,q,r])
    return np.dot(bodytoworldrotmat(theta,phi,si),a)

class kalman:

    def __init__(self,Qang,Qbias,Rmeas):
        self.Qang = Qang #0.001
        self.Qbias = Qbias #0.003
        self.Rmeas = Rmeas #0.03
        self.p = anp.array([[0, 0], [0,0]])

    def acctoang(self, acc):
        self.accx = arcsin(acc[0]/9.81)
        self.accy = arctan(acc[1]/acc[2])
        
    def getanglex(self, gyro, freq):

        k = np.array([0, 0])

        dt = 1/freq
        rate = gyro - bias
        gyroang += dt * rate

        self.p[0,0] += dt * (dt*self.p[1,1] - self.p[0,1] - self.p[1,0] + self.Qang)
        self.p[0,1] -= dt * self.p[1,1]
        self.p[1,0] -= dt * self.p[1,1]
        self.p[1,1] += Qbias * dt

        S = self.p[0,0] + Rmeas
        k[0] = self.p[0,0] / S
        k[1] = self.p[1,0] / S

        y = self.accx - gyroang
        
        gyroang += k[0] * y
        bias += k[1] * y

        self.p[0,0] -= k[0] * self.p[0,0]
        self.p[0,1] -= k[0] * self.p[0,1]
        self.p[1,0] -= k[1] * self.p[0,0]
        self.p[1,1] -= k[1] * self.p[0,1]

        return gyroang

    def getangley(self, gyro, freq):

        k = np.array([0, 0])

        dt = 1/freq
        rate = gyro - bias
        gyroang += dt * rate

        self.p[0,0] += dt * (dt*self.p[1,1] - self.p[0,1] - self.p[1,0] + self.Qang)
        self.p[0,1] -= dt * self.p[1,1]
        self.p[1,0] -= dt * self.p[1,1]
        self.p[1,1] += Qbias * dt

        S = self.p[0,0] + Rmeas
        k[0] = self.p[0,0] / S
        k[1] = self.p[1,0] / S

        y = self.accy - gyroang
        
        gyroang += k[0] * y
        bias += k[1] * y

        self.p[0,0] -= k[0] * self.p[0,0]
        self.p[0,1] -= k[0] * self.p[0,1]
        self.p[1,0] -= k[1] * self.p[0,0]
        self.p[1,1] -= k[1] * self.p[0,1]

        return gyroang

class xm:
    def __init__(self,kp, ki, kd, xdes):
        self.kp = kp
        self.kd = kd
        self.ki = ki
        self.xdes = xdes

    def integ(self, xmes):
        self.q += (self.xdes - xmes)

    def x(self, xmes, xmesdot):
        x = self.kp*(xdes - xmes) + self.ki*self.q - self.kd*xmesdot
        return x

class ym:
    def __init__(self,kp, ki, kd, ydes):
        self.kp = kp
        self.kd = kd
        self.ki = ki
        self.ydes = ydes

    def integ(self, ymes):
        self.q += (self.ydes - ymes)

    def y(self, xmes, xmesdot):
        y = self.kp*(ydes - ymes) + self.ki*self.q - self.kd*ymesdot
        return y

class uthrust:
    def __init__(self,kp, ki, kd, zdes, mass):
        self.kp = kp
        self.kd = kd
        self.ki = ki
        self.zdes = zdes
        self.mass = mass

    def integ(self, qmes):
        self.q += (self.zdes - qmes)

    def thrust(self, zmes, zmesdot):
        u = self.kp*(zdes - zmes) + self.ki*self.q - self.kd*zmesdot + self.mass*9.81
        return u

class uroll:
    def __init__(self,kp, ki, kd, rolldes):
        self.kp = kp
        self.kd = kd
        self.ki = ki
        self.rolldes = rolldes

    def integ(self, rollmes):
        self.q += (self.rolldes - rollmes)

    def roll(self, rollmes, rollmesdot):
        r = self.kp*(rolldes - rollmes) + self.ki*self.q - self.kd*rollmesdot
        return r

class upitch:
    def __init__(self,kp, ki, kd, pitchdes):
        self.kp = kp
        self.kd = kd
        self.ki = ki
        self.pitchdes = pitchdes

    def integ(self, pitchmes):
        self.q += (self.pitchdes - pitchmes)

    def pitch(self, pitchmes, pichmesdot):
        p = self.kp*(pitchdes - pitchmes) + self.ki*self.q - self.kd*pitchmesdot
        return p

class uyaw:
    def __init__(self,kp, ki, kd, yawdes):
        self.kp = kp
        self.kd = kd
        self.ki = ki
        self.yawdes = yawdes

    def integ(self, yawmes):
        self.q += (self.yawdes - qmes)

    def yaw(self, yawmes, yawmesdot):
        y = self.kp*(yawdes - yawmes) + self.ki*self.q - self.kd*yawmesdot
        return y

phides = ((sin(si)*xm.roll()) - (cos(si)*ym.pitch()))/9.81
thetades = ((cos(si)*xm.roll()) + (sin(si)*ym.pitch()))/9.81


def motin(kf,km,d):
    b = np.array([[kf, kf, kf, kf], [kf*d, 0, -kf*d, 0], [0, kf*d, 0, -kf*d], [km, km, km, km]])
    w = np.dot(inv(np.matrix(b)),b)









   
