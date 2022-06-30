# EE2703 Applied programming Lab-2022
#End_semster
#Name:Gugulothu Naveen
#Roll:EE20B039
from cmath import pi
from matplotlib.pylab import *
#The parameters of this problem as follows:

c=2.9979e8 #speed of light
l=0.5 #quater wavelength
mu0=4e-7*pi #permeability of free space
Im=1.0 #current injected into the antenna
a=0.01 #radius of wire 
lamda=l*4.0
f=c/lamda #frequency
k=2*pi/lamda #wave number
N=4
dz=l/N #spacing of current spamples

#Q1 Characterizing the areas of all known and obscure currents
z = linspace(-l,l,2*N+1)
z_1=z[1:-1]
# Characterizing an cluster containing the areas of all obscure currents
u = delete(z,[0,N,2*N])
print((z).round(2))
print((u).round(2))
#Q2 Defining a function which computes and returns a matrix to be used in understanding the Ampere's Law
def Matrix(N,a):
    M = (1/(2*pi*a))*identity(2*N -2)
    return M
    
print((Matrix(N,a)).round(2))

#Q3
U = meshgrid(u, u)
ui = U[0]
uj = U[1]
Ru = sqrt((ui-uj)**2 + ones([2*N-2,2*N-2],dtype=complex)*(a**2))
print((Ru).round(2))



Z = meshgrid(z,z)
zi = Z[0]
zj = Z[1]

Rz = sqrt((zi-zj)**2 + ones([2*N+1,2*N+1],dtype=complex)*(a**2))

print((Rz).round(2))
RN = delete(Rz[N],[0,N,2*N])
P = (mu0 / (4 * pi)) * (exp(-k * Ru * 1j) / Ru) * dz
P_B = (mu0 / (4 * pi)) * (exp(-k * RN * 1j) / RN) * dz
Q = -(a / mu0) * P * ((-k * 1j / Ru) + (-1 / Ru ** 2))
Q_B = -(a / mu0) * P_B * ((-k * 1j / RN) + (-1 / RN ** 2))
print((RN).round(2))
print((P*1e8).round(2))
print((P_B*1e8).round(2))
print((Q).round(2))
print((Q_B).round(2))


J = dot(inv(Matrix(N, a) - Q) , Q_B) * Im
I1 = concatenate(([0],J[:N-1],[Im],J[N-1:],[0]))
I2 = Im*sin(k*(l-abs(z)))
print((J).round(2))
print((I1).round(2))
print((I2).round(2))
plot(z,abs(I1),label = "calculated Current")
plot(z,I2,label = "Assumed Current")
xlabel(r"$z$")
ylabel(r"$I$")
title("antenna currents in a half-wave dipole antenna")
grid()
legend()
savefig('f1.png')
show()