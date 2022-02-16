from numpy import zeros

def linear_shooting(a,b,alpha,beta,N):
    
    h   = (b-a)/N
    
    u1 = zeros(N+1)        # try initialize 
    u2 = zeros(N+1)
    v1 = zeros(N+1)
    v2 = zeros(N+1)
    w1 = zeros(N)
    w2 = zeros(N)
    
    u1[0] = alpha
    u2[0] = 0
    v1[0] = 0
    v2[0] = 1
    print("x", "      " , " w1", "      " , "w2")
    for i in range(0,N):
        
        x = a + i*h
        
        k11 = h*u2[i]
        k12 = h*(p(x)*u2[i] + q(x)*u1[i] + r(x))
        k21 = h*(u2[i] + (1/2)*k12)
        k22 = h*(p(x + h/2)*(u2[i] + (1/2)*k12 + q(x + h)*(u1[i] + 0.5*k11) + r(x + h/2)))
        k31 = h*(u2[i] + (1/2)*k22)
        k32 = h*(p(x + h/2)*(u2[i] + (1/2)*k22) + q(x + h/2)*(u1[i] + (1/2)*k21) + r(x+h/2))
        k41 = h*(u2[i] + k32)
        k42 = h*(p(x + h)*(u2[i] + k32) + q(x + h)*(u1[i] + k31) + r(x + h))
        
        u1[i+1] = u1[i] + (1/6)*(k11 + 2*k21 + 2*k31 + k41)
        u2[i+1] = u2[i] + (1/6)*(k12 + 2*k22 + 2*k32 + k42)
        
        kp11 = h*v2[i]
        kp12 = h*(p(x)*v2[i] + q(x)*v1[i])
        kp21 = h*(v2[i] + (1/2)*kp12)
        kp22 = h*(p(x + h/2)*(v2[i] + (1/2)*kp12) + q(x + h/2)*(v1[i]+ (1/2)*kp11))
        kp31 = h*(v2[i] + (1/2)*kp22)
        kp32 = h*(p(x + h/2)*(v2[i] + (1/2)*kp22) + q(x + h/2)*(v1[i] + (1/2)*kp12))
        kp41 = h*(v2[i] + kp32)
        kp42 = h*(p(x + h)*(v2[i] + kp32) + q(x+ h)*(v1[i] + kp31))
        
        v1[i+1] = v1[i] + (1/6)*(kp11 + 2*kp21 + 2*k31 + kp41)
        v2[i+1] = v2[i] + (1/6)*(kp12 + 2*kp22 + 2*k32 + kp42)
        
    w1[0] = alpha
    w2[0] = (beta - u1[N])/(v1[N])
    print("%.2f"  "    "  "%.5f"    "    " "%.5f" %(a,w1[0],w2[0]))
    
    for i in  range(1, N+1):
        W1 = u1[i] + w2[0]*v1[i]
        W2 = u2[i] + w2[0]*v2[i]
        x = a + i*h
        print("%.2f"  "    "  "%.5f"  "    "  "%.5f"    %(x,W1,W2))
    return   

    #########################################################################
from numpy import sin, log

def p(x):
    px = -2/x
    return px

def q(x):
    qx = 2/x**2
    return qx

def r(x):
    rx = (sin(log(x)))/x**2
    return rx

a  = 1
b  = 2
alpha = 1
beta = 2
N = 10


#################################################
linear_shooting(a,b,alpha,beta,N)
