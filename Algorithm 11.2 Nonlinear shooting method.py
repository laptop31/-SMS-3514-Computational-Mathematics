from numpy import zeros, abs

def shoot_nonlinear(a,b,alpha, beta, n, tol, M):
    
    w1 = zeros(n+1)  
    w2 = zeros(n+1)
    h = (b-a)/n
    k = 1
    TK = (beta - alpha)/(b - a)
    
    print("i""  x" "       " "W1""        " "W2")
    while k <= M:
        
        w1[0] = alpha
        w2[0] = TK
        u1    = 0
        u2    = 1
        
        #print("w1 = ", w1[0])
        #print("w2 = ", w2[0])
        #print("u1 = ", u1)
        #print("u2 = ", u2)
        #print("w1 = ", w1)
        #print("w2 = ", w2)
        
       
        #print("i" " " " x")
        #print("i" " " " k1")
        for i in range(1,n+1):
            x = a + (i-1)*h     #step 5
            
            t = x + 0.5*(h)
            
            k11 = h*w2[i-1]     #step 6
    
            k12 = h*f(x,w1[i-1],w2[i-1])
            k21 = h*(w2[i-1] + (0.5)*k12)
            k22 = h*f(t, w1[i-1] + (0.5)*k11, w2[i-1] + (0.5)*k12)
            k31 = h*(w2[i-1] + (0.5)*k22)
            k32 = h*f(t, w1[i-1] + (0.5)*k21, w2[i-1] + (0.5)*k22)
            t   = x + h
            k41 = h*(w2[i-1]+k32)
            k42 = h*f(t, w1[i-1] + k31, w2[i-1] + k32)
            w1[i] = w1[i-1] + (1/6)*(k11 + 2*k21 + 2*k31 + k41)
            w2[i] = w2[i-1] + (1/6)*(k12 + 2*k22 + 2*k32 + k42)  
            kp11 = h*u2
            kp12 = h*(fy(x,w1[i-1],w2[i-1])*u1 + fyp(x,w1[i-1], w2[i-1])*u2)
            t    = x + 0.5*(h)
            kp21 = h*(u2 + (0.5)*kp12)
            kp22 = h*((fy(t, w1[i-1],w2[i-1])*(u1 + (1/2)*kp11)) + fyp(x+h/2, w1[i-1],w2[i-1])*(u2 +(0.5)*kp12))
            kp31 = h*(u2 + (0.5)*kp22)
            kp32 = h*((fy(t, w1[i-1],w2[i-1])*(u1 + (1/2)*kp21)) + fyp(x+h/2, w1[i-1],w2[i-1])*(u2 +(0.5)*kp22))
            t    = x + h
            kp41 = h*(u2 + kp32)
            kp42 = h*(fy(t, w1[i-1], w2[i-1])*(u1+kp31) + fyp(x + h, w1[i-1], w2[i-1])*(u2 + kp32))
            u1 = u1 + (1/6)*(kp11 + 2*kp21 + 2*kp31 + kp41)
            u2 = u2 + (1/6)*(kp12 + 2*kp22 + 2*kp32 + kp42)
           
        
        r = abs(w1[n]) - beta
        #print(r)
        if r < tol:
            for i in range(0,n+1):
                x = a + i*h
                print("%.2f %.2f %.4f %.4f" %(i,x,w1[i],w2[i]))
            return
        #else:
        TK = TK -(w1[n]-beta)/u1
        
        k = k+1
    
            
    print("Maximum number of iterations exceeded")   
    return


#############################################################################################

def f(x,y,yp):
    fx = (0.125)*(32 + 2*x**3 -y*yp)
    return fx

def fy(xp,z,zp):
    fyy = -(0.125)*(zp)
    return fyy

def fyp(xpp,zpp,zppp):
    fypp = -(0.125)*(zpp)
    return fypp

 ##########################################################################

 a = 1         # start point
b = 3         # end point
alpha = 17    # boundary condition
beta = 14.33333   # boundary condition
N = 20        # number of subintervals
M = 10        # maximum number of iterations
tol = 0.000001 # tolerance

###############################################################################
shoot_nonlinear(a,b,alpha,beta,N,tol,M)   