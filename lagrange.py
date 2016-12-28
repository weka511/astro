import math as m

def fn(r,mu):
    return 3*r*r*r*(1-r+r*r/3)/((1+r+r*r)*(1-r)*(1-r)*(1-r))-mu/(1-mu)

def solve(r0,r1,mu,f0=None,f1=None,depth=500,eps=1e-12):
    r =0.5*(r0+r1)
    f = fn(r,mu)
    print (depth,r0,r1,r1-r0,r,f)
    
    if depth==0 or r1-r0<eps:
        return r
    if f0==None:
        f0=fn(r0,mu)
        f1=fn(r1,mu)

    s0=m.copysign(1.0, f0)
    s1=m.copysign(1.0, f1)    
    s=m.copysign(1.0, f)
    if s0==s:
        return solve(r,r1,mu,f,f1,depth=depth-1)
    if s==s1:
        return solve(r0,r,mu,f0,f,depth=depth-1)
    
def get_alpha(mu):
    return (mu/(3*(1-mu)))**(1/3)

def get_u(r,mu):
    mu1=1-mu
    r1=1-r
    return mu1*(1/r1+r1*r1/2) + mu*(1/r+r*r/2)-mu1*mu/2

alpha=get_alpha(0.2)

r=solve(alpha-0.1,alpha+0.05,0.2)
u=get_u(r,0.2)
print (2*u)
