import math
import numpy as np
import matplotlib.pyplot as plot

#pass one set of values at a time for kv
#pass to this function scalars
#x->the front position at this time step
#y->previous x values not including the current time step
#t->the current time step
#tau->previous time values not including the current time step

def kv(x,y,t,tau):
  g=1/math.sqrt(4*math.pi)*math.exp(-(x-y)**2/(4*(t-tau)))
  return g

#pass one set of values at a time for kd
#pass to this function scalars
#x->the front position at this time step
#y->previous x values not including the current time step
#t->the current time step
#tau->previous time values not including the current time step

def kd(x,y,t,tau):
  g=math.exp(-(x-y)**2/(4*(t-tau)))*(x-y)/(4*math.sqrt(math.pi)*(t-tau))
  return g

#this is our given value of distribution of temperature through the medium
#pass to this function scalars
#x->location
#y->xo
#t->the current time step
#to->initial time value to!=0

def u(x,y,t,tau):
  uRet=1/math.sqrt(t+tau)*math.exp(-(x-y)**2/(4*(t+tau)))
  return uRet

#this is our given value of initial distribution of temperature through the medium
#pass to this function scalars
#x->the front position at this time step
#to->initial time value to!=0

def uo(x,xo,to):
  uRet=u(x,xo,0,to)
  return uRet

#this is the boundary temperature
#pass this function a scalar
#predetermined for the test

def rKnown(t):
  return 1
  
def dRKnown(t):
  return 0

#how to determine the single layer potential V
#pass this function scalars
#t->the time list not including the current time step
#tn->the current time step
#h->the time jump
#passer-> a list of old r values to be indexed to reduce calculation
#I could maybe insert to the passser the new q value before we throw it into this loop and therefore don't need the x from the rpasser

def sLay(t,tn,h,rpasser,qpasser):
  sum1=0
  for i in range(len(t)):
    if len(t)==1:
      sum1=sum1+kv(rKnown(tn),rpasser[i],tn,t[i])*qpasser[i]/(math.sqrt(tn-t[i]))
    else:
      if i==0:
        sum1=sum1+1/2*kv(rKnown(tn),rpasser[i],tn,t[i])*qpasser[i]/(math.sqrt(tn-t[i]))
      else:
        sum1=sum1+kv(rKnown(tn),rpasser[i],tn,t[i])*qpasser[i]/(math.sqrt(tn-t[i]))
  return sum1*(h)

#how to determine the single layer potential V
#pass this function scalars
#t->the time list not including the current time value
#tn->the current time value
#h->the time step

def munSLay(t,tn,h):
  sum2=0
  for i in range(len(t)):
    if len(t)==1:
      sum2=sum2+h/math.sqrt(tn-t[i])
    else:
      if i==0:
        sum2=sum2+1/2*h/math.sqrt(tn-t[i])
      else:
        sum2=sum2+h/math.sqrt(tn-t[i])
  return 2*math.sqrt(tn)-sum2

#how to determine the double layer potential V
#pass this function scalars
#t->the time list not including the current time step
#tn->the current time step
#h->the time jump
#passer-> a list of old r values to be indexed to reduce calculation

def dLay(t,tn,h,rpasser,upasser):
  sum3=0
  for i in range(len(t)):
    if len(t)==1:
      sum3=sum3+kd(rKnown(tn),rpasser[i],tn,t[i])*upasser[i]/(math.sqrt(tn-t[i]))
    else:
      if i==0:
        sum3=sum3+1/2*kd(rKnown(tn),rpasser[i],tn,t[i])*upasser[i]/(math.sqrt(tn-t[i]))
      else:
        sum3=sum3+kd(rKnown(tn),rpasser[i],tn,t[i])*upasser[i]/(math.sqrt(tn-t[i]))
  return sum3*h

#how to determine the initial potential
#pass this function scalars
#x->newest front point
#t->current time
  
def initPot(t,n,xo,to):
  b=5
  a=min(b,rKnown(t)/math.sqrt(4*t))
  zlist=list(np.linspace(-a,b,n))
  h=zlist[1]-zlist[0]
  sum4=0
  for i in range(len(zlist)):
    if i==0:
      sum4=sum4+1/2*h*math.exp(-zlist[i]**2)*uo(rKnown(t)+zlist[i]*math.sqrt(4*t),xo,to)*1/math.sqrt(math.pi)
    elif i==n-1:
      sum4=sum4+1/2*h*math.exp(-zlist[i]**2)*uo(rKnown(t)+zlist[i]*math.sqrt(4*t),xo,to)*1/math.sqrt(math.pi)
    else:
      sum4=sum4+h*math.exp(-zlist[i]**2)*uo(rKnown(t)+zlist[i]*math.sqrt(4*t),xo,to)*1/math.sqrt(math.pi)
  return sum4

#from the mathematica file names JL Research 1
# given rKnown=t  
def qKnowna(t,to,xo):
  ans=math.exp(-(rKnown(t)-xo)**2/(4*(t+to)))*(2*t-rKnown(t)+xo+2*to)/(4*math.sqrt(math.pi)*(t+to)**(3/2))
  return ans

#given rKnown=1

def qKnownb(t,to,xo):
  ans=math.exp(-(rKnown(t)-xo)**2/(4*(t+to)))*(-rKnown(t)+xo)/(2*(t+to)**(3/2))
  return ans

def algTester(tmax,n,qfun,to,xo):
  t=list(np.linspace(0,tmax,n))
  q=[qfun(0,to,xo)]
  qKnownList=[qfun(0,to,xo)]
  r=[rKnown(0)]
  uPasser=[u(rKnown(0),xo,0,to)]
  diff=[0]
  sLays=[0]
  dLays=[0]
  munSLays=[0]
  initPotent=[0]
  for i in range(1,len(t)):
    q.append((-1/2*u(rKnown(t[i]),xo,t[i],to)-sLay(t[0:i],t[i],t[1]-t[0],r,q)+dLay(t[0:i],t[i],t[1]-t[0],r,uPasser)+dRKnown(t[i])/(4*math.sqrt(math.pi))*munSLay(t[0:i],t[i],t[1]-t[0])*u(rKnown(t[i]),xo,t[i],to)+initPot(t[i],1000,xo,to))/(munSLay(t[0:i],t[i],t[1]-t[0])*1/math.sqrt(4*math.pi)))
    r.append(rKnown(t[i]))
    uPasser.append(u(rKnown(t[i]),xo,t[i],to))
    qKnownList.append(qfun(t[i],to,xo))
    sLays.append(sLay(t[0:i],t[i],t[1]-t[0],r,q))
    munSLays.append(munSLay(t[0:i],t[i],t[1]-t[0]))
    dLays.append(dLay(t[0:i],t[i],t[1]-t[0],r,uPasser))
    initPotent.append(initPot(t[i],1000,xo,to))
  
  for i in range(len(uPasser)):
    uPasser[i]=uPasser[i]*1/2
    
  return [q,qKnownList,r,uPasser,t,sLays,munSLays,initPotent,dLays]

qAns=algTester(5,6,qKnownb,2,2)

print('t')
print(qAns[4])

print('uPasser')
print(qAns[3])

print('init')
print(qAns[7])

print('slay')
print(qAns[5])

print("qknown")
print(qAns[1])

plot.figure(1)
mine,=plot.plot(qAns[0],label="mine")
exact,=plot.plot(qAns[1],label="exact")
uPass,=plot.plot(qAns[3],label="-uPass")
slayers,=plot.plot(qAns[5],label="-single")
initial,=plot.plot(qAns[7],label="+Init")
difference,=plot.plot(qAns[8],label="DLays+")
plot.legend(handles=[mine,exact,uPass,slayers,initial,difference])
plot.title("q")
plot.gcf()
plot.savefig('Figure 1')
plot.clf()