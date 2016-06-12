# -*- coding: utf-8 -*-

import math

from prettytable import PrettyTable

rgas_water=0.461526     # gas constant  KJ/(kg K)
tc_water=647.096        # critical temperature  K
pc_water=22.064         # critical pressure Mpa
dc_water=322.0          # critical density kg/m**3
                 
parameter={'p':[],'T':[],'v':[],'u':[],'s':[],'h':[],'Cp':[],'Cv':[],'w':[]}

# Table 2. Page2 ,34 I,J,N
I=0
J=1
N=2
number=0

IJN =[   [0,    -2,    0.14632971213167E+00]
    ,[0,    -1,    -0.84548187169114E+00]
    ,[0,    0,    -0.37563603672040E+01]
    ,[0,    1,    0.33855169168385E+01]
    ,[0,    2,    -0.95791963387872E+00]
    ,[0,    3,    0.15772038513228E+00]
    ,[0,    4,    -0.16616417199501E-01]
    ,[0,    5,    0.81214629983568E-03]
    ,[1,    -9,    0.28319080123804E-03]
    ,[1,    -7,    -0.60706301565874E-03]
    ,[1,    -1,    -0.18990068218419E-01]
    ,[1,    0,    -0.32529748770505E-01]
    ,[1,    1,    -0.21841717175414E-01]
    ,[1,    3,    -0.52838357969930E-04]
    ,[2,    -3,    -0.47184321073267E-03]
    ,[2,    0,    -0.30001780793026E-03]
    ,[2,    1,    0.47661393906987E-04]
    ,[2,    3,    -0.44141845330846E-05]
    ,[2,    17,    -0.72694996297594E-15]
    ,[3,    -4,    -0.31679644845054E-04]
    ,[3,    0,    -0.28270797985312E-05]
    ,[3,    6,    -0.85205128120103E-09]
    ,[4,    -5,    -0.22425281908000E-05]
    ,[4,    -2,    -0.65171222895601E-06]
    ,[4,    10,    -0.14341729937924E-12]
    ,[5,    -8,    -0.40516996860117E-06]
    ,[8,    -11,    -0.12734301741641E-08]
    ,[8,    -6,    -0.17424871230634E-09]
    ,[21,    -29,    -0.68762131295531E-18]
    ,[23,    -31,    0.14478307828521E-19]
    ,[29,    -38,    0.26335781662795E-22]
    ,[30,    -39,    -0.11947622640071E-22]
    ,[31,    -40,    0.18228094581404E-23]
    ,[32,    -41,    -0.93537087292458E-25]
     ]

# Table 4 ��Page8 
def gammareg1(pi,tau):
    """ Fundamental equation for region 1 """
    tau=tau-1.222
    pi=7.1-pi
    reg1=0
    for k in range(34):
      reg1 += IJN[k][N]*math.pow(pi,IJN[k][I])*math.pow(tau,IJN[k][J])
    return reg1

def gammapireg1(pi,tau):
    """ First derivative of fundamental equation in pi for region 1"""
    tau=tau-1.222
    pi=7.1-pi
    reg1=0
    for k in range(34):
        reg1 -= IJN[k][N]*IJN[k][I]*math.pow(pi,IJN[k][I]-1)*math.pow(tau,IJN[k][J])
    return reg1

def gammapipireg1(pi,tau):
    """ Second derivative of fundamental equation in pi for region 1"""
    tau=tau-1.222
    pi=7.1-pi
    reg1=0
    for k in range(34):
      reg1 +=IJN[k][N]*IJN[k][I]*(IJN[k][I]-1)*math.pow(pi,IJN[k][I]-2)*math.pow(tau,IJN[k][J])
    return reg1

def gammataureg1(pi,tau):
    """ First derivative of fundamental equation in tau for region 1""" 
    tau=tau-1.222
    pi=7.1-pi
    reg1=0
    for k in range(34):
       reg1 += IJN[k][N]*math.pow(pi,IJN[k][I])*IJN[k][J]*math.pow(tau,IJN[k][J]-1)
    return reg1

def gammatautaureg1(pi,tau):
    """ Second derivative of fundamental equation in tau for region 1 """
    tau=tau-1.222
    pi=7.1-pi
    reg1=0
    for k in range(34):
       reg1 +=IJN[k][N]*math.pow(pi,IJN[k][I])*IJN[k][J]*(IJN[k][J]-1)*math.pow(tau,IJN[k][J]-2);
    return reg1

def gammapitaureg1(pi,tau):
    """  Second derivative of fundamental equation in pi and tau for region 1""" 
    tau=tau-1.222
    pi=7.1-pi
    reg1=0
    for k in range(34):
       reg1 -= IJN[k][N]*IJN[k][I]*pow(pi,IJN[k][I]-1)*IJN[k][J]*pow(tau,IJN[k][J]-1)
    return reg1

def cal_pT(p,T):
    
    tau=1386.0/T
    pi=p/16.53 
     
    v= 0.001*rgas_water*T*pi*gammapireg1(pi,tau)/p
    
    u= rgas_water*T*(tau*gammataureg1(pi,tau)-pi*gammapireg1(pi,tau))
    
    s= rgas_water*(tau*gammataureg1(pi,tau)-gammareg1(pi,tau))
    
    h= rgas_water*T*tau*gammataureg1(pi,tau)
    
    Cp= -rgas_water*tau*tau*gammatautaureg1(pi,tau)
    
    a=-tau*tau*gammatautaureg1(pi,tau)
    b=gammapireg1(tau,pi)-tau*gammapitaureg1(pi,tau)
    b *= b
    
    Cv= rgas_water*(a+b/gammapipireg1(pi,tau))
    
    gammapi=gammapireg1(pi,tau)
    c=gammapi-tau*gammapitaureg1(pi,tau)
    c *= c
    d=c/(tau*tau*gammatautaureg1(pi,tau))
    d=d-gammapipireg1(pi,tau)
    
    w= rgas_water*(c+d/gammapipireg1(pi,tau))
   
    
    parameter['p'].append(float(p))
    parameter['T'].append(float(T))
    parameter['v'].append(float(v))
    parameter['u'].append(float(u))
    parameter['s'].append(float(s))
    parameter['h'].append(float(h))
    parameter['Cp'].append(float(Cp))
    parameter['Cv'].append(float(Cv))
    parameter['w'].append(float(w))
    

def processing_table_TextFiles(i):

    table = PrettyTable(["p",
                         "T", "v", "u", "s","h",
                         "Cp", "Cv", "w"])
    table.align= "m"
    table.padding_width = 1 # One space between column edges and contents (default)
    for i in parameter:
        table.add_row(parameter['p'][i],
                   "%.3f" % parameter['T'][i],
                   "%.3f" % parameter['v'][i],"%.3f"  % parameter['u'][i],
                   "%.3f" % parameter['s'][i],"%.3f"  % parameter['h'][i],
                   "%.3f" % parameter['Cp'][i],
                   "%.3f" % parameter['Cv'][i],"%.3f" % parameter['w'][i])
    print(table)

#main 
cal_pT( 3, 300)
cal_pT(80, 300)
cal_pT( 3, 500)

processing_table_TextFiles(1)












