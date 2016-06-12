"""#---------------------------------------------------------------------------
# IAPWS-IF97-Rev: August 2007：
#       Basic Equation for Region 1
# 
#        (p,T) -> volreg1, energyreg1,entropyreg1,enthalpyreg1
#                 cpreg1, cvreg1,spsoundreg1
# ---------------------------------------------------------------------------
"""
# -*- coding: utf-8 -*-
import math

rgas_water=0.461526     # gas constant  KJ/(kg K)
tc_water=647.096        # critical temperature  K
pc_water=22.064         # critical pressure Mpa
dc_water=322.0          # critical density kg/m**3

# Table 2. Page2 ,34 I,J,N
I=0
J=1
N=2

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

# Table 4 ，Page8 
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

# Table 3 Page8
#  p,T -> *

def volreg1(p,T):            #v(p,T)
    """  specific volume in region 1""" 
    tau=1386.0/T
    pi=p/16.53 
    return 0.001*rgas_water*T*pi*gammapireg1(pi,tau)/p

def energyreg1(p,T):         #u(p,T)
    """  specific internal energy in region 1""" 
    tau=1386.0/T
    pi=p/16.53
    return rgas_water*T*(tau*gammataureg1(pi,tau)-pi*gammapireg1(pi,tau))

def entropyreg1(p,T):        #s(p,T)
    """ specific entropy in region 1""" 
    tau=1386.0/T
    pi=p/16.53
    return rgas_water*(tau*gammataureg1(pi,tau)-gammareg1(pi,tau))

def enthalpyreg1(p,T):       #h(p,T)
    """ specific enthalpy in region 1""" 
    tau=1386.0/T
    pi=p/16.53
    return rgas_water*T*tau*gammataureg1(pi,tau)

def cpreg1(p,T):             #Cp(p,T)
    """ specific isobaric heat capacity in region 1""" 
    tau=1386.0/T
    pi=p/16.53
    return -rgas_water*tau*tau*gammatautaureg1(pi,tau)

def cvreg1(p,T):             #Cv(p,T)
    """  specific isochoric heat capacity in region 1
        cvreg1 in kJ/(kg K),T in K, p in MPa""" 
    tau=1386.0/T
    pi=p/16.53
    a=-tau*tau*gammatautaureg1(pi,tau)
    b=gammapireg1(tau,pi)-tau*gammapitaureg1(pi,tau)
    b *= b;
    return rgas_water*(a+b/gammapipireg1(pi,tau))

def spsoundreg1(p,T):        #w(p,T)
    """  speed of sound in region 1
        spsoundreg1 in m/s, T in K, p in Mpa"""  
    tau=1386.0/T
    pi=p/16.53
    gammapi=gammapireg1(pi,tau)
    a=gammapi-tau*gammapitaureg1(pi,tau)
    a *= a
    b=a/(tau*tau*gammatautaureg1(pi,tau))
    b=b-gammapipireg1(pi,tau)
    return rgas_water*(a+b/gammapipireg1(pi,tau))

def thetaPHreg1(pi,  eta):
    """  Dimensionless form Backward equation T(p,h) for region 1
            Initialize coefficients and exponents
              IAPWS-IF97-Rev :Page 11, Table6 :
    """
    IJN =[
      [0, 0, -0.23872489924521E+03]
     ,[0, 1, 0.40421188637945E+03]
     ,[0, 2, 0.11349746881718E+03]
     ,[0, 6, -0.58457616048039E+01]
     ,[0, 22, -1.528548241314E-04]

     ,[0, 32, -0.10866707695377E-05]
     ,[1, 0, -0.13391744872602E+02]
     ,[1, 1,  0.43211039183559E+02]
     ,[1, 2,  -0.54010067170506E+02]
     ,[1, 3,  0.30535892203916E+02]

     ,[1, 4, -0.65964749423638E+01]
     ,[1,10,  0.93965400878363E-02]
     ,[1,32,  0.11573647505340E-06]
     ,[2,10, -0.25858641282073E-04]
     ,[2,32, -0.40644363084799E-08]
          
     ,[3,10, 0.66456186191635E-07]
     ,[3,32, 0.80670734103027E-10]
     ,[4,32,-0.93477771213947E-12]
     ,[5,32,  0.58265442020601E-14]
     ,[6,32,  -0.15020185953503E-16]
        ]  
   
    eta=eta+1.0
    theta=0
    for k in range(20):
       theta += IJN[k][2]*pow(pi,IJN[k][0])*pow(eta,IJN[k][1])
    return theta


def phtoTreg1(p,h):         #T(p,h)
    pi=p/1.0
    eta=h/2500.0
    return (1.0*thetaPHreg1(pi,eta))

def thetaPStoT(pi,sigma):
    """  Dimensionless form Backward equation T(p,s) for region 1
            Initialize coefficients and exponents
               IAPWS-IF97-Rev :Page 12, Table 8 :
    """
    IJN=[ [ 0, 0, 0.17478268058307E+03]
         ,[ 0, 1, 0.34806930892873E+02]
         ,[ 0, 2, 0.65292584978455E+01]
         ,[ 0, 3, 0.33039981775489]
         ,[ 0,11, -0.19281382923196E-06]

         ,[ 0, 31, -0.24909197244573E-22]
         ,[ 1, 0,  -0.26107636489332]
         ,[ 1, 1,   0.22592965981526]
         ,[ 1, 2,  -0.64256463395226E-01]
         ,[ 1, 3,   0.78876289270526E-02]

         ,[ 1, 12,  0.35672110607366E-9]
         ,[ 1, 31,  0.17332496994895E-23]
         ,[ 2,  0,  0.56608900654837E-03]
         ,[ 2,  1,  -0.32635483139717E-03]
         ,[ 2,  2,   0.44778286690632E-04]

         ,[ 2, 9,   -0.51322156908507E-9]
         ,[ 2, 31,  -0.42522657042207E-25]
         ,[ 3, 10,   0.26400441360689E-12]
         ,[ 3, 32,   0.78124600459723E-28]
         ,[ 4, 32,   -0.30732199903668E-30]
        ]

    sigma=sigma+2.0
    theta=0.0
    for k in range(20):
      theta += IJN[k][2]*pow(pi,IJN[k][0])*pow(sigma,IJN[k][1])
    return theta

def pstoTreg1(p,s):         #T(p,s)
    pi=p/1.0
    sigma=s/1.0
    return (1.0*thetaPStoT(pi,sigma))


def piHS(eta, sigma):
    """  Dimensionless form Backward equation P(h,s) for region 1
          Initialize coefficients and exponents 
            Supp-PHS12-2014:  Table 2 Page 5:
    """
    IJN=[
          [ 0, 0,  -0.691997014660582]
         ,[ 0, 1,  -0.183612548787560E+02]
         ,[ 0, 2,  -0.928332409297335E+01]
         ,[ 0, 4,   0.659639569909906E+02]
         ,[ 0, 5,   -0.162060388912024E+02]

         ,[ 0, 6,  0.450620017338667E+03]
         ,[ 0, 8,  0.854680678224170E+03]
         ,[ 0, 14, 0.607523214001162E+04]
         ,[ 1, 0,  0.326487682621856E+02]
         ,[ 1, 1, -0.269408844582931E+02]

         ,[ 1, 4, -0.319947848334300E+03]
         ,[ 1, 6, -0.928354307043320E+03]
         ,[ 2, 0,  0.303634537455249E+02]
         ,[ 2, 1, -0.650540422444146E+02]
         ,[ 2, 10,-0.430991316516130E+04]

         ,[ 3, 4, -0.747512324096068E+03]
         ,[ 4, 1,  0.730000345529245E+03]
         ,[ 4, 4,  0.114284032569021E+04]
         ,[ 5, 0, -0.436407041874559E+03]
        ] 

    eta=eta+0.05
    sigma=sigma+0.05
    pi=0.0
    for k in range(19):
      pi += IJN[k][2]*pow(eta,IJN[k][0])*pow(sigma,IJN[k][1])
    return pi

def hstopreg1(h,s):    #p(h,s)
    eta=h/3400.0
    sigma=s/7.6
    return (100.0*piHS(eta,sigma))


