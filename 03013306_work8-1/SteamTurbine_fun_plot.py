# -*- coding: utf-8 -*- 
"""%matplotlib inline"""
import seuif97 as if97
import matplotlib.pyplot as plt
import numpy as np
import string as str
import SteamTurbine_fun_eff_calc as calc

def linepoint3(sstep,p,t,h,s,sis=None):
    if sis==None:
        s0 = s - sstep
    else:
        s0=sis-sstep  
    s2 = s + sstep
    point_h = np.zeros(shape=3)
    point_h[0] = if97.ps2h(p, s0)
    point_h[1] = h
    point_h[2] = if97.ps2h(p, s2)
    point_s = np.zeros(shape=3)
    point_s[0] = s0
    point_s[1] = s
    point_s[2] = s2
    return point_s,point_h

def linepoint2(h0,h1,s0,s1):
    point_h = np.zeros(shape=2)
    point_h[0] = h0
    point_h[1] = h1
    point_s = np.zeros(shape=2)
    point_s[0] = s0
    point_s[1] = s1
    return point_s,point_h

# p10��ѹ
point_p10_s,point_p10_h=linepoint3(0.02,calc.p10,calc.t10,calc.h10,calc.s10)
# p11��ѹ
point_p11_s,point_p11_h=linepoint3(0.02,calc.p11,calc.t11,calc.h11,calc.s11)

# p2��ѹ��
point_p2_s,point_p2_h=linepoint3(0.02,calc.p2,calc.t2,calc.h2,calc.s2,calc.s11)

# p11-p2���ؽ���
point_is_s11_2,point_is_h11_2=linepoint2(calc.h11,calc.h3,calc.s11,calc.s11)
# p11-p2 Expansion Line
point_hp_s11_2,point_hp_h11_2=linepoint2(calc.h11,calc.h2,calc.s11,calc.s2)
# p10-p2���ؽ���
point_is_s10_2,point_is_h10_2=linepoint2(calc.h10,calc.h4,calc.s10,calc.s10)
# p10-p2 Expansion Line
point_hp_s10_2,point_hp_h10_2=linepoint2(calc.h10,calc.h2,calc.s10,calc.s2)

plt.plot(point_p10_s,point_p10_h,'bs-')

plt.plot(point_p11_s,point_p11_h,'ys-')

plt.plot(point_p2_s,point_p2_h,'bs-')

plt.plot(point_is_s11_2,point_is_h11_2,'gs--')
plt.plot(point_hp_s11_2,point_hp_h11_2,'rs-',label=" ")

plt.plot(point_is_s10_2,point_is_h10_2,'gs--')
plt.plot(point_hp_s10_2,point_hp_h10_2,'rs-')
plt.minorticks_on()


_title =( 'The internal efficiency = ' + 
r'$\frac{h11-h2}{h11-h3}$' + '=' + '{:.2f}'.format(ieff) + '%'
+'\n'+
'The external efficiency = ' + 
r'$\frac{h10-h2}{h10-h4}$' + '=' + '{:.2f}'.format(eeff) + '%')
plt.legend(loc="best",  bbox_to_anchor=[0.5, 0.5],
        ncol=2, shadow=True,title=_title)
#annotate some interesting points using the annotate command
plt.annotate('H0',
         xy=(s10, h10), xycoords='data',
         xytext=(+5, +20), textcoords='offset points', fontsize=12,
         arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

    plt.annotate('H1',
         xy=(s11, h11), xycoords='data',
         xytext=(+10, +20), textcoords='offset points', fontsize=12,
         arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

    plt.annotate('H2',
         xy=(s2, h2), xycoords='data',
         xytext=(+10, +30), textcoords='offset points', fontsize=12,
         arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

    plt.annotate('H3',
         xy=(s11, h3), xycoords='data',
         xytext=(+10, +20), textcoords='offset points', fontsize=12,
         arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

    plt.annotate('H4',
         xy=(s10, h4), xycoords='data',
         xytext=(-10, +20), textcoords='offset points', fontsize=12,
         arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

plt.xlabel('s')
plt.ylabel('h')
plt.show()
