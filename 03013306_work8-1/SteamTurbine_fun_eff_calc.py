import seuif97 as if97

pam=0.10

p10=16.18+pam
t10=532.91
p11=16.15+pam
t11=535.21

p2=3.85+pam
t2=344.39

h10=if97.pt2h(p10,t10)
s10=if97.pt2s(p10,t10)

h11=if97.pt2h(p11,t11)
s11=if97.pt2s(p11,t11)

h2=if97.pt2h(p2,t2)
s2=if97.pt2s(p2,t2)

h3=if97.ps2h(p2,s11)
hdis13=h11-h3   # 等熵焓降
hd112=h11-h2     # 过程焓降

ieff=100*hd112/hdis13

h4=if97.ps2h(p2,s10)
hdis14=h10-h4   # 等熵焓降
hd102=h10-h2     # 过程焓降

eeff=100*hd102/hdis14

print(h10,s10)
print(h11,s11)

print('internal efficiency: ',ieff)
print('external efficiency: ',eeff)
