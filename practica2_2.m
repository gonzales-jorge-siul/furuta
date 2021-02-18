Ptf = tf([-18], [1 0 -2055])
Pss = ss(Ptf)
Pss = canon(Pss,'companion')
ctrb(Pss)
rank(ctrb(Pss))
wn=2
z=1
Pd=roots([1 2*z*wn wn^2])
K=acker(Pss.A, Pss.B, Pd)
N=-(Pss.C*inv(Pss.A-Pss.B*K)*Pss.B)\eye(1,1)
Pclss=ss(Pss.A-Pss.B*K,Pss.B*N,Pss.C,Pss.D);
stepinfo(Pclss)
step(Pclss)

wn=2*wn
Pdl=roots([1 2*z*wn wn^2])
L=acker(Pss.A', Pss.c', Pdl)'