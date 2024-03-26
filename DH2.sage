from ecdsa import NIST256p

E = NIST256p.curve
G = NIST256p.generator
p = 115792089210356248762697446949407573530086143415290314195533631308867097853951
a = -3
b = 41058363725152142129326129780047268409114441015993725554835256314039467401291
E = EllipticCurve(GF(p),[a,b])
G = E.lift_x(48439561293906451759052585252797914202762949526041747995844080717082404635286)
n = G.order()
def leak(Point,k):
    k = 256 - k
    t1 = (ZZ(Point.xy()[0]) >> k) << k
    t2 = ZZ(Point.xy()[0])-t1
    return t1,t2
def Get(a,R):
    return leak(a*R,164)
def Coe(H,Qx,h0,a,b,p):
    A = H*(h0-Qx)^2-2*h0^2*Qx-2*(a+Qx^2)*h0-2*a*Qx-4*b
    B = 2*(H*(h0-Qx)-2*h0*Qx-a-Qx^2)
    C = (H-2*Qx)
    D = (h0-Qx)^2
    E = 2*(h0-Qx)
    return A%p,B%p,C%p,D%p,E%p
aa = randint(0,n)
bb = randint(0,n)
AA = aa*G
BB = bb*G
SS = a*b*G
R1 = BB + 0*G
h0,e0 = Get(aa,R1)
Hs = []
EE = []
Qs = []
for i in range(1,6):
    h,e = Get(aa,BB+i*G)
    h_,e_ = Get(aa,BB-i*G)
    Q = i*AA
    h = h + h_
    e = e + e_
    Hs.append(h)
    EE.append(e)
    Qs.append(Q)
As,Bs,Cs,Ds,Es = [],[],[],[],[]
for i in range(5):
    Ai,Bi,Ci,Di,Ei = Coe(Hs[i],Qs[i].xy()[0],h0,a,b,p)
    As.append(ZZ(Ai))
    Bs.append(ZZ(Bi))
    Cs.append(ZZ(Ci))
    Ds.append(ZZ(Di))
    Es.append(ZZ(Ei))
PR.<x0,y1,y2,y3,y4,y5> = PolynomialRing(ZZ)
y = [y1,y2,y3,y4,y5]
Fs = []
for i in range(5):
    fi = As[i] + Bs[i]*x0 + Cs[i]*x0^2 + Ds[i]*y[i] + Es[i]*x0*y[i] + x0^2*y[i]
    assert fi(*([e0]+EE))%p == 0
    Fs.append(fi)
Exp_res =[]
for i in range(2^5):
    t = [int(j) for j in bin(i)[2:].zfill(5)]
    Exp_res.append(t)
def Get_js(exp):
    j = [exp[i]*(i+1) for i in range(len(exp)) if exp[i]!=0]
    return j
def get_row1(k,js):
    res = 1
    for u in range(len(js)):
        if u!=k:
            res*=(x0^2+Es[js[u]-1]*x0+Ds[js[u]-1])
    return res
def get_row2(k,js):
    res = x0
    for u in range(len(js)):
        if u!=k:
            res*=(x0^2+Es[js[u]-1]*x0+Ds[js[u]-1])
    return res
def get_M(js):
    G = Sequence([], Fs[0].parent())
    for i in range(len(js)):
        G.append(get_row1(i,js))
    for i in range(len(js)):
        G.append(get_row2(i,js))
    M , mons = G.coefficient_matrix()
    M = M[:,::-1]
    return M
def Get_prodF(u,js,y):
    ttf = y[js[u-1]-1]
    for i in range(len(js)):
        if i!=(u-1):
            ttf *= Fs[js[i]-1]
    return ttf
ff = []
n,d,t = 5,4,0
for i0 in range(2*d):
    for l in range(d+1):
#         f_temp=1
#         bounds = [(2^86)]*6
        if 1 <= l <= d and 0 <= i0 <= (2*l-1):
            f_temp = p^(d+1-l)
        if 0 <= l <= d and 2*l <= i0 <= 2*d:
            f_temp = p^(d-l)
        if l == 0 and 0 <= i0 <= 2*d:
            ft = x0^i0
            ft = f_temp*ft.change_ring(ZZ)
            ff.append(ft)
        if l == 1 and 0 <= i0 <= 1:
            for exps in Exp_res:
                if sum(exps) == l:
                    ft = x0^i0*prod(map(pow,y,exps))
                    ft = f_temp*ft.change_ring(ZZ)
                    ff.append(ft)
        if 1 <= l <= d  and 2*l <= i0 <= 2*d:
            for exps in Exp_res:
                if sum(exps) == l:
                    ft = x0^(i0-2*l)*prod(map(pow,Fs,exps))
                    ft = f_temp*ft.change_ring(ZZ)
                    assert ft(*([e0]+EE))%p^3 == 0
#                     print(,l,i0)
                    ff.append(ft)
        if 2 <= l <= d  and 0 <= i0 <= (2*l-1):
            for exps in Exp_res:
                if sum(exps) == l:
                    ft = 0
                    js = Get_js(exps)
                    M = get_M(js)
                    M = M.change_ring(Zmod(p^(l-1)))
                    W = M.inverse()
                    W = W.change_ring(ZZ)
                    for u in range(1,l+1):
                        for v in range(2):
                            ft += W[i0,u+l*v-1] * x0^v * Get_prodF(u,js,y)
                    ft = f_temp*(ft.change_ring(ZZ) %p^(l-1))
                    ff.append(ft)
def Get_H(u,js):
    ttf = y[js[u-1]-1]
    for i in range(len(js)):
        if i!=(u-1):
            ttf *= Fs[js[i]-1]
    return ttf
def Get_J(u,js):
    ttf = Cs[js[u-1]-1]
    for i in range(len(js)):
        if i!=(u-1):
            ttf *= Fs[js[i]-1]
    return ttf
def Get_K(u,js):
    ttf = Bs[js[u-1]-1]-Cs[js[u-1]-1]*Es[js[u-1]-1]
    for i in range(len(js)):
        if i!=(u-1):
            ttf *= Fs[js[i]-1]
    return ttf

for i0 in range(t+1):
    l = d+1
#     bounds = [(2^86)]*6
    for exps in Exp_res:
        if sum(exps) == l:
#             print(exps)
            ft = 0
            js = Get_js(exps)
#             print(js)
            M = get_M(js)
            M = M.change_ring(Zmod(p^(l-1)))
            W = M.inverse()
            W = W.change_ring(ZZ)
            H,J,K = 0,0,0
            for u in range(1,l+1):
                for v in range(2):
                    H += W[i0,u+l*v-1] * x0^v * Get_H(u,js)
            for u in range(1,l+1):
                for v in range(2):
                    J +=W[i0,u+l*v-1] * x0^v * Get_J(u,js)
            for u in range(1,l+1):
                K += W[i0,u+l-1] * Get_K(u,js)
#             print(H(*([e0]+EE)) %p^3,J(*([e0]+EE)) %p^3,K(*([e0]+EE)) %p^3)
            ft = (H+J+K) %p^d
            assert ft(*([e0]+EE)) %p^d == 0
            ff.append(ft)
import itertools
from sage.all import *
from Crypto.Util.number import bytes_to_long, getPrime
from re import findall
from subprocess import check_output
import fgb_sage
def flatter(M):
    # compile https://github.com/keeganryan/flatter and put it in $PATH
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in M) + "]]"
    ret = check_output(["flatter"], input=z.encode())
    return matrix(M.nrows(), M.ncols(), map(int, findall(b"-?\\d+", ret)))
def LLL_mononimals(ff,bounds,debug = False):
    G = Sequence([], ff[0].parent())
    for i in range(len(ff)):
        G.append(ff[i])
    B, monomials = G.coefficient_matrix()
    monomials = vector(monomials)
    factors = [monomial(*(bounds)) for monomial in monomials]
    for i, factor in enumerate(factors):
        B.rescale_col(i, factor)
    print("dim",B.nrows(),B.ncols())
    print('s flatter')
    # B = B.dense_matrix().LLL()
    B = flatter(B)
    print('e flatter')
    B = B.change_ring(QQ)
    for i, factor in enumerate(factors):
        B.rescale_col(i, 1 / factor)
    R = B*monomials
    if debug:
        for i in range(len(R)):
            if R[i](*([e0]+EE)) == 0:
                print(i)
    Res = []
    V = fgb_sage.groebner_basis(R[8:20])
    XX = [x0,y1,y2,y3,y4,y5]
    for i in range(len(XX)):
        Res.append(XX[i] - V[i])
    return Res
LLL_mononimals(ff,[2^93]*6,debug = True)