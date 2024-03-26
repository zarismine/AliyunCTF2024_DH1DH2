from itertools import product
from re import findall
from subprocess import check_output
import fgb_sage
def flatter(M):
    # compile https://github.com/keeganryan/flatter and put it in $PATH
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in M) + "]]"
    ret = check_output(["flatter"], input=z.encode())
#   我的发，findall破bug改半天
    from re import findall
    return matrix(M.nrows(), M.ncols(), map(int, findall(b"-?\\d+", ret)))
def small_roots(N, f, bounds, m = 4, t = 1, debug = False):
    f = f.change_ring(ZZ)
    G = Sequence([], f.parent())
    for k in range(m+1):
        i = 3*k
        for j in range(2*(m-k)+t+1): 
            G.append(y^j*f^k*N^(m-k))
        for i in range(3*k+1,3*k+3):
            for j in range(2*(m-k)+t-2+1): 
                G.append(x^(i-3*k)*y^j*f^k*N^(m-k))
    B, monomials = G.coefficient_matrix()
    monomials = vector(monomials)
    factors = [monomial(*bounds) for monomial in monomials]
    for i, factor in enumerate(factors):
        B.rescale_col(i, factor)
    print(B.nrows(),B.ncols())
    print('s flatter')
    if debug:
        from time import time
        t1 = time()
#     B = B.dense_matrix().LLL()
    try:
        B = flatter(B)
        if debug:
            print(time()-t1)
    except:
        B = B.dense_matrix().LLL()
    print('e flatter')
    B = B.change_ring(QQ)
    for i, factor in enumerate(factors):
        B.rescale_col(i, 1 / factor)
    R = B*monomials
    if debug:
        
        for i in range(len(R)):
            if R[i](xa,xb)==0:
                print(i)
    try:
        V = fgb_sage.groebner_basis(R[1:40])
        return ZZ(x - V[0]),ZZ(y - V[1])
    except:
        return [0,0]
from hashlib import sha256
from Crypto.Cipher import AES
from pwn import *
from tqdm import tqdm
context.log_level = 'debug'
p = remote("crypto0.aliyunctf.com",12346)

m = 2**211
E = EllipticCurve(GF(115792089210356248762697446949407573530086143415290314195533631308867097853951),[-3,41058363725152142129326129780047268409114441015993725554835256314039467401291])
G = E.lift_x(48439561293906451759052585252797914202762949526041747995844080717082404635286)
def proof():
    p.recvuntil(b'sha256(XXXX + ')
    temp=p.recvuntil(b') == ',drop=True)
    h=p.recvuntil(b'\n',drop=True)
    print(temp,h)
    for i in tqdm(string.ascii_letters+string.digits):
        for j in string.ascii_letters+string.digits:
            for m in string.ascii_letters+string.digits:
                for n in string.ascii_letters+string.digits:
                    temp1=(i+j+m+n).encode()+temp
                    # print(temp)
                    if sha256(temp1).hexdigest()==h.decode():
                        print('true')
                        p.sendline((i+j+m+n).encode())
                        return
proof()
p.recvuntil(b"my public key: ")
Alice = ZZ(int(p.recvuntil(b'\n',drop = True).decode(),16))
Alice = E.lift_x(Alice)
p.recvuntil(b"my public key: ")
Bob = ZZ(int(p.recvuntil(b'\n',drop = True).decode(),16))
Bob = E.lift_x(Bob)
sign = bytes.fromhex(p.recvuntil(b'\n',drop = True).decode())
p.recvuntil(b'Now, it is your turn:\n')
p.sendline(hex(Bob.xy()[0])[2:].encode())
p.sendline(hex(Bob.xy()[1])[2:].encode())
p.recvuntil(b'Leak: ')
x_l = int(p.recvuntil(b', ',drop = True).decode(),16)
y_l = int(p.recvuntil(b'\n',drop = True).decode(),16)
print(x_l,y_l)
pp = 115792089210356248762697446949407573530086143415290314195533631308867097853951
bb = 41058363725152142129326129780047268409114441015993725554835256314039467401291
PR.<x,y> = PolynomialRing(ZZ)
from time import time
t1 = time()
# for i in range(0,4):
#     for j in range(0,4):
f = (x+x_l)^3 - 3*(x+x_l) + bb - (y + y_l)^2
f = f(x*2^211,y*2^211)
f = f*ZZ(inverse_mod(ZZ(f.coefficients()[0]),pp)) %pp
xx,yy = small_roots(pp,f,[2^45,2^45],7,2)
print(time()-t1)
x = x_l+(xx<<211)
y = y_l+(yy<<211)
key = x, y

print(str(key))
key = sha256(str(key).encode()).digest()
tar = AES.new(key, AES.MODE_ECB).decrypt(sign)[15:15+64]
print(tar)
p.sendline(AES.new(key, AES.MODE_ECB).encrypt(tar).hex().encode())
p.recvline()
