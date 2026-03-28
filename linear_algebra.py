#!/usr/bin/env python3
"""Linear algebra: determinant, inverse, eigenvalues, SVD (2x2)."""
import math
def det(A):
    n=len(A)
    if n==1: return A[0][0]
    if n==2: return A[0][0]*A[1][1]-A[0][1]*A[1][0]
    d=0
    for j in range(n):
        minor=[[A[i][k] for k in range(n) if k!=j] for i in range(1,n)]
        d+=(-1)**j*A[0][j]*det(minor)
    return d
def inverse_2x2(A):
    d=det(A)
    if abs(d)<1e-15: return None
    return [[A[1][1]/d,-A[0][1]/d],[-A[1][0]/d,A[0][0]/d]]
def eigenvalues_2x2(A):
    a,b,c,d=A[0][0],A[0][1],A[1][0],A[1][1]
    tr=a+d;dt=a*d-b*c;disc=tr*tr-4*dt
    if disc>=0: sq=math.sqrt(disc);return (tr+sq)/2,(tr-sq)/2
    sq=math.sqrt(-disc);return complex(tr/2,sq/2),complex(tr/2,-sq/2)
def svd_2x2(A):
    a,b,c,d=A[0][0],A[0][1],A[1][0],A[1][1]
    s1_sq=((a*a+b*b+c*c+d*d)+math.sqrt((a*a+b*b-c*c-d*d)**2+4*(a*c+b*d)**2))/2
    s2_sq=((a*a+b*b+c*c+d*d)-math.sqrt((a*a+b*b-c*c-d*d)**2+4*(a*c+b*d)**2))/2
    return math.sqrt(max(s1_sq,0)),math.sqrt(max(s2_sq,0))
def solve(A,b):
    n=len(A);aug=[A[i][:]+[b[i]] for i in range(n)]
    for i in range(n):
        mx=max(range(i,n),key=lambda r:abs(aug[r][i]))
        aug[i],aug[mx]=aug[mx],aug[i]
        for j in range(i+1,n):
            f=aug[j][i]/aug[i][i]
            for k in range(i,n+1): aug[j][k]-=f*aug[i][k]
    x=[0]*n
    for i in range(n-1,-1,-1):
        x[i]=(aug[i][n]-sum(aug[i][j]*x[j] for j in range(i+1,n)))/aug[i][i]
    return x
if __name__=="__main__":
    A=[[3,1],[1,3]];d=det(A);assert d==8
    inv=inverse_2x2(A);print(f"Inverse: {inv}")
    ev=eigenvalues_2x2(A);print(f"Eigenvalues: {ev}")
    s=svd_2x2(A);print(f"Singular values: {s}")
    x=solve([[2,1,1],[4,3,3],[8,7,9]],[1,1,1])
    print(f"Solution: {x}");print("Linear algebra OK")
