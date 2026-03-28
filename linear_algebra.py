#!/usr/bin/env python3
"""Linear algebra — vector ops, dot, cross, norms, projections."""
import sys, math, json
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def cross(a,b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def norm(a): return math.sqrt(sum(x*x for x in a))
def normalize(a): n=norm(a); return [x/n for x in a] if n else a
def angle(a,b): return math.degrees(math.acos(min(1,max(-1,dot(a,b)/(norm(a)*norm(b))))))
def project(a,b): s=dot(a,b)/dot(b,b); return [s*x for x in b]
def add(a,b): return [x+y for x,y in zip(a,b)]
def scale(a,s): return [x*s for x in a]
def cli():
    if len(sys.argv)<2: print("Usage: linear_algebra <cmd> [vec1_json] [vec2_json]"); sys.exit(1)
    cmd=sys.argv[1]; a=json.loads(sys.argv[2]) if len(sys.argv)>2 else [1,2,3]
    b=json.loads(sys.argv[3]) if len(sys.argv)>3 else [4,5,6]
    if cmd=="dot": print(f"  {a} · {b} = {dot(a,b)}")
    elif cmd=="cross": print(f"  {a} × {b} = {cross(a,b)}")
    elif cmd=="norm": print(f"  |{a}| = {norm(a):.4f}")
    elif cmd=="angle": print(f"  angle = {angle(a,b):.2f}°")
    elif cmd=="project": print(f"  proj = {[round(x,4) for x in project(a,b)]}")
    elif cmd=="normalize": print(f"  unit = {[round(x,4) for x in normalize(a)]}")
if __name__=="__main__": cli()
