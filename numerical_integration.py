#!/usr/bin/env python3
"""Numerical integration — trapezoid, Simpson, Gauss-Legendre."""
import math

def trapezoid(f, a, b, n=1000):
    h = (b-a)/n; s = (f(a)+f(b))/2
    for i in range(1, n): s += f(a+i*h)
    return s*h

def simpson(f, a, b, n=1000):
    if n%2: n+=1
    h = (b-a)/n; s = f(a)+f(b)
    for i in range(1,n,2): s += 4*f(a+i*h)
    for i in range(2,n,2): s += 2*f(a+i*h)
    return s*h/3

def gauss_legendre(f, a, b, n=5):
    """n-point Gauss-Legendre quadrature (hardcoded for n=2,3,5)."""
    nodes = {2: [(-0.5773502692, 1.0), (0.5773502692, 1.0)],
             3: [(-0.7745966692, 0.5555555556), (0.0, 0.8888888889), (0.7745966692, 0.5555555556)],
             5: [(-0.9061798459, 0.2369268851), (-0.5384693101, 0.4786286705), (0.0, 0.5688888889),
                 (0.5384693101, 0.4786286705), (0.9061798459, 0.2369268851)]}
    pts = nodes.get(n, nodes[5])
    mid = (a+b)/2; half = (b-a)/2
    return half * sum(w * f(mid + half*x) for x, w in pts)

def main():
    f = lambda x: x**2
    print(f"∫x² dx [0,1] = trap:{trapezoid(f,0,1):.6f} simp:{simpson(f,0,1):.6f} GL:{gauss_legendre(f,0,1):.6f}")

if __name__ == "__main__": main()
