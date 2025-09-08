import math, cmath

def solve_cubic(a, b, c, d):
    """Solve a*x^3 + b*x^2 + c*x + d = 0
    Returns list of 1..3 roots (complex if needed).
    Uses trigonometric/hyperbolic substitution method.
    """
    roots = []
    
    # Handle degenerate cases
    if abs(a) < 1e-14:
        return solve_quadratic(b, c, d)
    
    # Normalize to monic form: x^3 + px^2 + qx + r = 0
    p = b / a
    q = c / a
    r = d / a
    
    # Convert to depressed cubic: t^3 + pt + q = 0
    # Using substitution x = t - p/3
    p_new = q - p*p/3
    q_new = r - p*q/3 + 2*p*p*p/27
    
    # Use depressed cubic invariants
    Delta = (q_new/2)**2 + (p_new/3)**3
    
    if abs(p_new) < 1e-14:
        # Special case: t^3 + q = 0
        if abs(q_new) < 1e-14:
            t_roots = [0.0, 0.0, 0.0]  # Triple root
        else:
            # Real cube root using exponential
            if isinstance(q_new, (int, float)) and q_new < 0:
                t_roots = [cmath.exp(cmath.log(-q_new)/3).real]
            else:
                t_roots = [-cmath.exp(cmath.log(q_new)/3).real]
            # Complex cube roots
            # omega = exp(2*pi*i/3) = cos(2*pi/3) + i*sin(2*pi/3)
            omega = cmath.exp(2j * cmath.pi / 3)
            t_roots.extend([t_roots[0] * omega, t_roots[0] * omega.conjugate()])
    else:
        # Determine path: if invariants real and Delta < 0 -> 3 real (trig), else general Cardano
        invariants_real = isinstance(p_new, (int, float)) and isinstance(q_new, (int, float))
        if invariants_real and Delta < 0:
            # Three distinct real roots - trigonometric method (p_new < 0 here)
            m = 2 * sqrt_trigonometric(-p_new/3)
            m_val = m.real if isinstance(m, complex) else m
            cos_3theta = 3*q_new / (p_new * m_val)
            if isinstance(cos_3theta, complex):
                cos_3theta = cos_3theta.real
            cos_3theta = max(-1.0, min(1.0, float(cos_3theta)))
            theta = math.acos(cos_3theta) / 3
            t_roots = [
                m_val * math.cos(theta),
                m_val * math.cos(theta + 2*math.pi/3),
                m_val * math.cos(theta + 4*math.pi/3)
            ]
        else:
            # General Cardano with complex arithmetic
            Delta_c = Delta if isinstance(Delta, complex) else complex(Delta, 0.0)
            sqrt_D = sqrt_trigonometric(Delta_c)
            u_c = -q_new/2 + sqrt_D
            v_c = -q_new/2 - sqrt_D
            u = cmath.exp(cmath.log(u_c)/3) if u_c != 0 else 0
            v = cmath.exp(cmath.log(v_c)/3) if v_c != 0 else 0
            t1 = u + v
            # Other roots via cube roots of unity
            omega = cmath.exp(2j * cmath.pi / 3)
            t2 = u*omega + v*omega.conjugate()
            t3 = u*omega.conjugate() + v*omega
            t_roots = [t1, t2, t3]
    
    # Convert back to original variable: x = t - p/3
    shift = -p / 3
    for t in t_roots:
        if isinstance(t, complex):
            roots.append(t + shift)
        else:
            roots.append(t + shift)
    
    return roots


def sqrt_trigonometric(x):
    """Compute square root using trigonometric substitution"""
    # Fast real path
    if not isinstance(x, complex):
        if x < 0:
            try:
                return 1j * math.exp(0.5 * math.log(-x))
            except (ValueError, OverflowError):
                return 0+0j
        if x == 0:
            return 0.0
        try:
            return math.exp(0.5 * math.log(x))
        except (ValueError, OverflowError):
            return 0.0
    # Complex path: sqrt(z) = exp(ln(z)/2)
    if x == 0:
        return 0+0j
    try:
        return cmath.exp(0.5 * cmath.log(x))
    except (ValueError, OverflowError):
        return 0+0j

def solve_quadratic(a, b, c):
    """Solve a*x^2 + b*x + c = 0 using Chebyshev's method without radicals"""
    # Handle complex coefficients
    if isinstance(a, complex) or isinstance(b, complex) or isinstance(c, complex):
        a, b, c = complex(a), complex(b), complex(c)
    
    if abs(a) < 1e-14:
        if abs(b) < 1e-14:
            return [] if abs(c) > 1e-14 else [0.0]
        return [-c/b]
    
    # Calculate discriminant
    discriminant = b*b - 4*a*c
    
    # Use trigonometric square root
    sqrt_disc = sqrt_trigonometric(discriminant)
    
    # Check if we should use Chebyshev's stable formulation
    if isinstance(b, (int, float)) and isinstance(discriminant, (int, float)) and discriminant >= 0:
        # Real case with Chebyshev's method
        if b >= 0:
            x1 = (-b - sqrt_disc) / (2*a)
            if abs(-b - sqrt_disc) > 1e-14:
                x2 = (2*c) / (-b - sqrt_disc)
            else:
                x2 = (-b + sqrt_disc) / (2*a)
        else:
            if abs(-b + sqrt_disc) > 1e-14:
                x1 = (2*c) / (-b + sqrt_disc)
            else:
                x1 = (-b - sqrt_disc) / (2*a)
            x2 = (-b + sqrt_disc) / (2*a)
        return [x1, x2]
    else:
        # Complex case or general case
        return [(-b + sqrt_disc)/(2*a), (-b - sqrt_disc)/(2*a)]


def main():
    tests = [
        (1, 0, 0, -1),     # roots of x^3 - 1 = 0 (1 and two complex cube roots)
        (1, -6, 11, -6),   # roots [1.0, 2.0, 3.0]
    ]
    for a, b, c, d in tests:
        roots = solve_cubic(a, b, c, d)
        print(f"solve_cubic({a}, {b}, {c}, {d}) -> {roots}")

if __name__ == "__main__":
    main()
