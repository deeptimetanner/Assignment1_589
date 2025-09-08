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
    
    # Discriminant
    discriminant = -(4*p_new*p_new*p_new + 27*q_new*q_new)
    
    if abs(p_new) < 1e-14:
        # Special case: t^3 + q = 0
        if abs(q_new) < 1e-14:
            t_roots = [0.0, 0.0, 0.0]  # Triple root
        else:
            # Real cube root using exponential
            if q_new < 0:
                t_roots = [cmath.exp(cmath.log(-q_new)/3).real]
            else:
                t_roots = [-cmath.exp(cmath.log(q_new)/3).real]
            # Complex cube roots
            # omega = exp(2*pi*i/3) = cos(2*pi/3) + i*sin(2*pi/3)
            omega = cmath.exp(2j * cmath.pi / 3)
            t_roots.extend([t_roots[0] * omega, t_roots[0] * omega.conjugate()])
    elif discriminant > 0:
        # Three distinct real roots - use trigonometric method
        # t^3 + pt + q = 0 with p < 0
        m = 2 * sqrt_trigonometric(-p_new/3)
        # cos(3*theta) = 3*q_new / (p_new * m)
        cos_3theta = 3*q_new / (p_new * m)
        
        # Clamp to valid range for acos
        cos_3theta = max(-1.0, min(1.0, cos_3theta))
        
        theta = math.acos(cos_3theta) / 3
        
        t_roots = [
            m * math.cos(theta),
            m * math.cos(theta + 2*math.pi/3),
            m * math.cos(theta + 4*math.pi/3)
        ]
    else:
        # One real root, two complex roots - use hyperbolic method
        if abs(p_new) > 1e-14:
            m = 2 * sqrt_trigonometric(abs(p_new)/3)
            cosh_3u = abs(3*q_new / (p_new * m))
            u = math.acosh(cosh_3u) / 3
            
            if q_new > 0:
                t_real = -m * math.cosh(u)
            else:
                t_real = m * math.cosh(u)
                
            t_roots = [t_real]
            
            # Complex roots using the cubic formula
            # For complex roots when discriminant < 0
            sqrt_disc = sqrt_trigonometric(discriminant + 0j)
            sqrt_27 = sqrt_trigonometric(27+0j)
            alpha = (-q_new + sqrt_disc/sqrt_27) / 2
            beta = (-q_new - sqrt_disc/sqrt_27) / 2
            
            # Cube root using exponential: x^(1/3) = exp(ln(x)/3)
            alpha_cbrt = cmath.exp(cmath.log(alpha)/3) if alpha != 0 else 0
            beta_cbrt = cmath.exp(cmath.log(beta)/3) if beta != 0 else 0
            
            # omega = exp(2*pi*i/3) = cos(2*pi/3) + i*sin(2*pi/3)
            omega = cmath.exp(2j * cmath.pi / 3)
            
            t_roots.extend([
                alpha_cbrt * omega + beta_cbrt * omega.conjugate(),
                alpha_cbrt * omega.conjugate() + beta_cbrt * omega
            ])
    
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
    if isinstance(x, complex):
        # For complex numbers, use exponential form
        if x == 0:
            return 0+0j
        # sqrt(z) = exp(ln(z)/2)
        try:
            return cmath.exp(0.5 * cmath.log(x))
        except (ValueError, OverflowError):
            return 0+0j
    else:
        if x < 0:
            # Return complex result for negative real numbers
            try:
                return 1j * cmath.exp(0.5 * cmath.log(-x)).real
            except (ValueError, OverflowError):
                return 0+0j
        elif x == 0:
            return 0.0
        else:
            # Use exp(ln(x)/2) instead of sqrt
            try:
                result = cmath.exp(0.5 * cmath.log(x))
                return result.real if abs(result.imag) < 1e-14 else result
            except (ValueError, OverflowError):
                return 0.0

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
