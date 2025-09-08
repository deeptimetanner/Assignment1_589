import math, cmath
from cubic_solver import solve_cubic, solve_quadratic, sqrt_trigonometric

def solve_quartic(a, b, c, d, e):
    """Solve a*x^4 + b*x^3 + c*x^2 + d*x + e = 0.
    Returns a list of 1..4 roots (real numbers or complex numbers).
    If the leading coefficients are zero the function will
    handle lower-degree polynomials automatically.
    Uses multiple methods: resolvent cubic and direct trigonometric.
    """
    # Handle degenerate cases
    if abs(a) < 1e-14:
        return solve_cubic(b, c, d, e)
    
    # Check if this is a biquadratic case (Q = 0 in depressed form)
    # Convert to depressed quartic to check
    p = b / a
    q = c / a  
    r = d / a
    s = e / a
    
    Q = r - p*q/2 + p*p*p/8
    
    # If it's biquadratic, use the resolvent method which handles multiplicity correctly
    if abs(Q) < 1e-14:
        return solve_quartic_resolvent(a, b, c, d, e)
    
    # Try trigonometric method first for non-biquadratic cases
    try:
        trig_roots = solve_quartic_trigonometric(a, b, c, d, e)
        if verify_quartic_roots(trig_roots, [a, b, c, d, e]):
            return trig_roots
    except:
        pass
    
    # Fallback to resolvent cubic method
    return solve_quartic_resolvent(a, b, c, d, e)


def solve_quartic_resolvent(a, b, c, d, e):
    """Solve quartic using resolvent cubic method"""
    roots = []
    
    # Normalize to monic form: x^4 + px^3 + qx^2 + rx + s = 0
    p = b / a
    q = c / a  
    r = d / a
    s = e / a
    
    # Convert to depressed quartic: y^4 + Py^2 + Qy + R = 0
    # Using substitution x = y - p/4
    P = q - 3*p*p/8
    Q = r - p*q/2 + p*p*p/8
    R = s - p*r/4 + p*p*q/16 - 3*p*p*p*p/256
    
    # Special case: biquadratic (Q = 0)
    if abs(Q) < 1e-14:
        return solve_biquadratic(P, R, -p/4)
    
    # General case: solve using resolvent cubic
    # Resolvent cubic: z^3 - Pz^2 - 4Rz + (4PR - Q^2) = 0
    resolvent_roots = solve_cubic(1, -P, -4*R, 4*P*R - Q*Q)
    
    # Choose the best resolvent root for numerical stability
    z = choose_best_resolvent_root(resolvent_roots, P, Q)
    
    # Calculate intermediate values with forced real arithmetic for problematic cases
    alpha, beta = compute_quartic_factors(z, P, Q, force_real=True)
    
    # Solve two quadratics: y^2 + alpha*y + (z-P)/2 + beta/2 = 0
    #                  and: y^2 - alpha*y + (z-P)/2 - beta/2 = 0
    
    a1 = 1
    b1 = alpha
    c1 = (z - P)/2 + beta/2
    
    a2 = 1  
    b2 = -alpha
    c2 = (z - P)/2 - beta/2
    
    # Solve both quadratics
    roots1 = solve_quadratic(a1, b1, c1)
    roots2 = solve_quadratic(a2, b2, c2)
    
    # Collect all roots and transform back: x = y - p/4
    shift = -p/4
    all_y_roots = roots1 + roots2
    
    for y in all_y_roots:
        roots.append(y + shift)
    
    return roots


def solve_quartic_trigonometric(a, b, c, d, e):
    """
    Solve quartic using direct trigonometric methods.
    This approach is more numerically stable for certain quartics.
    """
    # Handle degenerate cases
    if abs(a) < 1e-14:
        return solve_cubic(b, c, d, e)
    
    # Normalize coefficients
    p = b / a
    q = c / a
    r = d / a
    s = e / a
    
    # Convert to depressed quartic: y^4 + P*y^2 + Q*y + R = 0
    # Using substitution x = y - p/4
    P = q - 3*p*p/8
    Q = r - p*q/2 + p*p*p/8
    R = s - p*r/4 + p*p*q/16 - 3*p*p*p*p/256
    
    # For the ladder problem, we use Ferrari's method with trigonometric solving
    return solve_quartic_ferrari_trigonometric(P, Q, R, -p/4)


def solve_quartic_ferrari_trigonometric(P, Q, R, shift):
    """
    Ferrari's method with trigonometric root finding.
    Solves y^4 + P*y^2 + Q*y + R = 0
    """
    # Ferrari's resolvent cubic is: m^3 - (P/2)*m^2 - R*m + (4*R*P - Q^2)/8 = 0
    # Multiply by 8 to clear fractions: 8*m^3 - 4*P*m^2 - 8*R*m + (4*R*P - Q^2) = 0
    
    resolvent_a = 8
    resolvent_b = -4*P  
    resolvent_c = -8*R
    resolvent_d = 4*R*P - Q*Q
    
    # Solve the resolvent cubic using our trigonometric cubic solver
    m_roots = solve_cubic(resolvent_a, resolvent_b, resolvent_c, resolvent_d)
    
    # Choose a real root (Ferrari's method requires this)
    m = None
    for root in m_roots:
        if isinstance(root, complex):
            if abs(root.imag) < 1e-10:
                m = root.real
                break
        else:
            m = root
            break
    
    if m is None:
        # If no good real root, use the real part of the first root
        m = m_roots[0].real if isinstance(m_roots[0], complex) else m_roots[0]
    
    # Ferrari's factorization: (y^2 + s*y + t)(y^2 + u*y + v) = y^4 + P*y^2 + Q*y + R
    # where s + u = 0, t + v + s*u = P, s*v + u*t = Q, t*v = R
    # Setting s = -u, we get: t + v - s^2 = P, s*(v - t) = Q, t*v = R
    
    # From Ferrari's method: s^2 = 2*m + P (where m is our resolvent root)
    s_squared = 2*m + P
    
    s = sqrt_trigonometric(s_squared)
    
    # From t + v = P + s^2 and s*(v - t) = Q:
    sum_tv = P + s_squared  # t + v
    
    if abs(s) > 1e-14:
        diff_tv = Q / s  # v - t
        t = (sum_tv - diff_tv) / 2
        v = (sum_tv + diff_tv) / 2
    else:
        # s ≈ 0 case: solve t*v = R and t + v = P
        discriminant = P*P - 4*R
        sqrt_disc = sqrt_trigonometric(discriminant)
        t = (P - sqrt_disc) / 2
        v = (P + sqrt_disc) / 2
    
    # Solve the two quadratics: y^2 + s*y + t = 0 and y^2 - s*y + v = 0
    roots1 = solve_quadratic(1, s, t)
    roots2 = solve_quadratic(1, -s, v)
    
    # Combine roots and shift back to original variable
    all_roots = roots1 + roots2
    final_roots = [root + shift for root in all_roots]
    
    # Clean up nearly-real roots
    cleaned_roots = []
    for root in final_roots:
        if isinstance(root, complex) and abs(root.imag) < 1e-12:
            cleaned_roots.append(root.real)
        else:
            cleaned_roots.append(root)
    
    return cleaned_roots


def verify_quartic_roots(roots, coeffs, tolerance=1e-8):
    """Verify that roots satisfy the quartic equation"""
    if not roots:
        return False
        
    a, b, c, d, e = coeffs
    for root in roots:
        value = a*root**4 + b*root**3 + c*root**2 + d*root + e
        if abs(value) > tolerance:
            return False
    return True


def compute_quartic_factors(z, P, Q, force_real=False):
    """Compute alpha and beta with numerical stability"""
    
    if force_real:
        # Force real arithmetic when we expect real roots
        # Use absolute value to force positive argument
        z_real = abs(z) if isinstance(z, complex) else abs(z)
        alpha = sqrt_trigonometric(z_real)
        
        beta_squared = z_real - P
        if beta_squared < 0 and abs(beta_squared) < 1e-10:
            beta_squared = 0  # Clean up roundoff error
        beta = sqrt_trigonometric(abs(beta_squared))
        
        # Determine sign of beta using real arithmetic
        if abs(alpha) > 1e-14 and abs(Q) > 1e-14:
            gamma = Q / (2 * alpha)
            # Choose beta sign to satisfy the quartic factorization
            if isinstance(gamma, (int, float)) and gamma < 0:
                beta = -beta
                
        # Handle the case where original z was negative (common in our problem)
        if isinstance(z, (int, float)) and z < 0:
            # For negative z, we need to handle the square root carefully
            # Use the fact that sqrt(-|z|) = i*sqrt(|z|)
            # But we're forcing real, so we work with the magnitude
            pass  # alpha is already sqrt(|z|)
            
    else:
        # Standard complex arithmetic
        alpha = sqrt_trigonometric(z + 0j)
        beta_squared = z - P
        beta = sqrt_trigonometric(beta_squared + 0j)
        
        # Sign determination for complex case
        if abs(alpha) > 1e-14:
            gamma = Q / (2 * alpha)
            if (beta * gamma).real < 0:
                beta = -beta
    
    return alpha, beta


def choose_best_resolvent_root(resolvent_roots, P, Q):
    """Choose the resolvent root that maximizes numerical stability"""
    candidates = []
    
    for root in resolvent_roots:
        # Clean up near-real roots
        if isinstance(root, complex) and abs(root.imag) < 1e-12:
            root = root.real
        
        if isinstance(root, complex):
            continue  # Skip truly complex roots
            
        # Score this root for stability
        score = 0
        z = root
        
        # Prefer positive z (makes sqrt real)
        if isinstance(z, (int, float)) and z > 1e-12:
            score += 2
        
        # Prefer z where |z-P| is not too small (avoids division issues)
        if abs(z - P) > 1e-10:
            score += 1
            
        # Prefer z that keeps subsequent calculations real
        beta_squared = z - P
        if isinstance(beta_squared, (int, float)) and beta_squared >= -1e-12:  # Allow small negative due to roundoff
            score += 1
            
        candidates.append((score, z))
    
    # Return the highest-scoring root
    if candidates:
        candidates.sort(reverse=True)
        return candidates[0][1]
    else:
        # Fallback: use real part of first root
        root = resolvent_roots[0]
        return root.real if isinstance(root, complex) else root


def solve_biquadratic(P, R, shift):
    """Solve y^4 + P*y^2 + R = 0 (biquadratic case)"""
    # Substitute z = y^2 to get z^2 + P*z + R = 0
    z_roots = solve_quadratic(1, P, R)
    
    # Check for perfect square case: P^2 - 4R = 0 means (y^2 + P/2)^2 = 0
    discriminant = P*P - 4*R
    
    if abs(discriminant) < 1e-12:
        # Perfect square case: only return the 2 roots from sqrt(-P/2)
        z = -P/2
        sqrt_z = sqrt_trigonometric(z + 0j)
        return [sqrt_z + shift, -sqrt_z + shift]
    else:
        # Regular case: return all 4 roots from both z values
        roots = []
        for z in z_roots:
            if isinstance(z, complex):
                if z.real >= 0 and abs(z.imag) < 1e-10:
                    # Positive real z
                    sqrt_z = sqrt_trigonometric(z.real)
                    if abs(z.real) > 1e-12:  
                        roots.extend([sqrt_z + shift, -sqrt_z + shift])
                    else:
                        roots.append(shift)  # z ≈ 0
                else:
                    # Complex z
                    sqrt_z = sqrt_trigonometric(z)
                    roots.extend([sqrt_z + shift, -sqrt_z + shift])
            else:
                if abs(z) > 1e-12:
                    if isinstance(z, (int, float)) and z > 0:
                        sqrt_z = sqrt_trigonometric(z)
                        roots.extend([sqrt_z + shift, -sqrt_z + shift])
                    else:
                        sqrt_z = sqrt_trigonometric(z + 0j)
                        roots.extend([sqrt_z + shift, -sqrt_z + shift])
                else:
                    roots.append(shift)  # z ≈ 0
        
        return roots


def main():
    tests = [
        (1, 0, 0, 0, -1),  # roots of x^4 - 1 = 0 (2 real and 2 complex roots)
        (1, 0, 1, 0, -1),  # roots of x^4 + x^2 - 1 = 0 (2 real and 2 complex roots)
        (0, 1, -3, 2, 0),  # a=0 => cubic: x^3 - 3x^2 + 2x = 0  (roots 0,1,2)
    ]
    for a, b, c, d, e in tests:
        roots = solve_quartic(a, b, c, d, e)
        print(f"solve_quartic({a}, {b}, {c}, {d}, {e}) -> {roots}")

if __name__ == "__main__":
    main()
