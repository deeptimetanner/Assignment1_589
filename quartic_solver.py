import math, cmath
from cubic_solver import solve_cubic, solve_quadratic, sqrt_trigonometric

def get_distinct_roots(roots, tolerance=1e-10):
    """
    Extract distinct roots with their multiplicities from a list of roots.
    Returns list of (root, multiplicity) tuples.
    """
    from collections import defaultdict
    
    distinct_roots = defaultdict(int)
    
    for root in roots:
        # Find if this root is close to any existing distinct root
        found = False
        for distinct_root in distinct_roots:
            if abs(root - distinct_root) < tolerance:
                distinct_roots[distinct_root] += 1
                found = True
                break
        
        if not found:
            distinct_roots[root] += 1
    
    return [(root, mult) for root, mult in distinct_roots.items()]

def solve_quartic(a, b, c, d, e):
    """Solve a*x^4 + b*x^3 + c*x^2 + d*x + e = 0.
    Returns a list of 4 roots (real numbers or complex numbers) with multiplicity.
    If the leading coefficients are zero the function will
    handle lower-degree polynomials automatically.
    Uses multiple methods: resolvent cubic and direct trigonometric.
    """
    # Handle degenerate cases
    if abs(a) < 1e-14:
        return solve_cubic(b, c, d, e)
    
    # Pre-check: handle biquadratic pattern early using depressed form
    try:
        p = b / a
        q = c / a
        r = d / a
        s = e / a
        P = q - 3*p*p/8
        Q = r - p*q/2 + p*p*p/8
        R = s - p*r/4 + p*p*q/16 - 3*p*p*p*p/256
        if abs(Q) < 1e-9 * (1.0 + abs(P) + abs(R)):
            return finalize_roots(solve_biquadratic(P, R, -p/4))
    except Exception:
        pass

    # Prefer resolvent method; if it produces poor residuals, try trig variant
    try:
        res_roots = solve_quartic_resolvent(a, b, c, d, e)
        if verify_quartic_roots(res_roots, [a, b, c, d, e]) and len(res_roots) == 4:
            return finalize_roots(res_roots)
    except:
        pass

    try:
        trig_roots = solve_quartic_trigonometric(a, b, c, d, e)
        if verify_quartic_roots(trig_roots, [a, b, c, d, e]) and len(trig_roots) == 4:
            return finalize_roots(trig_roots)
    except:
        pass

    # As a last resort, return whatever resolvent produced (possibly fewer roots)
    return finalize_roots(res_roots if 'res_roots' in locals() else [])


def solve_quartic_resolvent(a, b, c, d, e):
    """Solve quartic using resolvent cubic method"""
    # Normalize to monic form: x^4 + p x^3 + q x^2 + r x + s = 0
    p = b / a
    q = c / a
    r = d / a
    s = e / a

    # Depressed quartic y^4 + P y^2 + Q y + R = 0 via x = y - p/4
    P = q - 3.0*(p*p)/8.0
    Q = r - 0.5*p*q + (p*p*p)/8.0
    R = s - 0.25*p*r + (p*p*q)/16.0 - 3.0*(p**4)/256.0

    # Biquadratic shortcut
    if abs(Q) < 1e-12 * (1.0 + abs(P) + abs(R)):
        return solve_biquadratic(P, R, -p/4.0)

    # Standard Ferrari resolvent cubic: z^3 - P z^2 - 4 R z + (4 P R - Q^2) = 0
    z_candidates = solve_cubic(1.0, -P, -4.0*R, 4.0*P*R - Q*Q)

    best_set = None
    best_res = float('inf')

    def residual_sum(roots_set):
        # Compute polynomial via Horner for slight speed
        sres = 0.0
        for xr in roots_set:
            p = ((a*xr + b)*xr + c)*xr + d
            p = p*xr + e
            sres += abs(p)
        return sres

    # From each resolvent root, build candidate and early-verify
    for z0 in z_candidates:
        # Use nearly-real z0 only
        if isinstance(z0, complex):
            if abs(z0.imag) > 1e-10:
                continue
            z0 = z0.real

        # S = sqrt(z0 - P) per Ferrari factorization
        S = sqrt_trigonometric(z0 - P)
        if abs(S) < 1e-14 and abs(Q) > 1e-14:
            # Avoid division by zero when Q != 0
            continue

        if abs(S) > 1e-14:
            e_const = (z0 - Q / S) / 2.0
            f_const = (z0 + Q / S) / 2.0
        else:
            # Q ~ 0 handled by biquadratic path earlier; split evenly
            e_const = z0 / 2.0
            f_const = z0 / 2.0

        # Solve the two quadratics in y
        y_roots_1 = solve_quadratic(1.0, S, e_const)
        y_roots_2 = solve_quadratic(1.0, -S, f_const)
        x_roots = [y + (-p/4.0) for y in (y_roots_1 + y_roots_2)]
        # Early accept if they verify
        if verify_quartic_roots(x_roots, [a, b, c, d, e], tolerance=1e-8):
            return x_roots
        # Track best residual
        res = residual_sum(x_roots)
        if res < best_res:
            best_res = res
            best_set = x_roots

    # Fallback: try trig variant only if none verified
    trig_candidates = solve_quartic_ferrari_trigonometric(P, Q, R, -p/4.0)
    if trig_candidates and verify_quartic_roots(trig_candidates, [a, b, c, d, e], tolerance=1e-8):
        return trig_candidates
    # Choose the best of what we built
    return best_set if best_set is not None else []


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
    """Compute alpha and beta per Ferrari: alpha = sqrt(2*z - P), beta = Q/alpha."""
    
    z_c = complex(z)
    P_c = complex(P)
    Q_c = complex(Q)

    alpha_sq = 2*z_c - P_c
    alpha = sqrt_trigonometric(alpha_sq)

    if abs(alpha) > 1e-14:
        beta = Q_c / alpha
    else:
        # Fallback if alpha ~ 0: use alternative beta estimate and keep alpha as computed
        beta = sqrt_trigonometric(z_c - P_c)
        if abs(beta) > 1e-14 and abs(Q_c) > 0:
            gamma = Q_c / (2 * (alpha if abs(alpha) > 1e-14 else 1))
            if (beta * gamma.conjugate()).real < 0:
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
    
    roots = []
    discriminant = P*P - 4*R
    
    if abs(discriminant) < 1e-9:
        # Perfect square case: (y^2 + P/2)^2 = 0, double root with multiplicity 2 each
        z = -P/2
        sqrt_z = sqrt_trigonometric(z + 0j)
        
        # Clean up floating point noise in the square root
        if isinstance(sqrt_z, complex):
            real_part = sqrt_z.real if abs(sqrt_z.real) > 1e-12 else 0.0
            imag_part = sqrt_z.imag if abs(sqrt_z.imag) > 1e-12 else 0.0
            sqrt_z_clean = complex(real_part, imag_part)
        else:
            sqrt_z_clean = sqrt_z if abs(sqrt_z) > 1e-12 else 0.0
        
        # Return all 4 roots with multiplicity: each sqrt appears twice
        pos_root = sqrt_z_clean + shift
        neg_root = -sqrt_z_clean + shift
        
        # Ensure clean representation for pure imaginary numbers
        if isinstance(neg_root, complex) and abs(neg_root.real) < 1e-12:
            neg_root = complex(0.0, neg_root.imag)
        if isinstance(pos_root, complex) and abs(pos_root.real) < 1e-12:
            pos_root = complex(0.0, pos_root.imag)
            
        roots = [pos_root, pos_root, neg_root, neg_root]
    else:
        # Regular case: return all 4 roots from both z values
        for z in z_roots:
            if isinstance(z, complex):
                if z.real >= 0 and abs(z.imag) < 1e-10:
                    # Positive real z
                    sqrt_z = sqrt_trigonometric(z.real)
                    if abs(z.real) > 1e-12:  
                        roots.extend([sqrt_z + shift, -sqrt_z + shift])
                    else:
                        # z ≈ 0, double root at origin
                        roots.extend([shift, shift])
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
                    # z ≈ 0, double root at origin  
                    roots.extend([shift, shift])
    
    # Clean up all roots to remove floating point noise
    cleaned_roots = []
    for root in roots:
        if isinstance(root, complex):
            real_part = root.real if abs(root.real) > 1e-12 else 0.0
            imag_part = root.imag if abs(root.imag) > 1e-12 else 0.0
            # Handle pure imaginary numbers properly
            if abs(real_part) < 1e-12 and abs(imag_part) > 1e-12:
                cleaned_roots.append(complex(0.0, imag_part))
            elif abs(imag_part) < 1e-12 and abs(real_part) > 1e-12:
                cleaned_roots.append(real_part)  # Return as real number
            elif abs(real_part) < 1e-12 and abs(imag_part) < 1e-12:
                cleaned_roots.append(0.0)  # Return as real zero
            else:
                cleaned_roots.append(complex(real_part, imag_part))
        else:
            cleaned_roots.append(root if abs(root) > 1e-12 else 0.0)
    
    return cleaned_roots


def _snap_close_roots(roots, tol=1e-9):
    """Snap numerically close roots to identical representatives to enforce multiplicities."""
    snapped = []
    reps = []
    for z in roots:
        # normalize near-real
        if isinstance(z, complex) and abs(z.imag) < tol:
            z = z.real
        if isinstance(z, (int, float)) and abs(z) < tol:
            z = 0.0
        # find representative
        found = False
        for i, r in enumerate(reps):
            if abs(z - r) < tol:
                snapped.append(reps[i])
                found = True
                break
        if not found:
            reps.append(z)
            snapped.append(z)
    return snapped


def finalize_roots(roots, tol=1e-12):
    """Convert nearly-real complex numbers to real floats and snap duplicates."""
    snapped = _snap_close_roots(roots, tol=max(tol, 1e-10))
    finalized = []
    for z in snapped:
        if isinstance(z, complex) and abs(z.imag) < tol:
            finalized.append(z.real)
        else:
            finalized.append(z)
    return finalized


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
