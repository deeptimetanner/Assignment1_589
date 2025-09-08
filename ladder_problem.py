"""
Solution to the Two-Ladder Problem using quartic equation solver
"""
import math
from quartic_solver import solve_quartic

def solve_ladder_problem(L1, L2, h):
    """
    Solve the two-ladder problem.
    
    Args:
        L1, L2: ladder lengths
        h: crossing height
        
    Returns:
        w: alley width
    """
    # First, try the direct algebraic approach
    valid_solutions = solve_ladder_algebraic(L1, L2, h)
    
    # If algebraic method fails, use numerical bisection as fallback
    if not valid_solutions:
        numerical_solution = solve_ladder_bisection(L1, L2, h)
        if numerical_solution is not None:
            valid_solutions = [numerical_solution]
            print(f"Note: Used numerical bisection method. Algebraic methods had numerical issues.")
    
    return valid_solutions


def solve_ladder_algebraic(L1, L2, h):
    """Solve using algebraic quartic methods"""
    # Convert to dimensionless form
    a = (L1/h)**2  # a = (L1/h)^2
    b = (L2/h)**2  # b = (L2/h)^2
    
    # From Exercise 3.23, the quartic equation for u = sqrt(a-ξ) is:
    # u^4 - 2*u^3 - (a-b)*u^2 + 2*(a-b)*u - (a-b) = 0
    
    # Coefficients for quartic equation in u
    coeff_u4 = 1
    coeff_u3 = -2
    coeff_u2 = -(a-b)
    coeff_u1 = 2*(a-b)
    coeff_u0 = -(a-b)
    
    # Solve quartic equation
    u_roots = solve_quartic(coeff_u4, coeff_u3, coeff_u2, coeff_u1, coeff_u0)
    
    # Convert back to ξ and then to w
    valid_solutions = []
    
    for u in u_roots:
        if isinstance(u, complex):
            if abs(u.imag) < 1e-6:  # Essentially real
                u = u.real
            else:
                continue  # Skip complex roots
        
        if u > 0:  # u must be positive
            xi = a - u**2
            if xi > 0:  # ξ must be positive for physical solution
                w = h * math.sqrt(xi)
                valid_solutions.append(w)
    
    return valid_solutions


def solve_ladder_bisection(L1, L2, h, tolerance=1e-10):
    """
    Solve ladder problem using bisection method as fallback.
    This ensures we always get a solution when one exists.
    """
    def ladder_equation(w):
        """The ladder equation that should equal zero at the solution"""
        if w >= min(L1, L2):
            return float('inf')
        term1 = h / math.sqrt(L1*L1 - w*w)  
        term2 = h / math.sqrt(L2*L2 - w*w)
        return term1 + term2 - 1.0
    
    # Set up bisection bounds
    a = 0.0
    b = min(L1, L2) - 0.1  # Slightly less than minimum ladder length
    
    # Check if there's a sign change (solution exists)
    fa = ladder_equation(a)
    fb = ladder_equation(b)
    
    if fa * fb > 0:
        return None  # No solution in interval
    
    # Bisection method
    while (b - a) / 2.0 >= tolerance:
        midpoint = (a + b) / 2.0
        fmid = ladder_equation(midpoint)
        
        if abs(fmid) < tolerance:
            return midpoint
        elif fmid * fa < 0:
            b = midpoint
            fb = fmid
        else:
            a = midpoint
            fa = fmid
    
    return (a + b) / 2.0

def main():
    L1, L2, h = 40, 30, 10
    solutions = solve_ladder_problem(L1, L2, h)
    print(f"Two-ladder problem: L1={L1}, L2={L2}, h={h}")
    print(f"Alley width solutions: {solutions}")

if __name__ == "__main__":
    main()