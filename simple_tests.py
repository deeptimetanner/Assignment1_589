
import math, cmath
from cubic_solver import solve_cubic, solve_quadratic
from quartic_solver import solve_quartic
from ladder_problem import solve_ladder_problem

def is_close(a, b, tol=1e-8):
    """Check if two complex/real numbers are close"""
    if isinstance(a, complex) or isinstance(b, complex):
        return abs(a - b) < tol
    return abs(a - b) < tol

def verify_roots(coeffs, roots, tol=1e-8):
    """Verify that roots satisfy the polynomial equation"""
    degree = len(coeffs) - 1
    for root in roots:
        value = 0
        for i, coeff in enumerate(coeffs):
            value += coeff * (root ** (degree - i))
        if abs(value) > tol:
            return False
    return True

def test_quadratic():
    print("Testing Quadratic Solver:")
    
    # Test 1: x^2 - 5x + 6 = 0, roots: 2, 3
    roots = solve_quadratic(1, -5, 6)
    expected = [2.0, 3.0]
    print(f"  x^2 - 5x + 6 = 0: {roots}")
    success = verify_roots([1, -5, 6], roots)
    print(f"  Verification: {'PASS' if success else 'FAIL'}")
    
    # Test 2: x^2 + 1 = 0, roots: ±i
    roots = solve_quadratic(1, 0, 1)
    print(f"  x^2 + 1 = 0: {roots}")
    success = verify_roots([1, 0, 1], roots)
    print(f"  Verification: {'PASS' if success else 'FAIL'}")

def test_cubic():
    print("\nTesting Cubic Solver:")
    
    # Test 1: x^3 - 6x^2 + 11x - 6 = 0, roots: 1, 2, 3
    roots = solve_cubic(1, -6, 11, -6)
    print(f"  x^3 - 6x^2 + 11x - 6 = 0: {roots}")
    success = verify_roots([1, -6, 11, -6], roots)
    print(f"  Verification: {'PASS' if success else 'FAIL'}")
    
    # Test 2: x^3 - 1 = 0, roots: 1, ω, ω²
    roots = solve_cubic(1, 0, 0, -1)
    print(f"  x^3 - 1 = 0: {roots}")
    success = verify_roots([1, 0, 0, -1], roots)
    print(f"  Verification: {'PASS' if success else 'FAIL'}")

def test_quartic():
    print("\nTesting Quartic Solver:")
    
    # Test 1: (x-1)(x-2)(x-3)(x-4) = x^4 - 10x^3 + 35x^2 - 50x + 24
    roots = solve_quartic(1, -10, 35, -50, 24)
    print(f"  (x-1)(x-2)(x-3)(x-4) = 0: {roots}")
    success = verify_roots([1, -10, 35, -50, 24], roots)
    print(f"  Verification: {'PASS' if success else 'FAIL'}")
    
    # Test 2: x^4 - 1 = 0, roots: ±1, ±i
    roots = solve_quartic(1, 0, 0, 0, -1)
    print(f"  x^4 - 1 = 0: {roots}")
    success = verify_roots([1, 0, 0, 0, -1], roots)
    print(f"  Verification: {'PASS' if success else 'FAIL'}")

def test_ladder():
    print("\nTesting Ladder Problem:")
    
    # Known test case
    L1, L2, h = 40, 30, 10
    solutions = solve_ladder_problem(L1, L2, h)
    print(f"  L1={L1}, L2={L2}, h={h}: {solutions}")
    
    # Check if we get approximately the correct answer (around 26)
    expected_w = 26.0
    if solutions:
        success = any(abs(w - expected_w) < 1.0 for w in solutions)
        print(f"  Expected w ≈ {expected_w}: {'PASS' if success else 'FAIL'}")
    else:
        print(f"  Expected w ≈ {expected_w}: FAIL (no solutions found)")

if __name__ == "__main__":
    test_quadratic()
    test_cubic() 
    test_quartic()
    test_ladder()