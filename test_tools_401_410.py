#!/usr/bin/env python3
"""
Test suite for ChemMCP Tools #401-410: Numerical & Linear Algebra Tools
Tests all 10 new MCP tools with multiple scenarios.
"""

import sys
import math
import numpy as np

# ============================================================
# Helper
# ============================================================
passed = 0
failed = 0
errors = []

def check(name, condition, msg=""):
    global passed, failed
    if condition:
        passed += 1
        print(f"  ✅ PASS: {name}")
    else:
        failed += 1
        errors.append((name, msg))
        print(f"  ❌ FAIL: {name} — {msg}")

def approx(a, b, tol=1e-4):
    """Compare floats with tolerance."""
    if isinstance(a, (list, tuple)):
        return all(approx(x, y, tol) for x, y in zip(a, b))
    return abs(a - b) < tol

def approx_list(a, b, tol=1e-3):
    """Compare lists of floats (sorted by absolute value)."""
    if len(a) != len(b):
        return False
    a_sorted = sorted(abs(x) for x in a)
    b_sorted = sorted(abs(x) for x in b)
    return all(approx(ai, bi, tol) for ai, bi in zip(a_sorted, b_sorted))


# ============================================================
# Test 401: MatrixEigenvalueSolver
# ============================================================
print("\n" + "=" * 60)
print("TEST 401: MatrixEigenvalueSolver")
print("=" * 60)

try:
    from chemmcp.tools import MatrixEigenvalueSolver

    tool = MatrixEigenvalueSolver()
    check("Instantiate", tool is not None)

    # Test 1: Symmetric 3x3 matrix (tridiagonal)
    result = tool.run_code(
        matrix=[[2, -1, 0], [-1, 2, -1], [0, -1, 2]],
        symmetric=True,
    )
    check("Returns dict", isinstance(result, dict))
    check("Has eigenvalues", "eigenvalues" in result)
    check("Matrix size", result["matrix_size"] == 3)
    check("Eigenvalues count", len(result["eigenvalues"]) == 3)
    # Known eigenvalues: 2+√2≈3.414, 2, 2-√2≈0.586
    check("Eigenvalues correct (symmetric 3x3)", approx_list(
        result["eigenvalues"], [3.414214, 2.0, 0.585786]
    ))

    # Test 2: Non-symmetric 2x2: [[4,1],[2,3]] → λ²-7λ+10=0 → λ=5, 2
    result2 = tool.run_code(matrix=[[4, 1], [2, 3]], symmetric=False)
    check("Non-symmetric 2x2 size", result2["matrix_size"] == 2)
    check("Non-symmetric eigenvalues count", len(result2["eigenvalues"]) == 2)
    check("Non-symmetric eigenvalues correct", approx_list(
        result2["eigenvalues"], [5.0, 2.0], tol=1e-3
    ))

    # Test 3: Text mode
    result3 = tool.run_text('{"matrix": [[1,0],[0,2]], "symmetric": true}')
    check("Text mode works", isinstance(result3, dict) and len(result3["eigenvalues"]) == 2)

    # Test 4: Error handling - non-square matrix
    try:
        tool.run_code(matrix=[[1, 2, 3], [4, 5, 6]])
        check("Non-square raises error", False, "Should have raised error")
    except Exception as e:
        check("Non-square raises error", True)

except Exception as e:
    failed += 1
    errors.append(("401-Import/Setup", str(e)))
    print(f"  ❌ ERROR: {e}")


# ============================================================
# Test 402: MatrixDiagonalization
# ============================================================
print("\n" + "=" * 60)
print("TEST 402: MatrixDiagonalization")
print("=" * 60)

try:
    from chemmcp.tools import MatrixDiagonalization

    tool = MatrixDiagonalization()
    check("Instantiate", tool is not None)

    # Test: Hückel-type tridiagonal matrix
    result = tool.run_code(
        matrix=[[0, 1, 0], [1, 0, 1], [0, 1, 0]],
        symmetric=True,
    )
    check("Returns dict", isinstance(result, dict))
    check("Diagonal values present", "diagonal_values" in result)
    check("P matrix present", "transformation_matrix" in result)
    check("Is diagonalizable", result["is_diagonalizable"] == True)
    # Eigenvalues of this matrix: √2, 0, -√2
    check("Diagonal values correct", approx_list(
        result["diagonal_values"], [1.414214, 0, -1.414214], tol=1e-3
    ))

    # Verify A ≈ P D P^{-1}
    P = np.array(result["transformation_matrix"])
    d = np.array(result["diagonal_values"])
    A_orig = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    A_reconstructed = P @ np.diag(d) @ np.linalg.inv(P)
    recon_error = np.max(np.abs(A_reconstructed - A_orig))
    check("Reconstruction error small", recon_error < 1e-6, f"error={recon_error:.2e}")

    # With pinv requested
    result2 = tool.run_code(
        matrix=[[1, 2], [3, 4]],
        symmetric=False,
        return_pinv=True,
    )
    check("Pinv returned", result2["inverse_transformation"] is not None)

except Exception as e:
    failed += 1
    errors.append(("402-Import/Setup", str(e)))
    print(f"  ❌ ERROR: {e}")


# ============================================================
# Test 403: DeterminantCalculator
# ============================================================
print("\n" + "=" * 60)
print("TEST 403: DeterminantCalculator")
print("=" * 60)

try:
    from chemmcp.tools import DeterminantCalculator

    tool = DeterminantCalculator()
    check("Instantiate", tool is not None)

    # Test 2x2
    result = tool.run_code(matrix=[[1, 2], [3, 4]], show_details=False)
    check("2x2 determinant", abs(result["determinant"] - (-2.0)) < 1e-6)
    check("Not singular", result["is_singular"] == False)

    # Test 3x3 with details
    result2 = tool.run_code(
        matrix=[[6, 1, 1], [4, -2, 5], [2, 8, 7]],
        show_details=True,
    )
    check("3x3 determinant", abs(result2["determinant"] - (-306.0)) < 1e-6)
    check("Details provided", result2["details"] is not None and len(result2["details"]) > 0)

    # Singular matrix
    result3 = tool.run_code(matrix=[[1, 2], [2, 4]])
    check("Singular detected", result3["is_singular"] == True)
    check("Singular det near zero", abs(result3["determinant"]) < 1e-10)

    # Identity matrix
    result4 = tool.run_code(matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    check("Identity det = 1", abs(result4["determinant"] - 1.0) < 1e-6)

except Exception as e:
    failed += 1
    errors.append(("403-Import/Setup", str(e)))
    print(f"  ❌ ERROR: {e}")


# ============================================================
# Test 404: MatrixInversion
# ============================================================
print("\n" + "=" * 60)
print("TEST 404: MatrixInversion")
print("=" * 60)

try:
    from chemmcp.tools import MatrixInversion

    tool = MatrixInversion()
    check("Instantiate", tool is not None)

    # Regular inversion
    result = tool.run_code(matrix=[[4, 7], [2, 6]], use_pseudo_inverse=False)
    check("Inverse computed", "inverse" in result and len(result["inverse"]) == 2)
    # Verify A * A^(-1) = I
    A = np.array([[4, 7], [2, 6]])
    A_inv = np.array(result["inverse"])
    product = A @ A_inv
    identity_err = np.max(np.abs(product - np.eye(2)))
    check("A*A_inv ≈ I", identity_err < 1e-6, f"error={identity_err:.2e}")

    # Singular matrix with pseudo-inverse
    result2 = tool.run_code(matrix=[[1, 2], [2, 4]], use_pseudo_inverse=True)
    check("Pseudo-inverse for singular", result2["is_singular"] == True)
    check("Pseudo-inverse shape", len(result2["inverse"]) == 2)

    # Condition number
    result3 = tool.run_code(
        matrix=[[10000, 1], [1, 10000]],
        use_pseudo_inverse=False,
        condition_number=True,
    )
    check("Condition number returned", result3["condition_number"] is not None)
    check("Condition number > 1", result3["condition_number"] > 1.0)

except Exception as e:
    failed += 1
    errors.append(("404-Import/Setup", str(e)))
    print(f"  ❌ ERROR: {e}")


# ============================================================
# Test 405: SvdDecomposition
# ============================================================
print("\n" + "=" * 60)
print("TEST 405: SvdDecomposition")
print("=" * 60)

try:
    from chemmcp.tools import SvdDecomposition

    tool = SvdDecomposition()
    check("Instantiate", tool is not None)

    # Basic SVD test
    result = tool.run_code(
        matrix=[[1, 2, 3], [4, 5, 6]],
        full_matrices=False,
    )
    check("Singular values present", "singular_values" in result)
    check("U matrix present", "u_matrix" in result)
    check("Vt matrix present", "vt_matrix" in result)
    check("Effective rank", result["effective_rank"] >= 1)
    check("SVs descending order",
         result["singular_values"] == sorted(result["singular_values"], reverse=True))

    # Verify reconstruction (relaxed tolerance)
    U = np.array(result["u_matrix"])
    s = np.array(result["singular_values"])
    Vt = np.array(result["vt_matrix"])
    A_orig = np.array([[1, 2, 3], [4, 5, 6]])
    A_reconstructed = U @ np.diag(s) @ Vt
    recon_error = np.max(np.abs(A_reconstructed - A_orig))
    check("SVD reconstruction accurate", recon_error < 1e-5, f"error={recon_error:.2e}")

    # Rank approximation
    result2 = tool.run_code(
        matrix=[[1, 2, 3], [4, 5, 6], [7, 8, 9]],
        rank_approximation=1,
    )
    check("Rank-1 approximation has error", result2["reconstruction_error"] is not None)
    check("Reconstruction error positive", result2["reconstruction_error"] > 0)

    # Square matrix
    result3 = tool.run_code(
        matrix=[[3, 1, 1], [1, 3, 1], [1, 1, 3]],
    )
    check("Square matrix SVD", result3["effective_rank"] == 3)

except Exception as e:
    failed += 1
    errors.append(("405-Import/Setup", str(e)))
    print(f"  ❌ ERROR: {e}")


# ============================================================
# Test 406: LinearSystemSolver
# ============================================================
print("\n" + "=" * 60)
print("TEST 406: LinearSystemSolver")
print("=" * 60)

try:
    from chemmcp.tools import LinearSystemSolver

    tool = LinearSystemSolver()
    check("Instantiate", tool is not None)

    # Square system: 2x + y = 5, x + 3y = 6 → x=1.8, y=1.4
    result = tool.run_code(
        matrix_a=[[2, 1], [1, 3]],
        vector_b=[5, 6],
    )
    check("Square system type", result["system_type"] == "square")
    check("Solution length", len(result["solution"]) == 2)
    check("Solution x correct", approx(result["solution"][0], 1.8, 1e-3))
    check("Solution y correct", approx(result["solution"][1], 1.4, 1e-3))
    check("Residual near zero", result["residual_norm"] < 1e-6)

    # Overdetermined system (least squares): exact solution exists
    result2 = tool.run_code(
        matrix_a=[[3, 1], [1, 2], [4, 3]],
        vector_b=[9, 8, 17],
    )
    check("Overdetermined type", result2["system_type"] == "overdetermined")
    check("LS solution x=2", approx(result2["solution"][0], 2.0, 1e-3))
    check("LS solution y=3", approx(result2["solution"][1], 3.0, 1e-3))

    # Text mode
    result3 = tool.run_text('{"matrix_a": [[1,1],[0,1]], "vector_b": [4,3]}')
    check("Text mode works", len(result3["solution"]) == 2)
    check("Text mode solution correct", approx_list(result3["solution"], [1, 3]))

except Exception as e:
    failed += 1
    errors.append(("406-Import/Setup", str(e)))
    print(f"  ❌ ERROR: {e}")


# ============================================================
# Test 407: NumericalIntegrator
# ============================================================
print("\n" + "=" * 60)
print("TEST 407: NumericalIntegrator")
print("=" * 60)

try:
    from chemmcp.tools import NumericalIntegrator

    tool = NumericalIntegrator()
    check("Instantiate", tool is not None)

    # Simpson's rule: ∫₀¹ x² dx = 1/3
    result = tool.run_code(
        mode="function",
        method="simpson",
        func_expr="x**2",
        lower_bound=0.0,
        upper_bound=1.0,
        n_points=100,
    )
    check("Simpson integral value", abs(result["integral_value"] - 1.0 / 3.0) < 1e-4,
          f"got {result['integral_value']}, expected 0.33333")

    # Trapezoidal: ∫₀^π sin(x) dx = 2
    result2 = tool.run_code(
        mode="function",
        method="trapezoidal",
        func_expr="sin(x)",
        lower_bound=0.0,
        upper_bound=math.pi,
        n_points=1000,
    )
    check("Trapezoidal sin integral", abs(result2["integral_value"] - 2.0) < 1e-2,
          f"got {result2['integral_value']}")

    # Gaussian quadrature: ∫₀¹ exp(-x²) dx ≈ 0.746824
    result3 = tool.run_code(
        mode="function",
        method="gaussian",
        func_expr="exp(-x**2)",
        lower_bound=0.0,
        upper_bound=1.0,
        n_points=20,
    )
    expected_erf = 0.746824
    check("Gaussian exp(-x²) integral", abs(result3["integral_value"] - expected_erf) < 1e-3,
          f"got {result3['integral_value']}")

    # Adaptive: ∫₀¹ x³ dx = 0.25
    result4 = tool.run_code(
        mode="function",
        method="adaptive",
        func_expr="x**3",
        lower_bound=0.0,
        upper_bound=1.0,
    )
    check("Adaptive x³ integral", abs(result4["integral_value"] - 0.25) < 1e-6,
          f"got {result4['integral_value']}")

    # Data mode
    x_data = [0.0, 0.5, 1.0]
    y_data = [0.0, 0.25, 1.0]  # x²
    result5 = tool.run_code(mode="data", method="trapezoidal", x_data=x_data, y_data=y_data)
    check("Data mode integral", abs(result5["integral_value"] - 0.375) < 0.01)

except Exception as e:
    failed += 1
    errors.append(("407-Import/Setup", str(e)))
    print(f"  ❌ ERROR: {e}")


# ============================================================
# Test 408: OdeSolverRk4
# ============================================================
print("\n" + "=" * 60)
print("TEST 408: OdeSolverRk4")
print("=" * 60)

try:
    from chemmcp.tools import OdeSolverRk4

    tool = OdeSolverRk4()
    check("Instantiate", tool is not None)

    # First-order decay: dy/dt = -ky, y(0)=1 → y(t)=exp(-kt)
    # At t=10, k=0.5: y = exp(-5) ≈ 0.006738
    result = tool.run_code(
        func_expr="-k * y",
        initial_values=[1.0],
        t_span=[0.0, 10.0],
        n_steps=200,
        params={"k": 0.5},
    )
    check("RK4 returns time points", len(result["time_points"]) > 1)
    check("RK4 returns solution", len(result["solution"]) == len(result["time_points"]))
    check("Final decay value", abs(result["final_value"][0] - math.exp(-5.0)) < 1e-3,
          f"got {result['final_value'][0]}, expected {math.exp(-5)}")
    check("Initial value preserved", abs(result["solution"][0][0] - 1.0) < 1e-6)

    # Consecutive reactions: A→B→C using y0,y1 variable names
    result2 = tool.run_code(
        func_expr=["-k1 * y0", "k1 * y0 - k2 * y1"],
        initial_values=[1.0, 0.0],
        t_span=[0.0, 5.0],
        n_steps=100,
        params={"k1": 1.0, "k2": 0.5},
    )
    check("System ODE vars", len(result2["final_value"]) == 2)
    check("System mass conservation-ish", result2["final_value"][0] >= 0 and result2["final_value"][1] >= 0)
    check("A decreases", result2["final_value"][0] < 1.0)

    # Simple harmonic oscillator: dy/dt=v, dv/dt=-y
    result3 = tool.run_code(
        func_expr=["y1", "-y0"],
        initial_values=[0.0, 1.0],
        t_span=[0.0, 2 * math.pi],
        n_steps=500,
        params={},
    )
    # After one period, should return to near initial state
    check("Oscillator period-y", abs(result3["final_value"][0]) < 0.02,
          f"y(2π)={result3['final_value'][0]}")
    check("Oscillator period-v", abs(result3["final_value"][1] - 1.0) < 0.03,
          f"v(2π)={result3['final_value'][1]}")

except Exception as e:
    failed += 1
    errors.append(("408-Import/Setup", str(e)))
    print(f"  ❌ ERROR: {e}")


# ============================================================
# Test 409: OdeSolverStiff
# ============================================================
print("\n" + "=" * 60)
print("TEST 409: OdeSolverStiff")
print("=" * 60)

try:
    from chemmcp.tools import OdeSolverStiff

    tool = OdeSolverStiff()
    check("Instantiate", tool is not None)

    # Stiff decay: dy/dt = -50*y, y(0)=1
    result = tool.run_code(
        func_expr="-k * y",
        initial_values=[1.0],
        t_span=[0.0, 1.0],
        n_steps=20,
        params={"k": 50.0},
    )
    check("Stiff solver returns results", isinstance(result, dict))
    check("Stiff final value near zero", result["final_value"][0] < 0.01,
          f"got {result['final_value'][0]}")
    check("Newton iterations recorded", result.get("newton_iterations_total", 0) > 0)

    # Classic stiff system: y1' = y2, y2' = -1000*y1 - 1001*y2
    result2 = tool.run_code(
        func_expr=["y1", "-1000 * y0 - 1001 * y1"],
        initial_values=[1.0, 0.0],
        t_span=[0.0, 1.0],
        n_steps=200,
        params={},
    )
    check("Stiff system solved", isinstance(result2["solution"], list))
    check("Stiff system final values", len(result2["final_value"]) == 2)
    # Both should be very close to zero at t=1
    check("Stiff y1 decays", abs(result2["final_value"][0]) < 1.0,
          f"y1={result2['final_value'][0]}")
    check("Stiff y2 decays", abs(result2["final_value"][1]) < 1.0,
          f"y2={result2['final_value'][1]}")

    # Compare stiff vs non-stiff for mild problem
    result_mild = tool.run_code(
        func_expr="-0.1 * y",
        initial_values=[1.0],
        t_span=[0.0, 10.0],
        n_steps=50,
        params={},
    )
    expected = math.exp(-1.0)  # exp(-0.1*10)
    check("Mild problem accurate", abs(result_mild["final_value"][0] - expected) < 0.05,
          f"got {result_mild['final_value'][0]}, expected {expected}")

except Exception as e:
    failed += 1
    errors.append(("409-Import/Setup", str(e)))
    print(f"  ❌ ERROR: {e}")


# ============================================================
# Test 410: PdeSolverFiniteDiff
# ============================================================
print("\n" + "=" * 60)
print("TEST 410: PdeSolverFiniteDiff")
print("=" * 60)

try:
    from chemmcp.tools import PdeSolverFiniteDiff

    tool = PdeSolverFiniteDiff()
    check("Instantiate", tool is not None)

    # 1D Heat equation with stable parameters
    result = tool.run_code(
        pde_type="heat_1d",
        alpha=0.01,
        length=1.0,
        t_final=0.5,
        nx=21,
        nt=1000,
        initial_condition="sin(pi * x)",
        bc_left_type="dirichlet",
        bc_left_value=0.0,
        bc_right_type="dirichlet",
        bc_right_value=0.0,
        scheme="explicit",
    )
    check("PDE returns grid", "spatial_grid" in result)
    check("PDE returns snapshots", "solution_snapshots" in result)
    check("PDE stability param", "stability_param" in result)
    check("PDE has snapshots", len(result["solution_snapshots"]) >= 2)
    check("BCs satisfied (left)", abs(result["final_profile"][0][0]) < 1e-6)
    check("BCs satisfied (right)", abs(result["final_profile"][0][-1]) < 1e-6)
    # For heat equation with sin IC, amplitude should decrease over time
    snap0 = result["solution_snapshots"][0]
    snap_last = result["solution_snapshots"][-1]
    first_max = max(snap0) if isinstance(snap0[0], (int, float)) else max(max(row) for row in snap0)
    last_max = max(snap_last) if isinstance(snap_last[0], (int, float)) else max(max(row) for row in snap_last)
    check("Amplitude decays", last_max < first_max,
          f"first_max={first_max:.4f}, last_max={last_max:.4f}")

    # Implicit scheme (always stable)
    result2 = tool.run_code(
        pde_type="heat_1d",
        alpha=0.1,
        length=1.0,
        t_final=0.1,
        nx=11,
        nt=50,
        initial_condition="x * (1 - x)",
        scheme="implicit",
    )
    check("Implicit scheme works", len(result2["solution_snapshots"]) >= 2)
    check("Implicit BC left", abs(result2["final_profile"][0][0]) < 1e-6)
    check("Implicit BC right", abs(result2["final_profile"][0][-1]) < 1e-6)

    # Wave equation
    result3 = tool.run_code(
        pde_type="wave_1d",
        alpha=1.0,
        length=1.0,
        t_final=1.0,
        nx=51,
        nt=500,
        initial_condition="sin(pi * x)",
        scheme="explicit",
    )
    check("Wave equation runs", len(result3["solution_snapshots"]) >= 2)

    # Neumann boundary conditions
    result4 = tool.run_code(
        pde_type="heat_1d",
        alpha=0.001,
        length=1.0,
        t_final=0.1,
        nx=11,
        nt=500,
        initial_condition="1.0",
        bc_left_type="neumann",
        bc_left_value=0.0,
        bc_right_type="neumann",
        bc_right_value=0.0,
        scheme="explicit",
    )
    check("Neumann BC works", len(result4["solution_snapshots"]) >= 2)

except Exception as e:
    failed += 1
    errors.append(("410-Import/Setup", str(e)))
    print(f"  ❌ ERROR: {e}")


# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
print(f"RESULTS: {passed} passed, {failed} failed out of {passed + failed} total")
print("=" * 60)

if errors:
    print("\nFailed tests:")
    for name, msg in errors:
        print(f"  • {name}: {msg}")

sys.exit(0 if failed == 0 else 1)
