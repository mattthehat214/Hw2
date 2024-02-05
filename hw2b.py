import math


def Secant(fcn, x0, x1, maxiter=10, xtol=1e-5):
    """
    Secant method for finding the root of a function.

    Parameters:
    fcn: The function for which we want to find the root.
    x0, x1: Two initial guesses in the neighborhood of the root.
    maxiter: Maximum number of iterations.
    xtol: Tolerance for convergence.

    Returns:
    The estimated root of the function.
    """
    for _ in range(maxiter):
        fx0 = fcn(x0)
        fx1 = fcn(x1)
        if fx1 - fx0 == 0:  # Prevent division by zero
            return None
        # Secant method formula to find a new approximation
        x_new = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        if abs(x_new - x1) < xtol:  # Check for convergence
            return x_new
        x0, x1 = x1, x_new  # Update guesses
    return x1  # Return the last estimate if no convergence


# Now we define the main function that will use the Secant function with the given parameters
def main():
    """
    Main function that uses the Secant method to find roots of given functions.
    """

    # Define the functions for which we need to find roots
    def f1(x): return x - 3 * math.cos(x)

    def f2(x): return math.cos(2 * x) - x

    def f3(x): return math.cos(2 * x) - x ** 3

    # Define the parameters for the Secant method calls
    test_cases = [
        (f1, 0.1, 1.2, 5, 1e-4),
        (f2, 0.1, 1.2, 15, 1e-8),
        (f3, 0.1, 1.2, 3, 1e-8)
    ]

    # Call the Secant method for each case and print the results
    for f, x0, x1, maxiter, xtol in test_cases:
        root = Secant(f, x0, x1, maxiter, xtol)
        print(f"The root found for function {f.__name__} is: {root:.10f}")


# Call the main function to execute the Secant method on the test cases
main()
