import math


def PDF(args, c, GT=True):
    """
    Probability density function (PDF) for a Gaussian distribution.

    Parameters:
    - args: a tuple containing μ (mean) and σ (standard deviation)
    - c: the upper limit of integration (or lower if GT is False)
    - GT: boolean indicating if we want the probability of x being greater than c (True) or less than c (False)

    Returns:
    The probability density of x at the point c.
    """
    mu, sigma = args

    # Gaussian PDF function
    def gaussian_pdf(x):
        return (1 / (sigma * math.sqrt(2 * math.pi))) * math.exp(-0.5 * ((x - mu) / sigma) ** 2)

    return gaussian_pdf(c)


def simpsons_rule(pdf, mu, sigma, c, GT, n=1000):
    """
    Simpson's 1/3 rule for numerical integration to find the probability of x <= c or x > c.

    Parameters:
    - pdf: probability density function to integrate
    - mu: mean of the distribution
    - sigma: standard deviation of the distribution
    - c: integration limit
    - GT: if True, integrate from c to positive infinity, else from negative infinity to c
    - n: number of intervals (must be even)

    Returns:
    The probability calculated using numerical integration.
    """
    # Defining integration limits
    a = mu - 5 * sigma if not GT else c
    b = c if not GT else mu + 5 * sigma

    # Width of the intervals
    h = (b - a) / n

    # Calculating Simpson's Rule
    integration = pdf((mu, sigma), a, GT=False) + pdf((mu, sigma), b, GT=False)
    for i in range(1, n):
        x = a + i * h
        if i % 2 == 0:  # Even index terms
            integration += 2 * pdf((mu, sigma), x, GT=False)
        else:  # Odd index terms
            integration += 4 * pdf((mu, sigma), x, GT=False)
    integration *= h / 3

    # If we're looking for the probability of x > c, we subtract from 1
    if GT:
        return 1 - integration
    return integration

# Testing the functions with the given values
# Note: Since the Simpson's rule requires an even number of intervals, n must be even.
# For demonstration, we'll use n = 1000 for a fine approximation.

# We will not print results here, as we'll define the main function next to handle the printing.
def main():
    """
    Main function to calculate the probability of x using the PDF and Simpson's rule.
    It prints the results for specified values of mu, sigma, and c.
    """
    # Test cases as provided in the task
    test_cases = [
        (100, 12.5, 105),  # (mu, sigma, c)
        (100, 3, 106)  # (mu, sigma, c + 2*sigma)
    ]

    # Iterate over the test cases and calculate probabilities
    for mu, sigma, c in test_cases:
        probability = simpsons_rule(PDF, mu, sigma, c, GT=False)
        print(f'P(x<{c} | N({mu},{sigma}^2)) = {probability:.2f}')


# Call main function to print the results
main()
