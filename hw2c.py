def gauss_seidel(aug, x, n_iter=15):
    """
    Performs the Gauss-Seidel iteration for solving the linear equations Ax = b.

    :param aug: An augmented matrix [A|b] for the system of equations.
    :param x: Initial guess for the solution vector.
    :param n_iter: Number of iterations to perform.
    :return: The estimated solution vector.
    """
    n = len(aug)
    for _ in range(n_iter):
        x_new = x[:]
        for i in range(n):
            sum_ax = sum(aug[i][j] * x_new[j] for j in range(n) if j != i)
            x_new[i] = (aug[i][-1] - sum_ax) / aug[i][i]
        x = x_new[:]
    return x


def main():
    # System 1
    aug1 = [
        [3, 1, -1, 2],
        [1, 4, 1, 12],
        [2, 1, 2, 10]
    ]
    initial_guess1 = [0, 0, 0]
    result1 = gauss_seidel(aug1, initial_guess1)
    print("Solution for system 1:", result1)

    # System 2
    aug2 = [
        [1, -10, 2, 4, 2],
        [3, 1, 4, 12, 12],
        [9, 2, 3, 4, 21],
        [-1, 2, 7, 3, 37]
    ]
    initial_guess2 = [0, 0, 0, 0]
    result2 = gauss_seidel(aug2, initial_guess2)
    print("Solution for system 2:", result2)


if __name__ == "__main__":
    main()
