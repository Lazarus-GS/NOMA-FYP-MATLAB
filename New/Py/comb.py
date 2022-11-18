import math


# Counts the number of grids where the first i rows and first j columns are
# monochrome. All of the monochrome rows and columns must have the same color if
# and only if both i and j are nonzero. Otherwise, the color choices are
# independent.
def count_monochrome(n, i, m, j):
    return 3 ** (1 if i and j else i + j) * 3 ** ((n - i) * (m - j))


# Uses inclusion-exclusion to count the number of grids with no monochrome row
# or column.
def count(n, m=3):
    # There are math.comb(n, i) * math.comb(m, j) intersections with i
    # monochrome rows and j monochrome columns.
    return sum(
        math.comb(n, i)
        * math.comb(m, j)
        * (-1) ** (i + j)
        * count_monochrome(n, i, m, j)
        for i in range(n + 1)
        for j in range(m + 1)
    ) % (10 ** 9 + 7)


print(count(4))