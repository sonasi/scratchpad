from fractions import Fraction as Fr
 
def bernoulli(n):
    A = [0] * (n+1)
    for m in range(n+1):
        A[m] = Fr(1, m+1)
        for j in range(m, 0, -1):
          A[j-1] = j*(A[j-1] - A[j])
	A[0] # (which is Bn)
    return float(A[0].numerator)/float(A[0].denominator)
 
