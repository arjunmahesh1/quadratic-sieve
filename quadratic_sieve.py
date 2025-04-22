import math
import random
from typing import List, Tuple

def primes_up_to(n: int) -> List[int]:
    sieve = [True] * (n + 1)
    sieve[0:2] = [False, False]
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            for multiple in range(p*p, n+1, p):
                sieve[multiple] = False
    return [p for p, is_prime in enumerate(sieve) if is_prime]

def estimate_optimal_B(n: int) -> int:
    """
    Heuristic to estimate the optimal factor base bound B
    B = exp(0.5 * sqrt(ln(n) * ln(ln(n))))
    (run with just n to auto-generate B)
    """
    ln_n = math.log(n)

    return max(10, int(math.exp(0.5 * math.sqrt(ln_n * math.log(ln_n)))))

def trial_division(n: int, bound: int = 100) -> Tuple[int, List[int]]:
    """
    - remove all prime factors <= bound from n by trial division.
    - returns the remaining cofactor and a list of removed small primes.
    (this is good if the number is large and has small prime factors)
    """
    # trial division
    cofactor = n
    small_factors: List[int] = []
    for p in primes_up_to(bound):
        while cofactor % p == 0:
            cofactor //= p
            small_factors.append(p)

    return cofactor, small_factors


def is_probable_prime(n: int, k: int=5) -> bool:
    """
    - perform a strong probable-prime test on n.
    - returns True if n is (probably) prime, False otherwise
    (this is for optimization, faster than whole algorithm and if works then speed-up)
    """
    # Miller-Rabin test
    if n < 2:
        return False
    
    # small primes trial
    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23]
    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False
        
    # write n - 1 as 2^s * d
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # witness loop
    for _ in range(k):
        a = random.randrange(2, n - 1)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for __ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
        
    return True


# factor-base generation

def legendre_symbol(a: int, p: int) -> int:
    """
    - this is to check if x^2 = a mod p has a solution 
    - compute the Legendre symbol (a|p) for odd prime p.
    """
    # Euler criterion from article or use built-in
    ls = pow(a % p, (p - 1) // 2, p)
    if ls == p - 1:
        return -1
    
    return ls


def generate_factor_base(n: int, B: int) -> List[int]:
    """
    - build the factor base: all primes p <= B such that n is a quadratic residue mod p.
    - uses legendre_symbol to check if (n|p) = 1 and keeps if so.
    """
    # sieve primes up to B, check (n|p) = 1

    fb: List[int] = []
    for p in primes_up_to(B):
        if p == 2:
            # include 2 if n is odd
            if n % 2 == 1:
                fb.append(2)
        else:
            if legendre_symbol(n, p) == 1:
                fb.append(p)

    return fb


# Sieve part

def find_relations(n: int,
                   factor_base: List[int],
                   required: int) -> List[Tuple[int, List[int]]]:
    """
    - collect relations of the form: x^2 - n = product_{i}(p_i^{e_i}) over the factor_base.  
    - returns a list of tuples (x, exponent_vector).

    """
    relations: List[Tuple[int, List[int]]] = []
    # this is how i think it will go
    #   for x from ceil(sqrt(n)) up:
    #   compute y = x^2 - n
    #   for each p in factor_base, divide y by p repeatedly
    #   if remainder is +-1, record exponent vector
    #   stop when len(relations) >= required

    x = math.isqrt(n) + 1
    while len(relations) < required:
        y = x*x - n
        val = abs(y)
        exponents = [0] * len(factor_base)
        for i, p in enumerate(factor_base):
            while val % p == 0:
                val //= p
                exponents[i] += 1
        if val == 1:
            relations.append((x, exponents))
        x += 1

    return relations

# matrix and linear algebra stuff

def build_matrix(relations: List[Tuple[int, List[int]]]) -> List[List[int]]:
    """
    - build a binary matrix (mod 2) from the exponent vectors in relations 
    - rows correspond to relations, columns to factor_base primes
    """
    # extract exponent vectors, reduce mod 2
    
    return [[e % 2 for e in exponents] for _, exponents in relations]


def gaussian_elimination_mod2(matrix: List[List[int]]) -> List[int]:
    """
    - perform Gaussian elimination over GF(2).
    - return a nontrivial null space vector indicating which rows (relations)
    combine to give an even exponent combination.
    - MAKE SURE this is done mod 2
    """
    # implement bit-packed elimination or simple list-of-lists version
    rows = len(matrix)
    cols = len(matrix[0]) if rows > 0 else 0

    # Transpose to equations A * sol = 0, where A is (cols x rows)
    A = [[matrix[r][c] for r in range(rows)] for c in range(cols)]
    pivot_rows = {}
    r = 0

    # forward eliminate to RREF
    for c in range(rows):
    #for c in range(cols):
        if r >= cols:
            break
        sel = None
        for i in range(r, cols):
            if A[i][c] == 1:
                sel = i
                break
        if sel is None:
            continue
        A[r], A[sel] = A[sel], A[r]
        pivot_rows[r] = c

        # eliminate in all other rows
        for i in range(cols):
        #for i in range(rows):
            if i != r and A[i][c] == 1:
                for j in range(rows):
                    A[i][j] ^= A[r][j]
        r += 1

    # find free variable
    pivot_cols = set(pivot_rows.values())
    free = None
    for var in range(rows):
        if var not in pivot_cols:
            free = var
            break
    if free is None:
        raise ValueError("No nontrivial nullspace found")
    
    # build solution
    sol = [0] * rows
    sol[free] = 1

    # back-substitute for pivots
    for pr, pc in reversed(pivot_rows.items()):
        s = 0
        for j in range(rows):
            if j != pc and A[pr][j] == 1 and sol[j] == 1:
                s ^= 1
        sol[pc] = s

    return sol



# square-Root & GCD Extraction

def recover_squares(n: int,
                    fb: List[int],
                    relations: List[Tuple[int, List[int]]],
                    combination: List[int]) -> Tuple[int, int]:
    """
    - this multiplies out the selected relations to compute:
      A = product x_i mod n,
      B = sqrt(prod y_i) mod n,
    where y_i = x_i^2 - n.
    - returns (A, B).
    """
    # multiply selected x_i, reconstruct B from exponent sums

    fb_size = len(relations[0][1])

    # A = prod x_i
    A = 1
    for include, (x, _) in zip(combination, relations):
        if include:
            A = (A * x) % n

    # aggregate exponents
    exp_sum = [0] * fb_size
    for include, (_, exponents) in zip(combination, relations):
        if include:
            for i in range(fb_size):
                exp_sum[i] += exponents[i]

    # B = prod p_i^(exp_sum[i]//2)
    B = 1
    for p, e in zip(fb, exp_sum):
        B = (B * pow(p, e // 2, n)) % n

    return A, B


def extract_factor(n: int, A: int, B: int) -> Tuple[int, int]:
    """
    - compute gcd(A - B, n) and gcd(A + B, n).
     -return any nontrivial factor pair (g, n//g).
    """
    
    g = math.gcd(A - B, n)
    if 1 < g < n:
        return g, n // g
    g = math.gcd(A + B, n)
    if 1 < g < n:
        return g, n // g
    
    raise ValueError("Failed to extract nontrivial factor")


# putting together

def quadratic_sieve(n: int,
                     B: int=1000,
                      extra: int = 5) -> Tuple[int, int]:
    """
    - main function:
      n    is integer to factor
      B:     factor‐base bound
      extra: this is slack for number of relations collected (saw online this is how to do it)

    returns a nontrivial factor pair (p, q) with p*q = n.
    """
    # 1) trial division and primality check
    # fixed commented out part with code below
    '''
    cofactor, smalls = trial_division(n, B)
    if cofactor == 1:
        return n, 1
    if is_probable_prime(cofactor):
        return cofactor, n // cofactor
    '''

    # 1) pull off all primes <= B
    cofactor, smalls = trial_division(n, B)

    # 1a) if trial division completely factors n, peel off one small prime and return it
    if cofactor == 1:
        p = smalls[0]
        return p, n // p

    # 1b) if whatever remains is prime, return (product of smalls, cofactor)
    if is_probable_prime(cofactor, k=20):
        small_prod = 1
        for p in smalls:
            small_prod *= p
        return small_prod, cofactor


    # 2) build factor base, and if it’s too tiny keep doubling B
    fb = generate_factor_base(cofactor, B)
    MIN_FB_SIZE = 5
    while len(fb) < MIN_FB_SIZE:
        B *= 2
        print(f"Factor base too small (size={len(fb)}); doubling B to {B}")
        fb = generate_factor_base(cofactor, B)

    # 3) collect relations
    num_relations = len(fb) + extra
    relations = find_relations(cofactor, fb, num_relations)

    # 4) build matrix and find null‐space
    M = build_matrix(relations)
    combination = gaussian_elimination_mod2(M)

    # 5) get squares and extract factors
    A, Bval = recover_squares(cofactor, fb, relations, combination)

    return extract_factor(cofactor, A, Bval)


if __name__ == "__main__":
    import sys
    import time
    if len(sys.argv) < 2:
        print("Usage: python3 quadratic_sieve.py <n> [B]")
        sys.exit(1)

    n = int(sys.argv[1])
    if len(sys.argv) >= 3:
        B = int(sys.argv[2])
    else:
        B = estimate_optimal_B(n)
        # don’t let B be too small
        MIN_B = 100
        if B < MIN_B:
            print(f"Auto‐selected B = {B} is too small; bumping up to MIN_B = {MIN_B}")
            B = MIN_B
        else:
            print(f"Auto‐selected factor‐base bound B = {B}")

    # 1) pull off all small primes <= B
    cofactor, smalls = trial_division(n, B)
    small_prod = 1
    for p in smalls:
        small_prod *= p

    # 2) If trial division completely factors n, split into two factors
    if cofactor == 1:
        if not smalls:
            print(f"Factors of {n}: 1 * 1")
        else:
            # Split into first prime and the product of the rest
            p = smalls[0]
            q = n // p
            print(f"Factors of {n}: {p} * {q}")
        sys.exit(0)

    # 3) if what remains is prime, combine
    if is_probable_prime(cofactor, k=20):
        print(f"Factors of {n}: {small_prod} * {cofactor}")
        sys.exit(0)

    # 4) otherwise sieve to split the cofactor
    start = time.perf_counter()
    try:
        p_co, q_co = quadratic_sieve(cofactor, B)
        elapsed = time.perf_counter() - start

        f1 = small_prod * p_co
        f2 = q_co

        if f1 > f2:
            f1, f2 = f2, f1
        print(f"Factors of {n}: {f1} * {f2}")
        print(f"Elapsed time: {elapsed:.3f} seconds")
    except Exception as e:
        elapsed = time.perf_counter() - start
        print(f"Failed after {elapsed:.3f} seconds: {e}")

    
    