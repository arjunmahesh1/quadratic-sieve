import math
from typing import List, Tuple

def trial_division(n: int, bound: int = 100) -> Tuple[int, List[int]]:
    """
    - remove all prime factors <= bound from n by trial division.
    - returns the remaining cofactor and a list of removed small primes.
    (this is good if the number is large and has small prime factors)
    """
    # TODO: just trial division
    pass


def is_probable_prime(n: int) -> bool:
    """
    - perform a strong probable-prime test on n.
    - returns True if n is (probably) prime, False otherwise
    (this is for optimization, faster than whole algorithm and if works then speed-up)
    """
    # TODO: implement a Miller-Rabin test
    pass


# factor-base generation

def legendre_symbol(a: int, p: int) -> int:
    """
    - this is to check if x^2 = a mod p has a solution 
    - compute the Legendre symbol (a|p) for odd prime p.
    """
    # TODO: implement Euler criterion from article or use built-in
    pass


def generate_factor_base(n: int, B: int) -> List[int]:
    """
    - build the factor base: all primes p <= B such that n is a quadratic residue mod p.
    - uses legendre_symbol to check if (n|p) = 1 and keeps if so.
    """
    factor_base = []
    # TODO: sieve primes up to B, check (n|p) = 1
    return factor_base


# Sieve part

def find_relations(n: int,
                   factor_base: List[int],
                   required: int) -> List[Tuple[int, List[int]]]:
    """
    - collect relations of the form: x^2 - n = product_{i}(p_i^{e_i}) over the factor_base.  
    - returns a list of tuples (x, exponent_vector).

    """
    relations = []
    # TODO: this is how i think it will go
    #   for x from ceil(sqrt(n)) up:
    #   compute y = x^2 - n
    #   for each p in factor_base, divide y by p repeatedly
    #   if remainder is +-1, record exponent vector
    #   stop when len(relations) >= required
    return relations

# matrix and linear algebra stuff

def build_matrix(relations: List[Tuple[int, List[int]]]) -> List[List[int]]:
    """
    - build a binary matrix (mod 2) from the exponent vectors in relations 
    - rows correspond to relations, columns to factor_base primes
    """
    # TODO: extract exponent vectors, reduce mod 2
    pass


def gaussian_elimination_mod2(matrix: List[List[int]]) -> List[int]:
    """
    - perform Gaussian elimination over GF(2).
    - return a nontrivial null space vector indicating which rows (relations)
    combine to give an even exponent combination.
    - MAKE SURE this is done mod 2
    """
    # TODO: implement bit-packed elimination or simple list-of-lists version
    pass



# square-Root & GCD Extraction

def recover_squares(n: int,
                    relations: List[Tuple[int, List[int]]],
                    combination: List[int]) -> Tuple[int, int]:
    """
    - this multiplies out the selected relations to compute:
      A = product x_i mod n,
      B = sqrt(prod y_i) mod n,
    where y_i = x_i^2 - n.
    - returns (A, B).
    """
    # TODO: multiply selected x_i, reconstruct B from exponent sums
    pass


def extract_factor(n: int, A: int, B: int) -> Tuple[int, int]:
    """
    - compute gcd(A - B, n) and gcd(A + B, n).
     -return any nontrivial factor pair (g, n//g).
    """
    # TODO: probably use math.gcd if we are allowed
    pass


# putting together

def quadratic_sieve(n: int,
                     B: int,
                      extra: int = 5) -> Tuple[int, int]:
    """
    - main function:
      n    is integer to factor
      B:     factor‐base bound
      extra: this is slack for number of relations collected (saw online this is how to do it)

    returns a nontrivial factor pair (p, q) with p*q = n.
    """
    # 1) trial division and primality check
    cofactor, small_primes = trial_division(n)
    if is_probable_prime(cofactor):
        return (cofactor, 1)

    # 2)  factor base
    fb = generate_factor_base(cofactor, B)

    # 3) collect relations
    num_relations = len(fb) + extra
    relations = find_relations(cofactor, fb, num_relations)

    # 4) build matrix and find null‐space
    M = build_matrix(relations)
    combination = gaussian_elimination_mod2(M)

    # 5) get squares and extract factors
    A, Bval = recover_squares(cofactor, relations, combination)
    return extract_factor(cofactor, A, Bval)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python3 quadratic_sieve.py <n> [B]")
        sys.exit(1)

    n = int(sys.argv[1])
    B = int(sys.argv[2]) if len(sys.argv) >= 3 else 1000

    p, q = quadratic_sieve(n, B)
    print(f"factors of {n}: {p} * {q}")
    