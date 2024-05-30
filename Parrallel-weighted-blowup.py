import sympy
import numpy as np
from fractions import Fraction
from sortedcollections import ValueSortedDict
import tqdm

def ideal_to_monomial_terms(ideal:sympy.polys.agca.ideals.Ideal):
    monomial_terms = np.array([term[0] for gen in ideal.gens for term in gen.terms()])
    return monomial_terms

def prepros(n_vars):
    out = ValueSortedDict(sum)
    for i in range(n_vars):
        endpoint = [0]*n_vars
        endpoint[i] = np.inf
        out[i] = tuple(endpoint)
    return out

def project_used(point, outpoints:ValueSortedDict):
    for index, w_var in outpoints.items():
        if len(np.nonzero(point)[0]) == 1:
            return point
        if point[index] == 0 or w_var[index] >= np.inf:
            continue
        if abs(w_var[index] - point[index])<0.1e-10 or w_var[index] < point[index]:
            scale = np.inf
        else:
            scale = w_var[index] / (w_var[index] - point[index])
        #print(f"Scale was {scale}, with point {point} and {w_var}")
        point[index] = 0
        with np.errstate(invalid='ignore'):
            point = point * scale
        point[np.isnan(point)] = 0
    return point


def handlepoints(points, verbose = False):
    n = points.shape[1]
    outpoints = prepros(n)
    for point in points:
        if verbose > 2:print("Point", point, "outpoints", outpoints)
        point = project_used(point, outpoints)
        index = np.nonzero(point)[0]
        weight = np.sum(point)
        for i in index:
            if weight >= outpoints[i][i]:
                continue
            newpoint = [0]*n
            newpoint[i] = weight
            outpoints[i] = tuple(newpoint)
    out = [outpoints[i][i] for i in range(n)]
    return out

def deg_sort(points):
    index = np.lexsort(
        (np.sum(points!=0, axis=1),
         np.sum(points, axis=1)))
    points = points[index]
    return points

def inv_finder(ideal:sympy.polys.agca.ideals.Ideal, verbose = False):
    vars = ideal.ring.symbols
    n = len(vars)
    points = ideal_to_monomial_terms(ideal) # returns (flat) list like [(a,b,c),...]
    points = deg_sort(points)
    if verbose >= 2:print(points)
    points1 = handlepoints(points, verbose)
    if verbose >= 2:print(points1)
    inv = [Fraction(i).limit_denominator(1000) if i < np.inf else np.inf for i in points1]
    return inv

def make_ideal(*inputs, gen_out = False):
    vars = set()
    generators = []
    for s in inputs:
        f = sympy.sympify(s)
        vars = vars.union(f.free_symbols)
        generators.append(f)
    vars = sorted(list(vars), key=lambda x: str(x))
    ring = sympy.QQ.old_poly_ring(*vars)
    ideal = ring.ideal(*generators)
    if gen_out:
        return ideal, generators
    return ideal

def test_admisibility(ideal:sympy.polys.agca.ideals.Ideal, inv):
    gens = list(ideal.gens)
    for poly in gens:
        for term in poly.terms():
            if sum(term[0][i]/inv[i] for i in range(len(inv))) < 1-0.1e-10:
                print(f"AssertionError {sympy.pprint(poly), term} did not work with {inv}.")
                return False
    return True

def make_weights(inv):
    denoms_lcm = sympy.lcm([frac.denominator for frac in inv if frac != np.inf])
    temp_inv = [(frac*denoms_lcm) if frac != np.inf else np.inf for frac in inv]
    nume_lcm = sympy.lcm([frac.numerator for frac in temp_inv if frac != np.inf])
    out = [round(1/(frac/nume_lcm)) if frac != np.inf else 0 for frac in temp_inv]
    return tuple(out)

def blow_up(ideal, generators, weights, chart_index = 0, verbose = False, keep_vars = False, new_var = "u", fast = False):
    vars = ideal.ring.symbols
    n = len(vars)
    # Charts by default x0 chart with u corresponding to x0
    if keep_vars:
        #new_vars = vars[:chart_index] + (sympy.sympify(new_var),) + vars[chart_index+1:]
        new_vars = vars[:]
    else:
        new_vars = sympy.symbols(' '.join([f"y{i}" if i != chart_index else new_var for i in range(n)]))
    u = new_vars[chart_index]
    u_var = {vars[chart_index]:u**weights[chart_index]}
    loc_vars = {vars[i]:u**weight * new_vars[i] if i != chart_index 
                else u for i, weight in enumerate(weights)}
    if verbose >= 1:
        print("====", "chart ", vars[chart_index], "====")
        print(loc_vars)

    
    new_variety = [gen.subs(u_var).subs(loc_vars) for gen in generators]
    all_powers = [term.as_coeff_exponent(u)[1] if term.has(u) 
                  else 0 
                  for eq in new_variety
                  for term in eq.as_ordered_terms()] # Finds the number of u's one can pull out

    # Calculate the minimum of the maximum powers
    max_power = min(all_powers)

    # Divide each term by u raised to the max_power
    simplified_expressions = [poly / (u**max_power) 
                           for poly in new_variety]
    if verbose >= 1:print(simplified_expressions)
    if fast:return simplified_expressions
    # Simplify the divided expression
    simplified_expressions = [sympy.simplify(divided_expression) 
                              for divided_expression in simplified_expressions]
    results = [(u**max_power) * simplified_expression 
               for simplified_expression in simplified_expressions]
    # Multiply the simplified expression by u raised to the max_power
    if verbose >= 1:
        print(results)
    return results, simplified_expressions
    
    
def multiblowup(ideal:str, ring = None, charts = [], verbose = False, secure = True, path = None, stop_depth = 0):
    if ring == None:
        ideal, generators = make_ideal(*ideal.split(","), gen_out=True)
        ring = ideal.ring
    else:
        generators = ideal
        ideal = ring.ideal(*ideal)
    #if generators[0] == 1: THIS IS FOR THE SPECIFIC APPLICATION!!
    #    return charts + ["success"]
    if ideal.contains(1):
        if verbose > 0.5:print("with", ring.symbols, "though ", charts, "0 no longer singular")
        return charts
    inv = inv_finder(ideal, verbose = verbose)
    if all(i == 1 for i in inv):
        if verbose > 0.5:print("with", ring.symbols, "though ", charts, "0 no longer singular")
        return charts
    if verbose > 0.5:print("invariant =", inv, f" at chart {charts}"*(len(charts)>1))
    if secure: 
        if not test_admisibility(ideal, inv):
            print("with", ring.symbols, "though ", charts, "ERROR")
            if secure>1:
                assert test_admisibility(ideal, inv)
    weights = make_weights(inv)
    if verbose > 0.5:print("Weights^{-1}", weights, ring.symbols)
    tree = []
    for index, weight in enumerate(weights):
        if path and len(charts) < len(path):
            path_index = path[len(charts)]
            if index != path_index:
                continue
            if weight == 0:
                return
        if stop_depth and len(charts) >= stop_depth:
            return charts + ["end"]
        if weight == 0:
            continue
        out = blow_up(ideal, generators, weights, chart_index=index, verbose = verbose, keep_vars=True, fast=True)
        chart_var = f"{ideal.ring.symbols[index]}**{weight}"
        leaf = multiblowup(out, ring, charts = charts + [chart_var], verbose = verbose, secure = secure, path=path, stop_depth=stop_depth)
        tree.append(leaf)
    return tree
    
    
    
import multiprocessing
def multiblowup_helper(args):
    return multiblowup(*args)
def multiblowup_parallel(ideal, verbose = 0, stop_depth = 0):
    ideal, generators = make_ideal(*ideal.split(","), gen_out=True)
    ring = ideal.ring
    n_cores = multiprocessing.cpu_count()//2
    n_paths = len(ring.symbols)
    #args = [(generators, ring, [], verbose, True, [8,i], stop_depth) for i in range(n_paths)]  # Only t 
    args = [(generators, ring, [], verbose, True, [i], stop_depth) for i in range(n_paths)]
    pool = multiprocessing.Pool(processes=min(n_paths, n_cores))
    out = list(tqdm.tqdm(pool.imap_unordered(multiblowup_helper, args), total=len(args)))
    return out

#==========For calculating thom k=4. See https://arxiv.org/pdf/2012.06425.pdf for more details========
import numpy as np
w0 = np.zeros(4)
w1 = np.array([1,0,0,0])
w2 = np.array([*sympy.symbols("b12,b22"), 0, 0])
w3 = np.array([*sympy.symbols("b13,b23,b33"), 0])
w4 = np.array([*sympy.symbols("b14,b24,b34,b44")])
tt = np.array([sympy.symbols("t")])
T = np.hstack((w1, np.outer(w0, w0).flatten() , np.outer(w0, np.outer(w0,w0)).flatten()))
U = np.hstack((w2, np.outer(tt*w1, w1).flatten() , np.outer(w0, np.outer(w0,w0)).flatten()))
V = np.hstack((w3, np.outer(tt*w1, w2).flatten() , np.outer(tt*tt*w1, np.outer(w1,w1)).flatten()))
W = np.hstack((w4, np.outer(tt*w1, w3).flatten() + np.outer(tt*w2, w2).flatten(), np.outer(tt*tt*w1, np.outer(w1,w2)).flatten()))
M = np.vstack((T,U,V,W))
MM = np.hstack((M,tt*tt*tt*np.array([0,0,0,1]).reshape(4,1)))
ring = [x for x in np.hstack((w0,w1,w2,w3,w4,tt)) if type(x)==sympy.Symbol]
from itertools import combinations
def minors2d(minor_size:int, matrix:np.ndarray, reduce = True, nonzero = True, asideal = True):
    if reduce: #reduce complexity by removing 0 cols
        matrix = matrix[:, ~np.all(matrix == 0, axis = 0)]
    nrows, ncols = matrix.shape
    assert nrows >= minor_size and ncols >= minor_size
    out = []
    row_choices = combinations(range(nrows), minor_size)
    col_choices = combinations(range(ncols), minor_size)
    for minor_row in row_choices:
        for minor_col in col_choices:
            minor = sympy.det(sympy.Matrix(matrix[minor_row,:][:,minor_col]))
            if minor != 0:
                out.append(minor)
    return out
def make_monomial(ideal:list):
    out = ",".join(str(s) for gen in ideal for s in sympy.Add.make_args(sympy.expand(gen)))
    return out
#============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Completely blow up some ideal.')
    parser.add_argument('filename', type=str,
                        help='file_name to output to'
                    )
    parser.add_argument('--ideal', 
                        metavar='I', 
                        type=str, 
                        default=make_monomial(minors2d(4,MM)),
                        help='Remember ""... The ideal to blow up. Like "x**2,y**3" (default: k=4 thom)',
                        )
    parser.add_argument('--stop_depth', type=int,
                    help='add if you want to stop at max depth (int)',
                    default=0)
    parser.add_argument('--verbose', type=int,
                    help='verbose \in 0,..,3',
                    default=0)
    

    args = parser.parse_args()
    print(args.filename, args.ideal)

    out = multiblowup_parallel(args.ideal, args.verbose, args.stop_depth) # CHANGE THIS
    # out = multiblowup_parallel(make_monomial(minors2d(4,MM)))
    with open(f"{args.filename}.txt", "w") as text_file:
        text_file.write(str(out))
    print(out)
