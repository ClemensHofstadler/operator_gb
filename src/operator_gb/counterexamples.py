from sage.sat.converters.polybori import CNFEncoder
from sage.sat.solvers.dimacs import DIMACS
from collections import Counter
from .quiver import Quiver
from itertools import chain, combinations
from sage.all import *

import sys
import time
import subprocess
import tempfile
import random
import copy
from pathlib import Path

#### LEAVE AS IS ########
_verbosity_ = 0
_TSEITIN_NAME_ = 'tseit__'
_RABINOWITSCH_NAME_= 'rab__'
_LIFTING_NAME_ = 'lift__'
_REPLACE_VARS_ = 'new__'


def most_common_pair(polies):
    mons = [m for poly in polies for c,m in poly]
    all_pairs = [w[i:i+2] for w in mons for i in range(len(w)-1)]
    value, count = Counter(all_pairs).most_common(1)[0]
    if count > 1:
        return value
    else:
        return None
        
##############################################################################        

# replace v by x
def replace_in_mon(m,v,x):
    i = m.find(v)
    if i < 0: return m
    return Word(m[:i]) + Word(x) + Word(m[i+len(v):])  

##############################################################################

def replace_word(polies,v,i):
    x = Word([_REPLACE_VARS_ + str(i)])  
    for i in range(len(polies)):
        polies[i] = [(c,replace_in_mon(m,v,x)) for c,m in polies[i]]           
    return polies, x

##############################################################################

def longest_mon(polies):
    mons = [m for poly in polies for c,m in poly]
    m = max(mons, key=lambda w: len(w))
    if len(m) > 2:
        return m
    else:
        return None
      
##############################################################################
    
def simplify_input(assumptions,f):    
    X = [str(v) for v in f.parent().gens()]    
    polies = [[(c, m.to_word()) for m,c in g.monomial_coefficients().items()] for g in [f] + assumptions]
    
    i = 0
    shortcuts = []
    
    # replace duplicate terms by new variables
    while True:
        v = most_common_pair(polies)
        if v is None:
            break
        polies, xi = replace_word(polies,v,i)
        shortcuts.append( [(1,xi), (-1,v)] )
        X.append(xi)
        i += 1 
        
    # replace all terms of length > 2 by new variables
    while True:
        v = longest_mon(polies)
        if v is None:
            break
        #  len = 3 -> only need 1 new var
        if len(v) == 3:
            polies, xi = replace_word(polies,v[:2],i)
            shortcuts.append( [(1,xi), (-1,v[:2])] )
            X.append(xi)
            i += 1 
        # len > 3 -> need 2 vars
        else:
            v1 = v[:len(v)//2]
            v2 = v[len(v)//2:]
            polies, x1 = replace_word(polies,v1,i)
            polies, x2 = replace_word(polies,v2,i+1)
            shortcuts.append( [(1,x1), (-1,v1)] )
            shortcuts.append( [(1,x2), (-1,v2)] )
            X += [x1,x2]
            i += 2 
    
    polies = ['+'.join(str(c) + '*' + '*'.join(list(m) + ['1']) for c,m in g) for g in polies]
    shortcuts = ['+'.join(str(c) + '*' + '*'.join(list(m)) for c,m in g) for g in shortcuts]
    
    F = FreeAlgebra(QQ,X)
    assumptions = [F(g) for g in polies[1:]]
    shortcuts = [F(g) for g in shortcuts]
    f = F(polies[0])
    
    return assumptions, shortcuts, f
    
##############################################################################    

def preprocess_input(assumptions, f, Q, dims):

    assumptions, shortcuts, f = simplify_input(assumptions, f)
        
    # update quiver
    for i,s in enumerate(shortcuts):
        lm = s.leading_monomial()
        U,V = Q.signature(lm).pop()
        Q.add_edge(U,V,_REPLACE_VARS_ + str(i)) 
        
    X = [str(v) for v in f.parent().gens()]
      
    matrix_sizes = {v : (dims[Q.target(v)[0]], dims[Q.source(v)[0]]) for v in X}
    U,V = Q.signature(f).pop()
    matrix_sizes.update({f : (dims[U], dims[V])})
       
    F = FreeAlgebra(ZZ,X)
    assumptions = [F(g) for g in  shortcuts + assumptions]
    f = F(f)    
    
    return assumptions, f, Q, matrix_sizes 

##############################################################################

#takes a list of polynomials obtained by the matrix, converts it to boolean and applies sat solver
def system_to_CNF(polynomials):
    r"""
    Solves a system of Polynomial equations over Z2.

    Input:
        - ``polynomials`` -- list of polynomials
    """
    if _verbosity_ > 0:
        print('Conversion to CNF')

    s = time.process_time() 
    
    if 1 in polynomials:
        raise ValueError("system is not solvable over Z2")
    polynomials = [f for f in polynomials if f != 0]
        
    # Encode to CNF
    R = polynomials[0].parent()
    solver = DIMACS()
    encoder = CNFEncoder(solver, R)
    for poly in polynomials:
        encoder.clauses_dense(poly)
            
    if _verbosity_ > 0:
        print(f'Conversion to CNF in {time.process_time() - s:.4f}')
    
    return solver, encoder

##############################################################################

def generate_commutative_system(assumptions, f, Q, matrix_sizes, symbolic=False):    
    if _verbosity_ > 0:
        print('Setting up commutative system')
    
    s = time.process_time()

    F = assumptions[0].parent()
    X = [str(v) for v in F.gens()]
    
    symbolic_vars = [f"{v}_{i}_{j}" for v in X if not v.endswith('_adj') for i in range(matrix_sizes[v][0]) for j in range(matrix_sizes[v][1])]
    rab_vars = [f"{_RABINOWITSCH_NAME_}{i}" for i in range(prod(matrix_sizes[f]))]
    tseitin_vars = [f"{_TSEITIN_NAME_}{i}" for i in range(prod(matrix_sizes[f]))]
    B = BooleanPolynomialRing(names=symbolic_vars + rab_vars + tseitin_vars)
    
    # make symbolic matrices
    if symbolic:
        symbolic_matrices = {v : create_symbolic_matrix(v, *matrix_sizes[v]) for v in X if not v.endswith('_adj')}
    else:
        symbolic_matrices = {v : create_symbolic_matrix(v, *matrix_sizes[v], B) for v in X if not v.endswith('_adj')}
        
    # make equations from assumptions    
    eqs = [evaluate(p, symbolic_matrices).dense_coefficient_list() for p in assumptions]
    eqs = [eq for eq_list in eqs for eq in eq_list]
   
    # make inequalities for claim
    ineqs = evaluate(f, symbolic_matrices).dense_coefficient_list()
    ineqs = [ineq for ineq in ineqs]
     
    # Rabinowitsch + Tseitin variables 
    rab_vars = [B(r) for r in rab_vars]
    tseitin_vars = [B(t) for t in tseitin_vars]
    eqs += [1 - r*ineq - t for t,r,ineq in zip(tseitin_vars, rab_vars, ineqs)]
    if symbolic:
        eqs.append( SR(prod(tseitin_vars)) )      
    else:
        eqs.append( prod(tseitin_vars) )  
    
    if _verbosity_ > 0:
        print(f'Commutative system set up in {time.process_time() - s:.4f}')

    return eqs, symbolic_matrices
    
##############################################################################

   
def run_cadical_with_cnf(cnf):
    """
    Run the system-installed CaDiCaL binary on a PySAT CNF object
    and return a model identical to PySAT's Cadical195.get_model().
    """
    
    s = time.process_time()
    
    def cadical_path():
        return Path(__file__).resolve().parent / "cadical"
        
    # Write CNF to a temporary DIMACS file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".cnf", delete=True) as tmp:
        cnf.write(tmp.name)

        # Run cadical on the file
        result = subprocess.run(
            [str(cadical_path()),"--seed=1234", tmp.name],
            capture_output=True,
            text=True
        )

    # Parse output
    status = None
    model = []
    for line in result.stdout.splitlines():
        if line.startswith("s "):
            status = line.split()[1]
        elif line.startswith("v "):
            lits = list(map(int, line.split()[1:]))
            if lits and lits[-1] == 0:
                lits = lits[:-1]
            model.extend(lits)
    
    if _verbosity_ > 0:
        if status == "SATISFIABLE":
            print(f'Found SAT solution in {time.process_time() - s:.4f}')
        else: 
            print(f'Found problem to be UNSAT in {time.process_time() - s:.4f}')
    if status == "SATISFIABLE":
        return model
    else:
        return None

##############################################################################

def construct_counterexample(assumptions,claim,Q=None,dims=None,n=2,random_tries=0,verbosity=1):
    r"""
    Input:
        - ``assumptions`` -- A list of noncommutative polynomials
        - ``claim`` -- a polynomial
        - ``Q`` -- a Quiver
        - ``dims`` -- a dictionary defining sizes of the matrices
        - ``n`` -- a natural number,indicates if we want to lift to 2^(2^n)
        - ``random_tries`` -- and integer number. If 0 we try Hensel lifting deterministic, if >0 we try hensel lifting not deterministic for random_tries number of                 tries
        - ``verbosity`` -- an integer number if >0 we output what the program does

    Output:
        - A counterexample in rational numbers if there exists one
    Example:
        sage: A.<x,y> = FreeAlgebra(QQ)
        sage: assumptions = [x*y]
        sage: claim = x
        sage: Q = Quiver([('U','U',x), ('U','U',y)])
        sage: dictionary = {'U':2}
        sage: n = 2
        sage: construct_counterexample(assumptions,claim,Q,dictionary,n)
    """
    start_time = time.process_time()
    
    global _verbosity_ 
    _verbosity_ = verbosity
    
    if Q is None and dims is not None:
        raise ValueError(f"Dimensions given but no quiver")
    if Q is None:
        Q, dims = generate_quiver(assumptions[0].parent().gens())    
    for g in assumptions + [claim]:
        if not Q.is_compatible(g):
            raise ValueError(f"{g} is not compatible with the quiver")
            
    # preprocessing - make copies for final check
    orig_assumptions = copy.deepcopy(assumptions)
    orig_claim = copy.deepcopy(claim)
    assumptions, claim, Q, dims = preprocess_input(assumptions, claim, Q, dims)
                    
    # generate commutative system
    eqs, symbolic_matrices = generate_commutative_system(assumptions, claim, Q, dims)
                    
    # convert system to CNF
    cnf, encoder = system_to_CNF(eqs)
    
    attempt = 0
    while True:
        sat_sol = run_cadical_with_cnf(cnf)
        if sat_sol is None:
            if _verbosity_ > 0:
                print("\nAll SAT solutions exhausted.")
                print('\nNo counterexample found.')
            return None
                        
        sat_dict = extract_sat_dict(sat_sol, encoder.phi[1:])                    
        sol = build_matrix_model(sat_sol, symbolic_matrices, sat_dict)
                        
        success = False
        if check_counterexample(orig_assumptions, sol, orig_claim):
            success = True            
        else: 
            if _verbosity_ > 0: 
                print('SAT solution NOT a counterexample. Needs lifting.') 
            
            if n == 2:
                sol = lift_mod_4(assumptions, sol, claim)
            else:
                sol = multi_lift_solution(assumptions, n, sol, claim, 1 + max(0, random_tries))
                  
            if sol is not None and check_counterexample(orig_assumptions, sol, orig_claim):
                if _verbosity_ > 0: 
                    print('Lifting was successful!')
                success = True
              
        if success:
            if _verbosity_ > 0:
                print(f'\nCounterexample found on attempt #{attempt + 1}')
                print_solution(sol, orig_claim.parent().gens())
                
            return sol
        else:
            if _verbosity_ > 0:
                print(f'Attempt #{attempt + 1}: Failed. Trying another SAT solution...')
            attempt += 1
            cnf.add_clause([-lit for lit in sat_sol])
 
##############################################################################            

def matrices_for_lifting(matrices, n, only_ones=False):
    new_matrices = dict()
    for i, (key, mat) in enumerate(matrices.items()):
        M = create_symbolic_matrix(f"{_LIFTING_NAME_}{i}", mat.nrows(), mat.ncols())
        if only_ones:
            for i in range(mat.nrows()):
                for j in range(mat.ncols()):
                    if mat[i,j] == 0:
                        M[i,j] = 0
        new_matrices[key] = mat + 2**n * M
    return new_matrices

##############################################################################            
    
def equations_for_lifting(matrices, n, X):
    ZX = PolynomialRing(ZZ, names=X)
    R = PolynomialRing(Integers(2**n), names=X)    
    
    eqs = set()
    for mat in matrices:
        if mat == 0: continue    
        mat = mat.change_ring(ZX) % 2**(2*n)
        mat = [R(str(c //  2**n)) for c in mat.coefficients()]
        eqs.update(mat) 
            
    return list(eqs)   

##############################################################################            

def solve_lifting_system(eqs):

    C,m = Sequence(eqs).coefficients_monomials()
    C = C.dense_matrix()    

    vars_in_sys = [str(var) for var in m if var != 1]
    is_1 = (m[-1] == 1)
        
    if is_1:
        b = -C[:, -1]              # RHS vector
        C = C[:, :-1]             # Coefficient matrix
        try:
            particular_sol = vector(C.solve_right(b))
        except:
            return None, None, None
    else:
        particular_sol = zero_vector(C.base_ring(),C.ncols())
    basis = C.right_kernel().gens()
                
    return particular_sol, basis, vars_in_sys
    
##############################################################################           

def extend_matrices(new_matrices, sol, all_vars, vars_in_sys, integer=False):
    vars_not_in_sys = {SR(v) : 0 for v in all_vars if v not in vars_in_sys}
  
    # plug into matrices
    var_to_val = {SR(var): int(val) for var, val in zip(vars_in_sys, sol)}
    var_to_val.update(vars_not_in_sys)
        
    lift = dict()        
    for key,mat in new_matrices.items():
        M = mat.subs(var_to_val).change_ring(ZZ)
        if integer:
            for i in range(M.nrows()):
                for j in range(M.ncols()):
                    if M[i,j] == 3:
                        M[i,j] = -1        
        lift[key] = M 
    return lift

##############################################################################          

def lift_solution(assumptions, n, matrices):

    new_matrices = matrices_for_lifting(matrices, n)    
    all_vars = [str(v) for M in new_matrices.values() for v in M.variables()]
                
    matrix_array = [evaluate(p, new_matrices) for p in assumptions]  
        
    eqs = equations_for_lifting(matrix_array, n, all_vars)
    if len(eqs) == 0:
        return None
    
    particular_sol, basis, vars_in_sys = solve_lifting_system(eqs)
    if particular_sol is None:
        return None
                  
    # make random linear combinations
    k = randint(0,len(basis))
    sol = particular_sol + sum(randint(0,2**(2*n-1)-1) * b for b in basis) 
                
    lift = extend_matrices(new_matrices, sol, all_vars, vars_in_sys)
    
    return lift
 
##############################################################################    

def lift_mod_4(assumptions, matrices, f=None):

    new_matrices = matrices_for_lifting(matrices,1,only_ones=True)
    all_vars = [str(v) for M in new_matrices.values() for v in M.variables()]
    
    matrix_array = [evaluate(p, new_matrices) for p in assumptions]  
        
    eqs = equations_for_lifting(matrix_array, 1, all_vars)  
    if len(eqs) == 0:
        return None

    particular_sol, basis, vars_in_sys = solve_lifting_system(eqs) 
    if particular_sol is None:
        return None  
        
    # make all linear combinations
    for subset in powerset(basis):
        sol = sum(subset) + particular_sol        
        lift = extend_matrices(new_matrices, sol, all_vars, vars_in_sys, integer=True)
        if check_counterexample(assumptions, lift, f):
            return lift
    
    return None
##############################################################################    
    
def multi_lift_solution(assumptions, n, matrices, f, random_tries=1):    
    s = time.process_time()
        
    if _verbosity_ > 0:
        print('Lifting solution...')
    
    failures = 0
                   
    for _ in range(random_tries):             
        lift = matrices
        good_lift = False
        for i in range(1, n):
            print(f'Lifting to 2^{2**i}...')
            lift = lift_solution(assumptions, 2**(i-1), lift)
            if lift is None:
                break
                                                                   
            mod = 2 ** (2**i)
            if evaluate(f,lift) % mod != 0:
                try: 
                    rat_lift = ratrec(lift, mod)
                    if check_counterexample(assumptions, rat_lift):
                        return rat_lift
                except Exception as e: 
                    pass                        
        
    return None

##############################################################################

def check_counterexample(assumptions, candidate, claim=None, Q=None):
    if Q is not None:
        for p in assumptions:
            if not Q.is_compatible(p):
                raise ValueError(f"{p} is not compatible with the quiver")
        if claim is not None:
            if not Q.is_compatible(claim):
                raise ValueError(f"{claim} is not compatible with the quiver")
       
    # Check assumptions: all must be zero
    if any(evaluate(p, candidate) != 0 for p in assumptions):
        return False
        
    if claim is None:
        return True
        
    return evaluate(claim, candidate) != 0

#takes the variable name, rows and cols and returns a mastrix of the given size where variables name depens on the input name
def create_symbolic_matrix(name, rows, cols, ring=None):
    #Create a symbolic matrix with entries named as name[i][j].
    if ring is None:
        ring = SR
    
    return Matrix(ring, rows, cols, lambda i, j: f"{name}_{i}_{j}")


def generate_quiver(X):
    Q = Quiver([('V','V',str(x)) for x in X])
    dims = {'V' : dim}
    return quiver, dims

def ratrec(lift, mod):
    return {
            key: Matrix(QQ, [
                [rational_reconstruction(matrix[i, j], mod) for j in range(matrix.ncols())]
                                                            for i in range(matrix.nrows())
            ])
            for key, matrix in lift.items()
        } 

def evaluate(p, matrices):
    r"""
    Input:
        - ``p`` -- a polynomial
        - ``matrix_dict`` -- a dictionary with matrices for each variable of the polynomial p (instead of dictionary of matrix sizes as in poly_to_matrix)

    Output:
        Evaluation of p under the given matrices

    """
    F = p.parent()
    X = [str(v) for v in F.gens()]
        
    subs_dict = dict()
    for v in X:
        if v.endswith('_adj'):
            subs_dict[F(v)] = matrices[v.removesuffix('_adj')].transpose()
        else:
            subs_dict[F(v)] = matrices[v]
                    
    # remove constant coefficient
    M = p.parent().basis().keys()
    c = p.coefficient(M(1))
    f = p - c 
                
    # plug in matrices               
    res = f.subs(subs_dict)
    # add constant coefficient back in
    if c != 0:
        res += c*identity_matrix(res.nrows())  
         
    return res

def build_matrix_model(model, symbolic_matrices, sat_dict):
    result = dict()
    for key, SM in symbolic_matrices.items():
        rows, cols = SM.nrows(), SM.ncols()
        M = zero_matrix(rows, cols)
        for i in range(rows):
            for j in range(cols):
                var_name = str(SM[i,j])
                if var_name in sat_dict:
                    M[i,j] = sat_dict[var_name]
                else:
                    M[i,j] = 0
        result[key] = M
    return result
    
#extracts important variables and its solution values
def extract_sat_dict(model, vars_list):
    return {
        str(var): (1 if model[i] > 0 else 0) 
        for i, var in enumerate(vars_list)
        if var != None and var.degree() == 1
    }

def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
    
    
def print_solution(sol,X):
    X = {str(v) for v in X}
    for key, matrix in sol.items():
        if str(key) in X:
            print(f'{key}: matrix([')
            for i,row in enumerate(matrix):
                if i < matrix.nrows()-1:
                    print(f'    {row},')
                else:
                    print(f'    {row}]),')