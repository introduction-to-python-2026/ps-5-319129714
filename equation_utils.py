import numpy as np
from string_utils import count_atoms_in_molecule 

def parse_chemical_reaction(reaction_equation):
    reaction_equation = reaction_equation.replace(" ", "")
    reactants, products = reaction_equation.split("->")
    return reactants.split("+"), products.split("+")


def count_atoms_in_reaction(molecules_list):
    molecules_atoms_count = []
    for molecule in molecules_list:
        molecules_atoms_count.append(count_atoms_in_molecule(molecule))
    return molecules_atoms_count


def balance_reaction(equation):
    
    reactants, products = parse_chemical_reaction(equation)
    all_molecules = reactants + products
    
    molecule_compositions = [count_atoms_in_molecule(m) for m in all_molecules]
    
    all_elements = set()
    for m in molecule_compositions:
        all_elements.update(m.keys())
    all_elements = sorted(list(all_elements))
    
    matrix = []
    for el in all_elements:
        row = []
        for i, mol_counts in enumerate(molecule_compositions):
            count = mol_counts.get(el, 0)
            if i >= len(reactants):
                count = -count
            row.append(count)
        matrix.append(row)
    
    matrix_np = np.array(matrix)
    
    a_matrix = matrix_np[:, :-1]
    b_vector = -matrix_np[:, -1]
    
    coeffs, residuals, rank, s = np.linalg.lstsq(a_matrix, b_vector, rcond=None)
    
    result = list(coeffs) + [1.0]
    
    return result

ELEMENTS = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'I', 'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo'
]

def generate_equation_for_element(compounds, coefficients, element):
    """Generates a symbolic equation for the given element from compounds and coefficients.  
    Example: For H in reactants [{'H': 2}, {'O': 4, 'H': 1}], coefficients [a0, a1], returns 2*a0 + a1."""
    equation = 0
    for i, compound in enumerate(compounds):
        if element in compound:
            equation += coefficients[i] * compound[element]
    return equation


def build_equations(reactant_atoms, product_atoms):
    """Builds a list of symbolic equations for each element to balance a chemical reaction.  
    Example: For H2 + O2 -> H2O, returns equations [2*a0 - 2*b0, a1 - b0]."""
    ## coefficients ##
    reactant_coefficients = list(symbols(f'a0:{len(reactant_atoms)}'))
    product_coefficients = list(symbols(f'b0:{len(product_atoms)}')) 
    product_coefficients = product_coefficients[:-1] + [1] # Ensure the last coefficient is 1

    ## equations ##
    equations = []
    for element in ELEMENTS:
        lhs = generate_equation_for_element(reactant_atoms, reactant_coefficients, element)
        rhs = generate_equation_for_element(product_atoms, product_coefficients, element)
        if lhs != 0 or rhs != 0:
            equations.append(Eq(lhs, rhs))

    return equations, reactant_coefficients + product_coefficients[:-1]


def my_solve(equations, coefficients):
    """Solves the system of equations for the coefficients of the reaction.  
    Example: For equations [2*a0 - 2*b0, a1 - b0], returns [1.0, 1.0]."""
    solution = sympy_solve(equations, coefficients)

    if len(solution) == len(coefficients):
        coefficient_values = list()
        for coefficient in coefficients:
            coefficient_values.append(float(solution[coefficient]))
        return coefficient_values





