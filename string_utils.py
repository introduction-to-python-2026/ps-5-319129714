import numpy as np

def split_before_uppercases(formula):
    parts = []
    current_part = ""
    for char in formula:
        if char.isupper() and current_part:
            parts.append(current_part)
            current_part = char
        else:
            current_part += char
    if current_part:
        parts.append(current_part)
    return parts


def split_at_digit(formula):
    element = ""
    number_str = ""
    for char in formula:
        if char.isdigit():
            number_str += char
        else:
            element += char
    
    if number_str == "":
        count = 1
    else:
        count = int(number_str)
        
    return element, count


def count_atoms_in_molecule(molecular_formula):
    atom_counts = {}
    for atom_segment in split_before_uppercases(molecular_formula):
        atom_name, atom_count = split_at_digit(atom_segment)
        if atom_name in atom_counts:
            atom_counts[atom_name] += atom_count
        else:
            atom_counts[atom_name] = atom_count
    return atom_counts


def balance_reaction(equation):
    """
    Balances a chemical equation using linear algebra.
    Example: "Fe2O3 + H2 -> Fe + H2O" -> [0.333, 1.0, 0.666, 1.0]
    """
    # 1. Parse the equation into reactants and products list
    reactants_str, products_str = equation.split('->')
    reactants = [m.strip() for m in reactants_str.split('+')]
    products = [m.strip() for m in products_str.split('+')]
    all_molecules = reactants + products
    
    # 2. Get atom counts for every molecule
    molecule_compositions = [count_atoms_in_molecule(m) for m in all_molecules]
    
    # 3. Find all unique elements in the reaction
    all_elements = set()
    for m in molecule_compositions:
        all_elements.update(m.keys())
    all_elements = sorted(list(all_elements))
    
    # 4. Build the Matrix for linear equation system
    # Rows = elements, Columns = molecules
    matrix = []
    for el in all_elements:
        row = []
        for i, mol_counts in enumerate(molecule_compositions):
            count = mol_counts.get(el, 0)
            # Products (right side) need to be negative to satisfy (Reactants - Products = 0)
            if i >= len(reactants):
                count = -count
            row.append(count)
        matrix.append(row)
    
    # 5. Solve using NumPy
    # We assume the coefficient of the LAST molecule is 1.0, and solve for the rest.
    # Ax = 0  =>  A_partial * x_partial = - (Last_Column * 1)
    
    matrix_np = np.array(matrix)
    
    # Split into coefficients to solve (A) and the dependent column (b)
    a_matrix = matrix_np[:, :-1]
    b_vector = -matrix_np[:, -1]
    
    # Use Least Squares to solve (works even if matrix is not square)
    coeffs, residuals, rank, s = np.linalg.lstsq(a_matrix, b_vector, rcond=None)
    
    # Add the last coefficient (which we fixed as 1.0)
    result = list(coeffs) + [1.0]
    
    return result
