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


def parse_chemical_reaction(reaction_equation):
    reaction_equation = reaction_equation.replace(" ", "")
    reactants, products = reaction_equation.split("->")
    return reactants.split("+"), products.split("+")


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


def count_atoms_in_reaction(molecules_list):
    molecules_atoms_count = []
    for molecule in molecules_list:
        molecules_atoms_count.append(count_atoms_in_molecule(molecule))
    return molecules_atoms_count
