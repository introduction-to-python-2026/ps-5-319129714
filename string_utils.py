


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
def count_atoms_in_molecule(molecular_formula):
    atom_counts = {}
    formula = split_before_uppercases(molecular_formula)

    for formula in formula:
        element, count = split_at_digit(formula)
        if element in atom_counts:
            atom_counts[element] += count
        else:
            atom_counts[element] = count

    return atom_counts


def parse_chemical_reaction(reaction_equation):
    """Takes a reaction equation (string) and returns reactants and products as lists.  
    Example: 'H2 + O2 -> H2O' → (['H2', 'O2'], ['H2O'])"""
    reaction_equation = reaction_equation.replace(" ", "")  # Remove spaces for easier parsing
    reactants, products = reaction_equation.split("->")
    return reactants.split("+"), products.split("+")

def count_atoms_in_reaction(molecules_list):
    """Takes a list of molecular formulas and returns a list of atom count dictionaries.  
    Example: ['H2', 'O2'] → [{'H': 2}, {'O': 2}]"""
    molecules_atoms_count = []
    for molecule in molecules_list:
        molecules_atoms_count.append(count_atoms_in_molecule(molecule))
    return molecules_atoms_count
