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


def count_atoms_in_reaction(molecules_list):
    molecules_atoms_count = []
    for molecule in molecules_list:
        molecules_atoms_count.append(count_atoms_in_molecule(molecule))
    return molecules_atoms_count
