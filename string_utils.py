def split_before_uppercases(formula):
    split_formula = []
    start = 0
    for i in range(1, len(formula)):
        if formula[i].isupper():
            split_formula.append(formula[start:i])
            start = i
    if formula:
        split_formula.append(formula[start:])
    return split_formula

def split_at_digit(formula):
    for i, char in enumerate(formula):
        if char.isdigit():
            prefix = formula[:i]
            numeric_part = int(formula[i:])
            return (prefix, numeric_part)
    return (formula, 1)

def count_atoms_in_molecule(molecular_formula):
    dict_of_atoms = {}
    split_formula = split_before_uppercases(molecular_formula)
    for atom in split_formula:
        atom_name, atom_count = split_at_digit(atom)
        dict_of_atoms[atom_name] = atom_count
    return dict_of_atoms
