import os
import re


def uppercase_bnd_node_names(file_path):
    with open(file_path, "r") as f:
        content = f.read()

    # 1. Uppercase node name after 'Node'
    content = re.sub(
        r"Node\s+([a-zA-Z0-9_]+)\s*\{",
        lambda m: f"Node {m.group(1).upper()} {{",
        content,
    )

    # 2. Uppercase all node names in logic expressions
    def logic_upper(match):
        expr = match.group(1)
        expr_up = re.sub(
            r"\b([a-zA-Z_][a-zA-Z0-9_]*)\b", lambda m: m.group(1).upper(), expr
        )
        return f"logic = {expr_up};"

    content = re.sub(r"logic\s*=\s*(.*?);", logic_upper, content)

    # 3. Uppercase all after the first underscore for $u_XXX, $d_XXX, u_XXX, d_XXX
    content = re.sub(
        r"(\$?[ud]_)([a-zA-Z0-9_]+)", lambda m: m.group(1) + m.group(2).upper(), content
    )

    with open(file_path, "w") as f:
        f.write(content)


def uppercase_cfg_node_names(file_path):
    with open(file_path, "r") as f:
        content = f.read()
    # Uppercase after $u_ or $d_ (with or without $)
    content = re.sub(
        r"(\$[ud]_)([a-zA-Z0-9_]+)", lambda m: m.group(1) + m.group(2).upper(), content
    )
    # Uppercase before .is_internal
    content = re.sub(
        r"([a-zA-Z0-9_]+)(\.is_internal)",
        lambda m: m.group(1).upper() + m.group(2),
        content,
    )
    # Uppercase before .istate
    content = re.sub(
        r"([a-zA-Z0-9_]+)(\.istate)", lambda m: m.group(1).upper() + m.group(2), content
    )
    # Uppercase inside brackets before .istate (for commented lines)
    content = re.sub(
        r"(\[)([a-zA-Z0-9_]+)(\]\.istate)",
        lambda m: m.group(1) + m.group(2).upper() + m.group(3),
        content,
    )
    with open(file_path, "w") as f:
        f.write(content)


# insert the new nodes
def create_new_bnd_nodes(file_path, original_node, new_node):
    with open(file_path, "r") as f:
        content = f.read()

    # Find the node block for the original node
    pattern = re.compile(
        rf"(Node\s+{original_node}\s*\{{.*?\n\}})", re.DOTALL | re.IGNORECASE
    )
    match = pattern.search(content)
    if not match:
        print(f"Node {original_node} not found.")
        return

    node_block = match.group(0)

    # Replace all whole-word occurrences of original_node (including $u_ and $d_ forms)
    def node_repl(m):
        # Handle $u_ORIGINAL, $d_ORIGINAL, u_ORIGINAL, d_ORIGINAL
        prefix = m.group(1) or ""
        return f"{prefix}{new_node}"

    # Replace $u_ORIGINAL, $d_ORIGINAL, u_ORIGINAL, d_ORIGINAL
    node_block_new = re.sub(
        rf"(\$?[ud]_){original_node}\b", rf"\1{new_node}", node_block
    )
    # Replace all other whole-word occurrences (e.g., in logic)
    node_block_new = re.sub(rf"\b{original_node}\b", new_node, node_block_new)

    # Insert the new node block after the original
    insert_pos = match.end()
    new_content = content[:insert_pos] + "\n\n" + node_block_new + content[insert_pos:]

    with open(file_path, "w") as f:
        f.write(new_content)


def create_new_cfg_nodes(file_path, original_node, new_node):
    with open(file_path, "r") as f:
        lines = f.readlines()

    new_lines = []
    # Patterns to match: $u_ORIGINAL, $d_ORIGINAL, ORIGINAL.is_internal
    patterns = [
        rf"\$u_{original_node}\b",
        rf"\$d_{original_node}\b",
        rf"\b{original_node}\.is_internal\b",
    ]
    for line in lines:
        new_lines.append(line)
        if any(re.search(pat, line) for pat in patterns):
            # Replace all occurrences of original_node with new_node in the line
            new_line = re.sub(rf"\b{original_node}\b", new_node, line)
            # Also handle $u_ORIGINAL and $d_ORIGINAL
            new_line = re.sub(
                rf"(\$[ud]_)" + original_node + r"\b", rf"\1" + new_node, new_line
            )
            new_lines.append(new_line)
    with open(file_path, "w") as f:
        f.writelines(new_lines)


def add_new_node_to_logic(file_path, original_node, new_node):
    with open(file_path, "r") as f:
        content = f.read()

    # Function to replace in logic expressions
    def logic_repl(match):
        expr = match.group(1)
        # Replace whole-word original_node with (original_node | new_node)
        expr_new = re.sub(
            rf"\b{original_node}\b", f"({original_node} | {new_node})", expr
        )
        return f"logic = {expr_new};"

    # Replace in all logic expressions
    content = re.sub(r"logic\s*=\s*(.*?);", logic_repl, content)

    with open(file_path, "w") as f:
        f.write(content)


# Main function
def replace_node_names_in_file(file_path, name_maps):
    ext = os.path.splitext(file_path)[1]
    if ext == ".cfg":
        uppercase_cfg_node_names(file_path)
    else:
        uppercase_bnd_node_names(file_path)

    # transform all the nodes names to genes/ proteins names depending on the type_models

    with open(file_path, "r") as f:
        content = f.read()
    for old_name, new_name in name_maps.items():
        # Use regex for case-insensitive replacement
        pattern = re.compile(re.escape(old_name), re.IGNORECASE)
        content = pattern.sub(new_name, content)

    with open(file_path, "w") as f:
        f.write(content)

    # create new nodes: montagud_nodes += ['MEK2', 'TSC2', 'MAP3K3', 'CHK2']

    nodes_add = {"MEK1": "MEK2", "TSC1": "TSC2", "MAP3K1": "MAP3K3", "CHK1": "CHK2"}

    for key, value in nodes_add.items():
        if ext == ".cfg":
            create_new_cfg_nodes(file_path, key, value)
            add_new_node_to_logic(file_path, key, value)
        else:
            create_new_bnd_nodes(file_path, key, value)
            add_new_node_to_logic(file_path, key, value)
