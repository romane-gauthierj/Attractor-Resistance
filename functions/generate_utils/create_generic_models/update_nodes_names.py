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

    # Check if the new node already exists to avoid duplicates
    new_node_pattern = re.compile(rf"Node\s+{re.escape(new_node)}\s*\{{", re.IGNORECASE)
    if new_node_pattern.search(content):
        print(f"Node {new_node} already exists. Skipping creation.")
        return
    

    # Find the node block for the original node
    pattern = re.compile(
        rf"(Node\s+{original_node}\s*\{{.*?\n\}})", re.DOTALL | re.IGNORECASE
    )
    match = pattern.search(content)
    if not match:
        print(f"Node {original_node} not found.")
        return

    node_block = match.group(0)

   # Create new node block by replacing node name and variable references
    node_block_new = re.sub(
        rf"Node\s+{re.escape(original_node)}\s*\{{", 
        f"Node {new_node} {{", 
        node_block
    )
    
    # Replace $u_ and $d_ references
    node_block_new = re.sub(
        rf"(\$[ud]_){re.escape(original_node)}\b", 
        rf"\1{new_node}", 
        node_block_new
    )


    # Replace all other whole-word occurrences (e.g., in logic)
    # node_block_new = re.sub(rf"\b{original_node}\b", new_node, node_block_new)


    # Insert the new node block after the original
    insert_pos = match.end()
    new_content = content[:insert_pos] + "\n\n" + node_block_new + content[insert_pos:]

    with open(file_path, "w") as f:
        f.write(new_content)
    
    print(f"Created new node {new_node} based on {original_node}")





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


def remove_cfg_node(file_path, node_to_remove):
    with open(file_path, "r") as f:
        lines = f.readlines()

    new_lines = []

    patterns_to_remove = [
        rf"\$u_{node_to_remove}\b",
        rf"\$d_{node_to_remove}\b",
        rf"\b{node_to_remove}\.is_internal\b",
        rf"^\s*//\s*\[{node_to_remove}\]\.istate=\s*(?:\d+\s*\[\d+\](?:\s*,\s*\d+\s*\[\d+\])*)?\s*;?\s*$",
    ]

    for line in lines:
        # Check if the line contains any of the patterns related to the node to remove.
        # If it does, we skip adding this line to new_lines, effectively removing it.
        if not any(re.search(pat, line) for pat in patterns_to_remove):
            cleaned_line = line

            cleaned_line = re.sub(rf"\s*\$u_{node_to_remove}\b", "", cleaned_line)
            cleaned_line = re.sub(rf"\s*\$d_{node_to_remove}\b", "", cleaned_line)

            new_lines.append(cleaned_line)

    with open(file_path, "w") as f:
        f.writelines(new_lines)


# Main function
def replace_node_names_in_file(file_path, name_maps, nodes_to_add):
    ext = os.path.splitext(file_path)[1]
    if ext == ".cfg":
        uppercase_cfg_node_names(file_path)
    else:
        uppercase_bnd_node_names(file_path)



    with open(file_path, "r") as f:
        content = f.read()
    for old_name, new_name in name_maps.items():
        # Use regex for case-insensitive replacement
        pattern = re.compile(re.escape(old_name), re.IGNORECASE)
        content = pattern.sub(new_name, content)

    with open(file_path, "w") as f:
        f.write(content)

    # Add duplicate tracking here
    if nodes_to_add:
        created_nodes = set()
        for key, value in nodes_to_add.items():
            if value not in created_nodes:
                if ext == ".cfg":
                    create_new_cfg_nodes(file_path, key, value)
                    add_new_node_to_logic(file_path, key, value)
                else:
                    create_new_bnd_nodes(file_path, key, value)
                    add_new_node_to_logic(file_path, key, value)
                created_nodes.add(value)
            else:
                print(f"Node {value} already created, skipping...")

