import sys
import os
import toytree
import toyplot.svg

def remove_single_child_nodes(tree):
    """
    Removes nodes with only one child from a phylogenetic tree.

    Args:
        tree (Toytree.Tree): The input tree from which to remove single-child nodes.

    Returns:
        Toytree.Tree: The tree with single-child nodes removed.
    """
    while True:
        singles = [node.idx for node in tree.idx_dict.values() if len(node.children) == 1]
        if not singles:
            break
        tree = tree.drop(treenodes=singles, axis=0)
    return tree

def tree_generator(infile, outfile):
    """
    Generates a phylogenetic tree from a Newick string and saves it as an SVG file.

    Args:
        infile (str): Path to the input file containing the Newick formatted tree string.
        outfile (str): Path to save the output SVG file.
    """
    # Check if input file exists
    if not os.path.isfile(infile):
        raise FileNotFoundError(f"Input file '{infile}' not found.")

    # Read the Newick tree string from the input file
    with open(infile, 'r') as f:
        newick = f.read().strip()

    # Create a Toytree object from the Newick string
    tree = toytree.tree(newick)

    # Root the tree with specified outgroup species
    outgroup_species = ["Macrobiotus_rybaki"]
    tree = tree.root(names=outgroup_species)

    # Remove nodes with only one child
    tree = remove_single_child_nodes(tree)

    # Define plot dimensions based on the number of taxa
    num_taxa = len(tree.get_tip_labels())
    width = max(2000, num_taxa * 50)
    height = max(1000, num_taxa * 30)

    # Get support values for internal nodes
    support_values = tree.get_node_values('support')
    labels = [f"{round(float(s), 2)}" if s else "" for s in support_values]

    # Create and save the phylogenetic tree plot
    canvas, axes, mark = tree.draw(node_labels=labels, node_sizes=30, width=width, height=height, tree_style="n")
    toyplot.svg.render(canvas, outfile)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    tree_generator(input_file, output_file)
