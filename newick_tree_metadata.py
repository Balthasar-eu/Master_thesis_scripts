color = {
"ERR1100909" : "orange",
"ERR1100911" : "orange",
"ERR1100912" : "orange",
"ERR1100913" : "orange",
"ERR1100923" : "orange",
"ERR1100924" : "orange",
"ERR1100925" : "orange",
"ERR1100926" : "orange",
"ERR1100928" : "orange",
"ERR1100930" : "orange",
"ERR1100935" : "orange",
"ERR1100936" : "orange",
"ERR1100937" : "orange",
"ERR1100938" : "orange",
"ERR1100939" : "orange",
"ERR1100940" : "orange",
"ERR1100941" : "orange",
"ERR1100942" : "orange",
"ERR1100944" : "orange",
"ERR1100946" : "orange",
"ERR1100947" : "orange",
"ERR1100949" : "orange",
"ERR1100950" : "orange",
"ERR1100952" : "orange",
"ERR1100954" : "orange",
"ERR1100955" : "orange",
"ERR1100972" : "orange",
"ERR1100973" : "orange",
"ERR1100975" : "orange",
"ERR1100976" : "orange",
"10CEB535LM" : "orange",
"10CEB540LM" : "orange",
"10CEB550LM" : "orange",
"10CEB552LM" : "orange",
"10CEB553LM" : "orange",
"10CEB554LM" : "orange",
"10CEB559LM" : "orange",
"10CEB560LM" : "orange",
}

Arla_results = [2, 3, 2, 3, 2, 2, 2, 2, 2, 2, 3, 2, 2, 4, 2, 2, 2, 2, 2, 4, 2, 3, 4, 2, 3, 4, 2, 2, 2, 2, 3, 3, 2, 2, 1, 3, 3,
 3, 2, 2, 2, 3, 2, 3, 3, 3, 3, 2, 2, 3]

Arla_country = [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

Arla_product = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1]

colordict = {0:"black", 1:"blue", 2:"green",3:"orange",4:"red"}

arla_class_color = {str(i):colordict[c] for i, c in enumerate(Arla_country,1)}

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def label_func(string):
    string = str(string)
    if string[:3] == "DTU":
        return string[30:]
    else:
        return string


from Bio import Phylo

fname = "data/arla_only_accessory_binary_genes.fa.newick"

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111)

patches = [mpatches.Patch(color=c, label=str(i)) for i, c in zip(["Not specified", "Germany", "Italy", "Finland"], ["black", "blue", "green", "orange"])]

ax.legend(handles=patches)

tree = Phylo.read(fname, "newick")
tree.root_with_outgroup({'name': '24'})

Phylo.draw(tree, label_func=label_func, label_colors=arla_class_color, axes=ax, show_confidence=False)

fig.savefig("plots/tree_arla_only_metadata_country.pdf", bbox_inches = 'tight', pad_inches = 0)

arla_class_color = {str(i):colordict[c] for i, c in enumerate(Arla_product,1)}

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111)

patches = [mpatches.Patch(color=c, label=str(i)) for i, c in zip(["Not specified", "Dairy/Cheese", "Plant based"], ["black", "blue", "green"])]

ax.legend(handles=patches)

tree2 = Phylo.read(fname, "newick")
tree2.root_with_outgroup({'name': '24'})

Phylo.draw(tree2, label_func=label_func, label_colors=arla_class_color, axes=ax, show_confidence=False)

fig.savefig("plots/tree_arla_only_metadata_product.pdf", bbox_inches = 'tight', pad_inches = 0)

