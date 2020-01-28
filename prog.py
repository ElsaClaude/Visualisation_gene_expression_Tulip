# Powered by Python 3.5
# To cancel the modifications performed by the script
# on the current gr, click on the undo button.
# Some useful keyboard shortcuts:
#   * Ctrl + D: comment selected lines.
#   * Ctrl + Shift + D: uncomment selected lines.
#   * Ctrl + I: indent selected lines.
#   * Ctrl + Shift + I: unindent selected lines.
#   * Ctrl + Return: run script.
#   * Ctrl + F: find selected text.
#   * Ctrl + R: replace selected text.
#   * Ctrl + Space: show auto-completion dialog.
from tulip import tlp
import pandas as pd
import csv

###### LES COMMENTAIRES : dire ce que fait la fonction, ce qu'elle prend en paramètre, ce qu'elle retourne

WD="/net/cremi/elclaude/espaces/travail/M2/DEA/Bourqui/projet/"

def create_interaction_graph(gr,data,viewLabel):
    nodes_dict={}
    for i in range(len(data["ID_locus1"])):
        if data["ID_locus1"][i] not in nodes_dict.keys():
            node_locus1 = gr.addNode()
            viewLabel[node_locus1] = data["ID_locus1"][i]
            nodes_dict[data["ID_locus1"][i]] = node_locus1     
        if data["ID_locus2"][i] not in nodes_dict.keys():
            node_locus2 = gr.addNode()
            viewLabel[node_locus2] = data["ID_locus2"][i]
            nodes_dict[data["ID_locus2"][i]] = node_locus2
        new_edge = gr.addEdge(nodes_dict[data["ID_locus1"][i]],nodes_dict[data["ID_locus2"][i]])
        gr.setEdgePropertiesValues(new_edge,{"Interaction":str(data["interaction_status"][i]),"Distance":str(data["distance"][i])})
    return nodes_dict 

def add_expression(gr,data,dico_nodes):
    for i in range(len(data["IDs"])):
        if data["IDs"][i] in dico_nodes.keys():
            gr.setNodePropertiesValues(dico_nodes[data["IDs"][i]],{"Expression":str(data["expression"][i])})

def read_symbols_csv(file):
    f = open(file,"r")
    dico = {}
    for line in f.readlines():
        line = line.split('\t')
        dico[line[0]]=line[2:]
    f.close()
    return dico

def voies_metaboliques(gr, viewLabel,data):
    for file in ["KEGG.symbols.csv","REACTOME.symbols.csv"]:
        voies_metabo = read_symbols_csv(WD+file)
        for (name, genes) in voies_metabo.items():
            intersection = list(set(genes) & set(data.keys()))
            print(intersection,">>>>>>>", name)
            if len(intersection) > 0 :
                currentSubgraph = gr.addSubGraph(name)
                for node_label in intersection:
                    gr.addNode(data[node_label])
                    for out_node in gr.getOutNodes(data[node_label]):
                        if viewLabel[out_node] in intersection:
                            currentSubgraph.addNode(out_node)
                            currentSubgraph.addEdge(data[node_label], out_node)

# The updateVisualization(centerViews = True) function can be called
# during script execution to update the opened views
# The pauseScript() function can be called to pause the script execution.
# To resume the script execution, you will have to click on the
# "Run script " button.
# The rungrScript(scriptFile, gr) function can be called to launch
# another edited script on a tlp.gr object.
# The scriptFile parameter defines the script name to call
# (in the form [a-zA-Z0-9_]+.py)
# The main(gr) function must be defined
# to run the script on the current gr

def main(gr):
    gr.clear()
    
    viewBorderColor = gr['viewBorderColor']
    viewBorderWidth = gr['viewBorderWidth']
    viewColor = gr['viewColor']
    viewFont = gr['viewFont']
    viewFontSize = gr['viewFontSize']
    viewIcon = gr['viewIcon']
    viewLabel = gr['viewLabel']
    viewLabelBorderColor = gr['viewLabelBorderColor']
    viewLabelBorderWidth = gr['viewLabelBorderWidth']
    viewLabelColor = gr['viewLabelColor']
    viewLabelPosition = gr['viewLabelPosition']
    viewLayout = gr['viewLayout']
    viewMetric = gr['viewMetric']
    viewRotation = gr['viewRotation']
    viewSelection = gr['viewSelection']
    viewShape = gr['viewShape']
    viewSize = gr['viewSize']
    viewSrcAnchorShape = gr['viewSrcAnchorShape']
    viewSrcAnchorSize = gr['viewSrcAnchorSize']
    viewTexture = gr['viewTexture']
    viewTgtAnchorShape = gr['viewTgtAnchorShape']
    viewTgtAnchorSize = gr['viewTgtAnchorSize']
     
    interaction_data=pd.read_csv(WD+"interactions_chromosome6.csv",sep="\t",header=0)
    expression_data=pd.read_csv(WD+"chromosome6_fragments_expressions.csv", sep="\t", header=0)
    
    dico_nodes=create_interaction_graph(gr,interaction_data,viewLabel)
    add_expression(gr,expression_data,dico_nodes)
    updateVisualization(centerViews = True)
    
    print("OK interactions")
    
    voies_metaboliques(gr,viewLabel,dico_nodes)
    print("OK voies métabo")
    
  
