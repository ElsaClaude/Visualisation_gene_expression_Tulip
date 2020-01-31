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
import os

###### LES COMMENTAIRES : dire ce que fait la fonction, ce qu'elle prend en paramètre, ce qu'elle retourne

#Amelie's path
#WD="/home/amelie/Documents/master/M2/DEA/Tulip/Visualisation_gene_expression_Tulip/"

#Elsa's path
#WD=""

#Salome's path
#WD=""

#Antoine's path
#WD="/net/cremi/alaporte006/espaces/travail/DEA_Bourqui/Visualisation_gene_expression_Tulip/"
#WD="~/Dropbox/Master/M2S2/DEA/R_Bourqui/Visualisation_gene_expression_Tulip/"
#WDopen=os.getcwd()+"/../Dropbox/Master/M2S2/DEA/R_Bourqui/Visualisation_gene_expression_Tulip/"
#WDopenDos=os.getcwd()+"\\..\\..\\Users\\antoi\\Dropbox\\Master\\M2S2\\DEA\\R_Bourqui\\Visualisation_gene_expression_Tulip\\"


def create_interaction_graph(gr,data,viewLabel):
    print("\nInteraction graph being constructed")
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
    print("\nAdding expression to graph")
    for i in range(len(data["IDs"])):
        if data["IDs"][i] in dico_nodes.keys():
            gr.setNodePropertiesValues(dico_nodes[data["IDs"][i]],{"Expression":str(data["expression"][i])})

def visu_algoFM(gr):
    params = tlp.getDefaultPluginParameters("FM^3 (OGDF)",gr)
    params["Unit Edge Length"] = gr["Distance"]
    gr.applyLayoutAlgorithm("FM^3 (OGDF)", params)
    
    params = tlp.getDefaultPluginParameters("Perfect aspect ratio",gr)
    gr.applyLayoutAlgorithm('Perfect aspect ratio', params)

def read_symbols_csv(files_symbols):
    dico = {}
    for file in files_symbols:
      f = open(WD+file,"r")
      for line in f.readlines():
          line = line.split('\t')
          dico[line[0]]=line[2:]
      f.close()
    return dico

def get_subgraphs_pathways(gr, viewLabel,data):
    print("\nCreation of subgraph for each pathway")
    voies_metabo = read_symbols_csv(["KEGG.symbols.csv","REACTOME.symbols.csv"])
    for (name, genes) in voies_metabo.items():
      intersection = list(set(genes) & set(data.keys()))
      if len(intersection) > 0 :
        currentSubgraph = gr.addSubGraph(name)
        edges_to_add = []
        for node_label in intersection:
          currentSubgraph.addNode(data[node_label])
          for edge in gr.getOutEdges(data[node_label]):
            if viewLabel[gr.target(edge)] in intersection:
              edges_to_add.append(edge)
            elif gr.getEdgePropertiesValues(edge)["Interaction"] == "gain" or gr.getEdgePropertiesValues(edge)["Interaction"] == "loss":
              currentSubgraph.addNode(gr.target(edge))
              edges_to_add.append(edge)
        currentSubgraph.addEdges(edges_to_add)

def node_custom(gr,dico_nodes,size,color):
    aspect = {
        "up": [tlp.Color(0,255,0),tlp.Size(2,2,2)],
        "down": [tlp.Color(255,0,0),tlp.Size(2,2,2)],
        "stable": [tlp.Color(105,105,105),tlp.Size(1,1,1)],
        "intergenic": [tlp.Color(200,200,200),tlp.Size(1,1,1)],
        "nan": [tlp.Color(255,255,255),tlp.Size(1,1,1)]
    }
    for node in dico_nodes.keys():
        color[dico_nodes[node]] = aspect[gr.getNodePropertiesValues(dico_nodes[node])["Expression"]][0]
        size[dico_nodes[node]] = aspect[gr.getNodePropertiesValues(dico_nodes[node])["Expression"]][1]

def visu_Edges(gr,viewBorderColor, viewBorderWidth, viewColor):
    interaction = gr["Interaction"]
    aspect = {
      "gain": [tlp.Color.Blue,10],
      "loss": [tlp.Color.Yellow,10],
      "stable": [tlp.Color.Gray,0]
    }
    for edge in gr.getEdges():
      viewBorderColor[edge] = aspect[interaction[edge]][0]
      viewBorderWidth[edge] = aspect[interaction[edge]][1]
      viewColor[edge] = aspect[interaction[edge]][0]

def create_interest_subgraph(gr):
    interest = gr.addSubGraph("Graph of interest")
    nodes_to_add = []
    edges_to_add = []
    for node in gr.getNodes():
        if gr["Expression"][node] in ["up","down"]:
            nodes_to_add.append(node)
        elif gr["Expression"][node] == "intergenic":
            count = 0
            for neighbor in gr.getInOutNodes(node):
                if gr["Expression"][neighbor] in ["up","down"]:
                    count+=1
            if count >=2:
                nodes_to_add.append(node)
    interest.addNodes(nodes_to_add)
    allSubNodes = interest.nodes()
    for edge in gr.getEdges():
        if gr["Interaction"][edge] in ["gain","loss","stable"]:
            if gr.target(edge) in allSubNodes and gr.source(edge) in allSubNodes:
                edges_to_add.append(edge)
    interest.addEdges(edges_to_add)

def get_regulators(gr,viewLabel):
  for node in gr.getNodes():
    if gr["Expression"][node] in ["up","down"] and list(set([gr["Interaction"][e] for e in list(gr.getInOutEdges(node))])) == ["stable"]:
      for neighbor_node in gr.getInOutNodes(node):
        if gr["Expression"][neighbor_node] in ["up","down","intergenic"]:
          for second_degree_edge in gr.getInOutEdges(neighbor_node):
            second_degree_neighbor_node = [n for n in gr.ends(second_degree_edge) if n != neighbor_node][0] 
            if gr["Expression"][second_degree_neighbor_node] in ["up","down"] and gr["Interaction"][second_degree_edge] != "stable":
              print(viewLabel[second_degree_neighbor_node], "influence", viewLabel[node],"via",viewLabel[neighbor_node])
              gr.setNodePropertiesValues(neighbor_node,{"Regulators": str(viewLabel[node])+" par "+str(viewLabel[second_degree_neighbor_node])})

def get_statistics(gr, viewLabel):
  statistics = {"genes": {}, "interactions": {}}
  genes_in_pathways = {}
  indirect_regulators = {}
  pathways_from_files = read_symbols_csv(["KEGG.symbols.csv", "REACTOME.symbols.csv"])
  genes_in_pathways_from_file = list(set(sum(list(pathways_from_files.values()), [])))
  print("read csv OK")
  i = 1
  
  for node in gr.getNodes(): 
    if gr["Expression"][node] in statistics["genes"].keys():
      statistics["genes"][gr["Expression"][node]] += 1
    else : 
      statistics["genes"][gr["Expression"][node]] = 1
      
    if viewLabel[node] in genes_in_pathways_from_file and gr["Expression"][node] in ["up","down"] :
      genes_in_pathways[viewLabel[node]] = []
      for (pathway, genes) in pathways_from_files.items():
        if viewLabel[node] in genes :
          genes_in_pathways[viewLabel[node]].append(pathway)
    print(i, gr.numberOfNodes())
    i+=1
    
  print("nodes OK")
  for edge in gr.getEdges():
    if gr["Interaction"][edge] in statistics["interactions"].keys():
      statistics["interactions"][gr["Interaction"][edge]] += 1
    else : 
      statistics["interactions"][gr["Interaction"][edge]] = 1
  print("edges OK")
  return statistics, genes_in_pathways

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
    
    #Reading of the interaction and expression files
    interaction_data=pd.read_csv(WD+"interactions_chromosome6.csv",sep="\t",header=0)
    expression_data=pd.read_csv(WD+"chromosome6_fragments_expressions.csv", sep="\t", header=0)
#    
#    #Functions to create the graph
    dico_nodes=create_interaction_graph(gr,interaction_data,viewLabel)
    add_expression(gr,expression_data,dico_nodes)
    updateVisualization(centerViews = True)
#    
#    ### temporaire : applique automatique FM³
    visu_algoFM(gr)
#    
    print("\nGraph constructed successfully")
#    
#    #Creating of subgraph for each pathway
#    get_subgraphs_pathways(gr,viewLabel,dico_nodes)
#    print("\nMetabolism done")

    #Customization of the nodes regarding their properties
    node_custom(gr,dico_nodes,viewSize,viewColor)
    visu_Edges(gr, viewBorderColor,viewBorderWidth, viewColor)
    
    create_interest_subgraph(gr)
    
#    statistics, genes_in_pathways = get_statistics(gr, viewLabel)
#    print(statistics)
#    print(genes_in_pathways)
    
#    get_regulators(gr,viewLabel)
#    print("OK")
    
#    for edge in gr.getEdges():
#      labels = [viewLabel[n] for n in gr.ends(edge)]
#      if "MCM3" in labels and ("IL17A" in labels or "EFHC1" in labels):
#        viewColor[edge] = tlp.Color.Blue
    
