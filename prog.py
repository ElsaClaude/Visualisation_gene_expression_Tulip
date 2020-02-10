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
WD="/autofs/unitytravail/travail/agruel/M2/DEA/tulip/Visualisation_gene_expression_Tulip/"

#Elsa's path
#WD=""

#Salome's path
#WD=""

#Antoine's path
#WD="/net/cremi/alaporte006/espaces/travail/DEA_Bourqui/Visualisation_gene_expression_Tulip/"
#WD="~/Dropbox/Master/M2S2/DEA/R_Bourqui/Visualisation_gene_expression_Tulip/"
#WDopen=os.getcwd()+"/../Dropbox/Master/M2S2/DEA/R_Bourqui/Visualisation_gene_expression_Tulip/"
#WDopenDos=os.getcwd()+"\\..\\..\\Users\\antoi\\Dropbox\\Master\\M2S2\\DEA\\R_Bourqui\\Visualisation_gene_expression_Tulip\\"


def create_interaction_graph(gr,viewLabel):
    data = pd.read_csv(WD+"interactions_chromosome6.csv",sep="\t",header=0)
    print("\nInteraction graph being constructed")
    nodes_dict = {}
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

def add_expression(gr,dico_nodes):
    data = pd.read_csv(WD+"chromosome6_fragments_expressions.csv", sep="\t", header=0)
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

def set_subgraphs_pathways(gr, viewLabel,data):
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

def visu_node_edge(gr,size,color,viewBorderColor,viewBorderWidth):
    interaction = gr["Interaction"]
    expression = gr["Expression"]
    aspect = {"node" : {
        "up": [tlp.Color(0,255,0),tlp.Size(2,2,2)],
        "down": [tlp.Color(255,0,0),tlp.Size(2,2,2)],
        "stable": [tlp.Color(105,105,105),tlp.Size(1,1,1)],
        "intergenic": [tlp.Color(200,200,200),tlp.Size(1,1,1)],
        "nan": [tlp.Color(255,255,255),tlp.Size(1,1,1)]}
    ,"edge" : {
        "gain": [tlp.Color.Blue,10],
        "loss": [tlp.Color.Yellow,10],
        "stable": [tlp.Color.Gray,0]
    }}
    for node in gr.getNodes():
        color[node] = aspect["node"][expression[node]][0]
        size[node] = aspect["node"][expression[node]][1]
    for edge in gr.getEdges():
        viewBorderColor[edge] = aspect["edge"][interaction[edge]][0]
        viewBorderWidth[edge] = aspect["edge"][interaction[edge]][1]
        color[edge] = aspect["edge"][interaction[edge]][0]

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
      if gr.target(edge) in allSubNodes and gr.source(edge) in allSubNodes:
        edges_to_add.append(edge)
    interest.addEdges(edges_to_add)

def set_secondary_regulators(gr,viewLabel):
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
  
  # outputs the results from the dictionary statistics into a CSV file
  file_stats = open("results_statistics.csv","w")
  for (type_of_element, number_of_element) in [(key+"_"+name, nb) for (key,item) in statistics.items() for (name,nb) in item.items()]:
    file_stats.write(type_of_element+"\t"+str(number_of_element)+"\n")
  file_stats.write("genes_TOTAL\t"+str(gr.numberOfNodes())+"\n"+"interactions_TOTAL\t"+str(gr.numberOfEdges())+"\n")
  file_stats.close()
  
  # outputs the results from the dictionary genes_in_pathways into a CSV file
  file_pathways = open("results_genes_in_pathways.csv","w")
  for line in [[key]+values for (key,values) in genes_in_pathways.items()]:
    file_pathways.write("\t".join(line)+"\n")
  file_pathways.close()
  
  return statistics, genes_in_pathways

def get_node_info(node_name,dico_nodes,viewLabel,gr):
  node = dico_nodes[node_name]
  info = {
    "expression": gr["Expression"][node],
    "gain": [],
    "loss": [],
    "stable": []
  }
  for neighbor_edge in gr.getInOutEdges(node):
    neighbor_node = [n for n in gr.ends(neighbor_edge) if n != node][0]
    info[gr["Interaction"][neighbor_edge]].append(viewLabel[neighbor_node]) 
  return info
  
  
def main(gr):
    gr.clear()
    
    viewBorderColor = gr['viewBorderColor']
    viewBorderWidth = gr['viewBorderWidth']
    viewColor = gr['viewColor']
    viewLabel = gr['viewLabel']
    viewSize = gr['viewSize']
    
    ###Functions to create the graph
    dico_nodes=create_interaction_graph(gr,viewLabel)
    add_expression(gr,dico_nodes)
    second_degree_reg = set_secondary_regulators(gr,viewLabel)
    updateVisualization(centerViews = True)
    print("\nGraph constructed successfully")
    
    ###Customization of the nodes and edge regarding their properties
    visu_node_edge(gr,viewSize,viewColor,viewBorderColor,viewBorderWidth)
    print("\nCustomization done")

    ###temporaire : applique automatique FM³
    visu_algoFM(gr)
    
    ###Creation of subgraphs for each pathway
    set_subgraphs_pathways(gr,viewLabel,dico_nodes)
    print("\nSubgraphs of pathways created")
    
    create_interest_subgraph(gr)
    print("\nInterest subgraph created")
    
#    statistics, genes_in_pathways = get_statistics(gr, viewLabel)

    info_SNX9 = get_node_info("SNX9",dico_nodes,viewLabel,gr)
    info_SYNJ2 = get_node_info("SYNJ2",dico_nodes,viewLabel,gr)
    
    for (info,node_name) in [(info_SNX9, "SNX9"),(info_SYNJ2, "SYNJ2")]:
        print("\n>",node_name,":") 
        for (key,values) in info.items():
            if values != [] and key != "expression":
                print(key,":",values)
