
ME_clusters = ['NM1_7', 'NM1_5', 'NM1_4', 'NM1_14', 'NM1_13', 'NM1_6', 'NM1_10', 'NM1_9', 'NM1_3', 'NM1_8', 'NM1_11', 'NM1_2', 'N1_5', 'N1_14', 'N1_6', 'N1_10', 'N1_9', 'N1_3', 'N1_8', 'N1_11', 'N1_2', 'NM2_4', 'NM2_14', 'NM2_6', 'NM2_10', 'NM2_9', 'NM2_3', 'NM2_8', 'NM2_11', 'NM2_2', 'N2_4', 'N2_13', 'N2_6', 'N2_10', 'N2_9', 'N2_3', 'N2_8', 'N2_11', 'N2_2', 'NM3_4', 'NM3_13', 'NM3_6', 'NM3_10', 'NM3_9', 'NM3_3', 'NM3_8', 'NM3_11', 'NM3_2', 'N3_6', 'N3_10', 'N3_9', 'N3_3', 'N3_8', 'N3_11', 'N3_2', 'N4_6', 'N4_10', 'N4_9', 'N4_3', 'N4_8', 'N4_11', 'N4_2', 'N5_4', 'N5_13', 'N5_6', 'N5_10', 'N5_9', 'N5_3', 'N5_8', 'N5_11', 'N5_2', 'NN1_13', 'NN1_10', 'NN1_9', 'NN1_3', 'NN1_8', 'NN1_11', 'NN1_2']

path = './New_report/Whippet_Delta/'


with  open(path + "TOTAL.ME.diff" , "w") as out_file :
    
    TOTAL_DIFF  = csv.writer(out_file, delimiter="\t") 
    
    header = ["Gene", "Node", "Coord", "Strand", "Type", "Psi_A", "Psi_B", "DeltaPsi", "Probability" ,"Complexity", "Entropy"]

    
    
    TOTAL_DIFF.writerow( ["ME_cluster", "Tissue_cluster"] + header + ["Potential_Exon", "Is_Annotated"])
    
    for c in ME_clusters:


        ME_cluster, Tissue_cluster =  c.split("_")

        with open(path + c+ ".ME.diff" ) as F: 

            reader = csv.DictReader(F, delimiter="\t")

            for row in reader:
                
                out =  [ME_cluster, Tissue_cluster] + [row[x] for x in header]
                
                #print( node_exons[ row["Coord"] + ":" + row["Strand"] ])
                
                if (row["Gene"], row["Node"] ) in node_exons:
                    
                    out =  [ME_cluster, Tissue_cluster] + [row[x] for x in header] + node_exons[(row["Gene"], row["Node"] )]
                
                    TOTAL_DIFF.writerow(out)
