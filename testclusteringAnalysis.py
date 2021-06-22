import os

#run Clustering.R to normalize and cluster data

def performClusteringFiltered():
    current_path = os.getcwd()                                  #create seurat results folder
    results_folder = "/seurat_output"
    os.mkdir(current_path + results_folder)


    seurat_cmd = 'Rscript testClusteringFiltered.R'
    os.system(seurat_cmd)

    print("Pipeline is done running! Open 'seurat_ouput' folder to view significant results!")
    print(' _  _')
    print('(o)(o)--. ')
    print(' \../ (  )')
    print(' m\/m--m--')


def performClusteringUnfiltered():

    seurat_cmd = 'Rscript testClusteringUnFiltered.R'
    os.system(seurat_cmd)

    print("Pipeline is done running! Open 'seurat_ouput' folder to view significant results!")
    print(' _  _')
    print('(o)(o)--. ')
    print(' \../ (  )')
    print(' m\/m--m--')