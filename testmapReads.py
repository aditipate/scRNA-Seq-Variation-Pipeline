import os

##############################################################################################
### ALL RUNS INCLUDED: testSRRs = ['SRR9290711', 'SRR9290712', 'SRR9290713', 'SRR9290714', 'SRR9290715', 'SRR9290716', 'SRR9290717', 'SRR9290718'] ###
##############################################################################################


#combined cell ranger - ALL
def run_CellRanger_test_ALLRuns():
    current_path = os.getcwd()                                                   #get current working directory
    id  = 'test_combined_cellranger_output'                                      #path designated for cellranger output
    fastqs_path = current_path + '/test_mouse_heart_SRA_data'                    #path to folder containing fastqs from SRA
    transcriptome = current_path + '/mouse_genome/refdata-gex-mm10-2020-A'       #path to mouse genome
    cellranger_cmd = 'cellranger count --id=' + id + ' --fastqs=' + fastqs_path + ' --sample=' + 'SAN,AVN,LPF,RPF' + ' --transcriptome=' + transcriptome
    os.system(cellranger_cmd)

#SAN cell ranger - 'SRR9290711', 'SRR9290712'
def runCellRanger_test_SAN():
    current_path = os.getcwd()                                                   #get current working directory
    id  = 'test_SAN_cellranger_output'                                           #path designated for cellranger output
    fastqs_path = current_path + '/test_mouse_heart_SRA_data'                    #path to folder containing fastqs from SRA
    transcriptome = current_path + '/mouse_genome/refdata-gex-mm10-2020-A'       #path to mouse genome
    cellranger_cmd = 'cellranger count --id=' + id + ' --fastqs=' + fastqs_path + ' --sample=' + 'SAN' + ' --transcriptome=' + transcriptome
    os.system(cellranger_cmd)


#AVN cell ranger - 'SRR9290713', 'SRR9290714'
def runCellRanger_test_AVN():
    current_path = os.getcwd()                                                   #get current working directory
    id  = 'test_AVN_cellranger_output'                                           #path designated for cellranger output
    fastqs_path = current_path + '/test_mouse_heart_SRA_data'                    #path to folder containing fastqs from SRA
    transcriptome = current_path + '/mouse_genome/refdata-gex-mm10-2020-A'       #path to mouse genome
    cellranger_cmd = 'cellranger count --id=' + id + ' --fastqs=' + fastqs_path + ' --sample=' + 'AVN' + ' --transcriptome=' + transcriptome
    os.system(cellranger_cmd)

#LPF cell ranger - 'SRR9290715', 'SRR9290716'
def runCellRanger_test_RPF():
    current_path = os.getcwd()                                                   #get current working directory
    id  = 'test_LPF_cellranger_output'                                           #path designated for cellranger output
    fastqs_path = current_path + '/test_mouse_heart_SRA_data'                    #path to folder containing fastqs from SRA
    transcriptome = current_path + '/mouse_genome/refdata-gex-mm10-2020-A'       #path to mouse genome
    cellranger_cmd = 'cellranger count --id=' + id + ' --fastqs=' + fastqs_path + ' --sample=' + 'LPF' + ' --transcriptome=' + transcriptome
    os.system(cellranger_cmd)

#RPF cell ranger - 'SRR9290717', 'SRR9290718'
def runCellRanger_test_LPF():
    current_path = os.getcwd()                                                   #get current working directory
    id  = 'test_RPF_cellranger_output'                                           #path designated for cellranger output
    fastqs_path = current_path + '/test_mouse_heart_SRA_data'                    #path to folder containing fastqs from SRA
    transcriptome = current_path + '/mouse_genome/refdata-gex-mm10-2020-A'       #path to mouse genome
    cellranger_cmd = 'cellranger count --id=' + id + ' --fastqs=' + fastqs_path + ' --sample=' + 'RPF' + ' --transcriptome=' + transcriptome
    os.system(cellranger_cmd)


