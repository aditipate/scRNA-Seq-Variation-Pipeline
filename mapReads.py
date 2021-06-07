import os

##############################################################################################
### ONE RUNS INCLUDED: testSRRs = ['SRR9290711', 'SRR9290713', 'SRR9290715', 'SRR9290717'] ###
##############################################################################################

#combined cell ranger - ONE run per zone
def runCellRanger_SNGRuns():
    current_path = os.getcwd()                                                   #get current working directory
    id  = 'combined_cellranger_output'                                           #path designated for cellranger output
    fastqs_path = current_path + '/mouse_heart_SRA_data'                         #path to folder containing fastqs from SRA
    transcriptome = current_path + '/mouse_genome/refdata-gex-mm10-2020-A'       #path to mouse genome
    cellranger_cmd = 'cellranger count --id=' + id + ' --fastqs=' + fastqs_path + ' --sample=' + 'SAN,AVN,LPF,RPF' + ' --transcriptome=' + transcriptome
    os.system(cellranger_cmd)

#SAN cell ranger - 'SRR9290711'
def runCellRanger_SAN():
    current_path = os.getcwd()                                                   #get current working directory
    id  = 'SAN_cellranger_output'                                                #path designated for cellranger output
    fastqs_path = current_path + '/mouse_heart_SRA_data'                         #path to folder containing fastqs from SRA
    transcriptome = current_path + '/mouse_genome/refdata-gex-mm10-2020-A'       #path to mouse genome
    cellranger_cmd = 'cellranger count --id=' + id + ' --fastqs=' + fastqs_path + ' --sample=' + 'SAN' + ' --transcriptome=' + transcriptome
    os.system(cellranger_cmd)


#AVN cell ranger - 'SRR9290713'
def runCellRanger_AVN():
    current_path = os.getcwd()                                                   #get current working directory
    id  = 'AVN_cellranger_output'                                                #path designated for cellranger output
    fastqs_path = current_path + '/mouse_heart_SRA_data'                         #path to folder containing fastqs from SRA
    transcriptome = current_path + '/mouse_genome/refdata-gex-mm10-2020-A'       #path to mouse genome
    cellranger_cmd = 'cellranger count --id=' + id + ' --fastqs=' + fastqs_path + ' --sample=' + 'AVN' + ' --transcriptome=' + transcriptome
    os.system(cellranger_cmd)

#LPF cell ranger - 'SRR9290715'
def runCellRanger_LPF():
    current_path = os.getcwd()                                                   #get current working directory
    id  = 'LPF_cellranger_output'                                                #path designated for cellranger output
    fastqs_path = current_path + '/mouse_heart_SRA_data'                         #path to folder containing fastqs from SRA
    transcriptome = current_path + '/mouse_genome/refdata-gex-mm10-2020-A'       #path to mouse genome
    cellranger_cmd = 'cellranger count --id=' + id + ' --fastqs=' + fastqs_path + ' --sample=' + 'LPF' + ' --transcriptome=' + transcriptome
    os.system(cellranger_cmd)

#RPF cell ranger - 'SRR9290717'
def runCellRanger_RPF():
    current_path = os.getcwd()                                                   #get current working directory
    id  = 'RPF_cellranger_output'                                                #path designated for cellranger output
    fastqs_path = current_path + '/mouse_heart_SRA_data'                         #path to folder containing fastqs from SRA
    transcriptome = current_path + '/mouse_genome/refdata-gex-mm10-2020-A'       #path to mouse genome
    cellranger_cmd = 'cellranger count --id=' + id + ' --fastqs=' + fastqs_path + ' --sample=' + 'RPF' + ' --transcriptome=' + transcriptome
    os.system(cellranger_cmd)

