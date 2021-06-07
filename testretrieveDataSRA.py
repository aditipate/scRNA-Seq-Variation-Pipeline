import os

#removed '-I' parameter, and downloading additional runs



#function to retrieve mouse heart data from NCBI's SRA database
def getSRAdata(SRRs):
    current_path = os.getcwd()                                                   #create and change to SRA data folder directory
    folder = "/test_mouse_heart_SRA_data"
    os.mkdir(current_path + folder)
    os.chdir(current_path + folder)

    for SRR in SRRs:
        if(SRR[-1] == '2' or SRR[-1] == '3' or SRR[-1] == '8'):
            SRR_data_address = 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-15/' + SRR + '/' + SRR + '.1'
            wget_SRR = 'wget' + ' ' + SRR_data_address
            fastq_dump_SRR = 'fastq-dump --split-files' + ' ' + SRR + '.1'
            os.system(wget_SRR)                                                #retrieve SRA data using wget command
            os.system(fastq_dump_SRR)                                          #uncompress data and convert to paired-end fastq files using fastq-dump command

        elif(SRR[-1] == '4'):
            SRR_data_address = 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/' + SRR + '/' + SRR + '.1'
            wget_SRR = 'wget' + ' ' + SRR_data_address
            fastq_dump_SRR = 'fastq-dump --split-files' + ' ' + SRR + '.1'
            os.system(wget_SRR)                                                #retrieve SRA data using wget command
            os.system(fastq_dump_SRR)                                          #uncompress data and convert to paired-end fastq files using fastq-dump command

        elif (SRR[-1] == '1' or SRR[-1] == '5' or SRR[-1] == '6' or SRR[-1] == '7' ):
            SRR_data_address = 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/' + SRR + '/' + SRR + '.1'
            wget_SRR = 'wget' + ' ' + SRR_data_address
            fastq_dump_SRR = 'fastq-dump --split-files' + ' ' + SRR + '.1'
            os.system(wget_SRR)                                                #retrieve SRA data using wget command
            os.system(fastq_dump_SRR)                                          #uncompress data and convert to paired-end fastq files using fastq-dump command

    os.chdir(current_path)                                                     #change to current directory
