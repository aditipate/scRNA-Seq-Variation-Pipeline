#set up environment
import testretrieveDataSRA
import testreformatDataSRA
import testmapReads

#GEO Accession  -- CCS Component -- SRA run (raw data)
#GSM3885058     -- Zone I: SAN region -- SRR9290711/2
#GSM3885059	-- Zone II: AVN/His region -- SRR9290713/4
#GSM3885060	-- Zone III (Left): BB/PF region (Left Ventricle) -- SRR9290715/6
#GSM3885061	-- Zone III (Right): BB/PF region (Right Ventricle) -- SRR9290717/8

#list of SRA sample runs (one of each)
testSRRs = ['SRR9290711', 'SRR9290712', 'SRR9290713', 'SRR9290714', 'SRR9290715', 'SRR9290716', 'SRR9290717', 'SRR9290718']

#retrieve mouse heart data from NCBI's SRA + retrieve mouse reference genome
testretrieveDataSRA.getSRAdata(testSRRs)
testretrieveDataSRA.getRefGenome()

#rename fastq files to Cell Ranger compatabile format
testreformatDataSRA.renameFastqs(testSRRs)

#map sample reads using Cell Ranger and create Cell Ranger output folder
testmapReads.run_CellRanger_test_ALLRuns()
testmapReads.runCellRanger_test_SAN()
testmapReads.runCellRanger_test_AVN()
testmapReads.runCellRanger_test_LPF()
testmapReads.runCellRanger_test_RPF()
