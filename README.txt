Step 0 : Get on the CMSSW_12_3_0_pre4

Step 1 : Get the GEN-SIM for run 3 : cmsRun EXO-HSCP-GENSIM.py

Step 2 : DIGI2RAW : You can modify the file to change input files, number of events processed, what collections you want to keep etc.. 

	 - cmsRun DIGI2RAW.py

Step 3 : Run the custom HLT menu
	- To obtain the menu : git clone https://github.com/DenkMybu/DIGI2NTUPLE.git
	- hltMenu.py is the custom menu, which keeps all the collections for now. It was generated using the following command :
 
	hltGetConfiguration /users/sharper/2021/hscp/CMSSW_1210_GRun --mc --globaltag auto:phase1_2021_realistic --output minimal --eras Run3 > hltMenu.py

	- To obtain the L1menu.root : 

	cmsRun Analysis/HLTAnalyserPy/test/dumpConditions_cfg.py outputFile=l1MenuMC.root globalTag=121X_mcRun3_2021_realistic_v15 startRunNr=1 endRunNr=1

	- To run the custom menu : 
	
	cmsRun hltMenu.py inputFiles=file:HSCP_Gluino_Mass1800_RAW.root outputFile=hltOutput.root >& hlt.log &

Step 4 : From HLT to AOD

	- cmsRun HLT2AOD.py (Change the file with the corresponding names)

Step 5 : From AOD to HSCP production
	- Visit https://github.com/enibigir/SUSYBSMAnalysis-HSCP/tree/dev




----------------------------------------------------------------------------------------------------------------------

This environment is the fork from git HSCP Analysis, commit from here all changes.







