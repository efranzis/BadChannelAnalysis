# BadChannelAnalysis

Here are some basic macros that are used for the
bad channel analysis.

First 
Put some AnalysisResults.root files that you have
copied e.g with the runCopyFromLegoTrain.sh script to your local
directory —>AnalysisInput.
Place a .txt file in the AnalysisInput folder that contains the run numbers.
It can ideally be called runList.txt

Second open root
run the BadChannelAnalysis.C
The convert function will use the AnalyisInput files and 
store a merged file in the ConvertOutput folder.

Third 
the bad channel analysis will be performed on this merged file
