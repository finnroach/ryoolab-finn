# ryoolab-finn

The code in this repository scans regulatory DNA sequences to determine potential transcription factors binding sites. Transcription factors have DNA binding domains which will bind to an optimal squence, and also to sequences which vary at specific positions in the binding site from that optimal sequence. The position frequency or weight matrix of the transcription factor quantifies this variation. These matrices allow the code to predict binding affinity by producing a "binding score" for any given sequence of DNA, with a higher binding score corresponding to higher putative affinity. The top binding scores – those above a certain percentage of the optimal binding score - are presented by the code in a graph, displayed within the DNA sequence, and printed alongside their percentage of the maximum binding score. 

The py file "gstd1_thor.py" contains a pared-down code that will run a specific file of genes for a specific transcription factor. The py file "template_code.py" contains a more flexible version of the code which allows the user to input any csv file with sequences, any transcsription factor matrix, and any cutoff without changing the script. The csv file "gstd1_thor.csv" is the spreadsheet imported by "gstd1_thor.py" and contains the sequences which the code scans. The pdf "gstd1_thor_output.pdf" contains the output generated when "gstd1_thor.py" is run. 
