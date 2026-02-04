# CBAS
Code to run Choice-Wide Behavioral Association Study

CBAS is written for IgorPro. The version in the repository was run on Igor Pro version 8.

There are two files, one that will run the comparative CBAS, and another that will run the correlatonal CBAS.

To run either version of CBAS, subject data has to be uploaded to data folders. The names of the data folders need to be written in a text file called 'expInfo.' Each data folder needs 4 waves, labeled 'sex,' 'types,' lesion,' 'implant.' Each of the waves need to be the number of subjects long, and can contain qualifiers about the subjects to be selected for.

Data needs to be structured in a specific way. The data for each subject should be in a number of trials x 5 matrix.
The first column contains the session number (0 ... N)
The second column contains the choice of the subject
The third column contains the presence (1) or absence (0) of reward
The fourth and fifth column contain information about the contingency. For the CBAS run in this manuscript there was only a single contingency; however, the code expects there to be at least two contingencies. The first contingency is not analyzed, and column four and five should have NaN. The second contingency is analyzed, and column four should have a 2 in it and column 5 should have a 1 in it. There can be any numbers of session per contingency.

The constants listed at the top of the file are as follows:
NumConting: number of contingencies to be analyzed.
NumArms: the number of different choices available to the subject
SeqLenMax: the maximum number of choices evaluated by CBAS (CBAS will evaluate all sequences in the dataset up to and including this length)
ContingList: a list of strings, separated by ';' that contains the contingency number to be evaluated. The number has to be greater than 0.
ResampleNumber: the total number of resamples that will be done by the Romano-Wolf method for multiple comparisons correction.

To run the comparative CBAS a wave that contains the two groups of subjets to be compared needs to be created. It needs to be an N x 2 matrix, where N is the total number of animals. The first column contains the subject number, and the second column contains a 0 or 1 demarcating which group the animal is in. 

To run the correlative CBAS a wave that contains the subjects for analysis needs to be created. It needs to be an N x 2 matrix, where N is the total number of animals. The first column contains the subject number, and the second column contains the value of the covariate of interest for that subject. 

To run CBAS run: 'runCBAS(anInfoWv)' in the command line of IGOR, where anInfoWv is the matrix described above that is different for the comparative or correlative CBAS.
