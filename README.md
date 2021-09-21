SordilloBargmann2021

This repository contains data files pertaining to Sordillo, A. and Bargmann, C. I. (2021). 

Questions about this repository can be directed towards Aylesse Sordillo at aylesse.sordillo@gmail.com or Cori Bargmann at cori.rockefeller.edu 

These data were collected in the Bargmann Lab at Rockefeller University in New York City from 2017-2021.

Funding was provided by HHMI and Chan Zuckerberg Science Initiative.

1) Data sets relevant to all of distributions in the main figures (excluding the optogenetics experiments) are included as both a .mat file that contains a MATLAB structure and a .csv file. Data sets from RIC supplements are also included for those interested.

Each genotype (or transgenic line) tested has its own files.  Files are named by experiment and then genotype/line.

All data provided are from 4-8 min off food during local search.
Corresponding rows for reversal and forward run parameters are from the same event.

Analysis was completed using custom MATLAB functions included here.
Nested MATLAB structures were converted to .csv files using the following:
Gero Nootz (2021). Nested Structure to table and or text file (https://www.mathworks.com/matlabcentral/fileexchange/51271-nested-structure-to-table-and-or-text-file). MATLAB Central File Exchange. Retrieved January 17, 2020.


Data Collection: 
Behavioral states were extracted from the State array generated by BargmannWormTracker. The tracker can be accessed at: https://github.com/navinpokala/BargmannWormTracker

Local search event frequencies per minute were calculated 4-8 min after removal from food. Global search frequencies per minute were calculated 36-40 min after removal from food. Only tracks that were continuous for the entire four-minute time interval were included in frequency analysis. When calculating frequencies, tracks taken on a single day from a single assay plate were averaged to give a single data point, e.g. in Figure 2B and 2D.

Distributions of reversal parameters and forward run durations were calculated using events observed during local search, 4-8 minutes after removal from food. All reversals were included; only forward runs over 2 s in length were included. Reversal length is the path length calculated using the X-Y coordinates, worm length, and pixel size extracted from the tracker. Reversal and forward run speed are the average of mean and median speed extracted from the tracker.

Tracks that were less than 5 minutes long, tracks approaching the copper barrier, and tracks that did not include a complete reversal or forward run were not included in reversal and forward run parameter analyses.
Custom MATLAB functions used to generate these data sets are included here (SordilloCode1.m, SordilloCode2.m).


Description of Variables:
	
	Genotype: A label describing the genotype or strain in the dataset

	Experiment: A label describing the experiment

		AllReversals: includes all reversals (Reversal Omegas, Pure Reversals, other)
	
		ReversalOmegas: reversals coupled to sharp omega turns

		PureReversals: reversals followed by immediate forward 	movement

		ForwardRuns: forward crawls (as opposed to reversals) Forward Runs with < 2 s duration are not included.

			Length: reversal length (body lengths)

			Speed: reversal or forward run speed (mm/s)

			Duration: reversal or forward run duration (s)

		Events per animal per minute: each value represents the count of a behavioral event per animal divided by the number of minutes (4)

			AllReversals: includes all reversals (Reversal Omegas, Pure Reversals, other)

			ReversalOmegas: reversals coupled to sharp omega turns
			
			PureReversals: reversals followed by immediate forward movement

			LongReversals: reversals > 0.5 body lengths

			ShortReversals: reversals < 0.5 body lengths

			ForwardRuns: forward crawls (as opposed to reversals, not filtered by duration).

			Pauses: animal is not moving forward or backward (based on low velocity, see Bargmann Lab Tracker).


2) A source data file (Sordilo_SourceData.xlxs) includes data pertaining to all dot plots is also included and labeled by figure panel. 

For behavioral data, data from individual tracks recorded on a single day from a single assay plate were averaged to give a single data point, e.g. in Figure 2B and 2D. Genotypes are as described in the figures/text.

For imaging data each dot is data from a single trace (Figure 5-figure supplement3). 
	
3) Imaging data (Figure 5-figure supplement3) was tracked using an ImageJ Macro originally used in in Larsch et al., 2013 (https://doi.org/10.1073/pnas.1318325110). Tracks were cleaned and processed using custom MATLAB and Python functions included here (SordilloCode3.m, SordilloCode4).



