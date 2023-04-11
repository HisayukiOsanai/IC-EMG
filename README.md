# IC-EMG
extract EMG signal from skull-screw-referenced LFP data using ICA

example data and results are placed in \example, and \results.


main_runica_230410_EMG.m:

This is an example code that shows how to extract EMG signals from
skull-screw-referenced LFP data using independent component analysis (ICA). 

Input:
  ICEMG_threshold: threshold of weight distribution uniformity for
                   detecting IC-EMG. 
                   default: ICEMG_threshold = 0.1;
  filename: path to LFP file (causion: ICA computaion may take long time for heavy input files)
  chn:    number of LFP channels
  Fs:     sampling rate [Hz]
  badch:  abnormal ch# you want to avoid in computing ICA. Including
          abnormal channels may cause unwanted results.
          example:    badch = []; if all channels are healthy
                      badch = [15, 27]; if you want to avoid ch#15 and #27
  
  (optional) EMGfilename: path to EMG file


Output:
  EEG.icaact = ICs obtained by ICA from LFP
  EEG.icaweights, EEG.icasphere: both are relating to IC weight distribution.   
      LFP and IC have relationship that  LFP = mixing*IC,
        where  mixing = pinv(EEG.icaweights * EEG.icasphere).
      Weight distribution of each IC is shown in each column of "mixing" matrix.          
  
  ICEMG:  IC-EMG, identified with uniformity of IC weight distribution.
  ICEMGi: Index of IC-EMG in ICs.
