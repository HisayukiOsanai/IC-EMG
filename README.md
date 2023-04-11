# IC-EMG
Extract EMG signal from skull-screw-referenced LFP data using ICA
Tha manuscript can be found here: https://www.biorxiv.org/content/10.1101/2022.11.15.516633v2.abstract


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
  
Examples:    
Raw LFP (skull-screw referenced):  
<img src="https://user-images.githubusercontent.com/60276754/231026752-40bd2b91-3727-4023-8b85-a0f9de476b3e.png"  width="300" height="500">

Independent components (red label indicates IC-EMG):  
<img src="https://user-images.githubusercontent.com/60276754/231026876-833dbd15-f43e-4453-83a8-895f3ea737ee.png"  width="300" height="500">

Weight distributions of ICs:  
<img src="https://user-images.githubusercontent.com/60276754/231026943-6198d1c6-d3cb-4778-bc6c-f4c37413c8c7.png"  width="300" height="500">
<img src="https://user-images.githubusercontent.com/60276754/231026971-1385883f-40e6-408e-afc8-36f15a46a8c4.png"  width="200" height="150">  
most-uniform (bottom-right) IC is identified as IC-EMG.

  
IC-EMG and real-EMGs:  
<img src="https://user-images.githubusercontent.com/60276754/231027233-59759fef-26f0-443f-be62-6539ae686964.png"  width="300" height="200">
