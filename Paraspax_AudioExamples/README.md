Supplementary material (audio examples) for:   
Arend, J. M., Amengual Garí, S. V., Schissler, C., Klein, F., & Robinson, P. W. (2021). Six-Degrees-of-Freedom Parametric Spatial Audio Based on One Monaural Room Impulse Response. *J. Audio Eng. Soc.*, 69(7/8), 557–575. https://doi.org/10.17743/jaes.2021.0009

The folder `BinauralRenderings` contains binaural renderings based on binaural room impulse responses (BRIRs) synthesized with the Paraspax method as described in the manuscript. The BRIRs were synthesized based on one monaural room impulse response (RIR) measured for each sound source at position 10 of the grid. The spatialization was done based on a 2nd order image source model, as described in the manuscript. 

The audio files are recordings of the PyBinSim output, i.e., the binaural output of the 6-DoF VAE presented in the manuscript. All examples are rendered for a static head orientation (frontal listener-related head orientation in the local coordinate system). Renderings are provided for 4 different source positions, 20 different receiver positions, and various audio content. The naming is as follows:  

`Src_1_2_3_4_Pos_2_PinkNoise`: Examples with pink noise, position 2 on the grid, sources 1 to 4 played sequentially.  
`Src_1_Pos_19_15_11_7_3_Speech`: Examples with speech, source 1, positions 19,15,11,7,3 played sequentially (movement towards the sound source in steps of 1 m).

Please note that the provided binaural renderings are without headphone compensation. We provide various headphone compensation filters (minimum phase FIR filters in `*.wav` format) for the employed KEMAR HRTFs in the folder `HPCF`. If no filter is available for the headphones used, we also provide the diffuse field filter (`DFC_KEMAR.wav`; inverse of the common transfer function) of the employed HRTF set as another option for equalization.

