# CRLB_T2MSE_EMCopt
CRLB_T2MSE_EMCopt it is a toolbox for sequence optimization with CRLB for T2 Multi Spin-Echo w/ Echo Modulating Curves &amp; Extended Phase Graphs - Open-Source Toolbox

With this toolbox, Multi Spin-Echo sequences (MSE) for T2 cartilage mapping with a dictionary-based approach were optimised comparing 
CRLB values. Candidate combinations of acquisition parameters can be identified for specific target T2 values
and validated through simulation.

Cramer-Rao Lower Bound values of Multi Spin-Echo sequences were compared to optimize dictionary-based T2 cartilage estimation. The refocusing flip angle,
inter-echo spacing and number of echoes were varied, while constraining the B1+rms to a maximum of 5μT and the exam time not to exceed 9 minutes,
keeping TR to a minimum.The effect of the acquisition time can be accounted for in the comparisons. Two different T2 values (8 and 45 ms), can be
representative of the cartilage is considered, and the optimization results can be tested through simulation, confirming improved T2 precision
for the sequence optimised for each target T2 value.

The present implementation aims to contribute towards the development of open-source tools for MRI, particularly of sequences based on MSE for cartilage imaging

Disclaimer: we use the following toolboxes:

   	https://www.jemris.org/

    	https://github.com/pulseq/pulseq

		

To run the matlab codes, download and add the following toolboxes to the Tools folder:

    JEMRIS - https://www.jemris.org/
    Pulseq - https://github.com/pulseq/pulseq
