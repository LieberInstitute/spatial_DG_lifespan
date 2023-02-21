%Load the output of refineVNS
load('*mask.mat')
%mask{n} == name of refineVNS output variable name in MATLAB workspace
BW = bwareafilt(mask_dark_blue,[100 2000]); % 'n' will the cluster number associated with nuclei from 2nd run
%Save as .mat file
save('*nuclei_final.mat','BW')