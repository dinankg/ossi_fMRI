function seq_params = getparams(seq_params)
%Baseline values for Low-field scanner.
seq_params.maxgrad = 15.5;%mT/m
seq_params.maxslew = 40;%T/m/s
seq_params.rfDtime = 100e-6;
seq_params.rfRtime = 60e-6;
seq_params.adcDtime = 40e-6;
seq_params.gradRasterTime = 10e-6;
seq_params.adcRasterTime = 4e-6;
% Acquisition parameters
seq_params.fov = [240e-3,240e-3,30e-3]; % FOV
seq_params.N = [128,128,10];%resolution
seq_params.slicethickness = seq_params.fov(3) ;
seq_params.alpha = 15; % flip angle (degrees)
% params for fieldmap
seq_params.deltak = 1./seq_params.fov;
seq_params.nCyclesSpoil = 2;
seq_params.rfSpoilingInc = 117;            % RF spoiling increment
seq_params.fatChemShift = 3.5e-6;          % 3.5 ppm
gamma_pulseq = 42576000;
seq_params.B0= .55;
seq_params.gamma    = 4.2576e3;       % Hz/Gauss
% seq_params.fatOffresFreq = gamma_pulseq*B0_pulseq*seq_params.fatChemShift;  % Hz
% seq_params.TE = 1/seq_params.fatOffresFreq*[1 2];     % fat and water in phase for both echoes
end