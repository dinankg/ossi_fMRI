% MATLAB code for Bloch Simulation in Oscillating Steady-State Imaging (OSSI)
clc; clear; close all;

% Choose field strength (in Tesla)
fieldStrength = 0.55; % Options: 0.55, 3
trSweep = 1; % 1 for TR sweep, 0 otherwise (tr vs. nc sweep)
spiralin = 0; % 0 otherwise (1 for Spiral in/ 0 for spiral out)
graymatter = 1;
% Define flip angles
flipAngles = 5:0.5:25;                % Flip angles in degrees
TE = 55;

% Define physical constants and parameters
gammaBar = 42570;                   % gamma/2pi in kHz/T
gamma = gammaBar * 2 * pi;          % gamma in radians/s/T
excitationPulseLength = 3.2;        % RF pulse duration


% Adjust T1, T2, T2*, and T2' based on field strength
switch fieldStrength
    case 3
        T1 = 1433.2; % T1 in ms for 3T
        T2 = 92.6;   % T2 in ms for 3T
        T2s1 = 55;   % 
        T2s2 = 1 / (1 / T2s1 - 0.01 / TE); % Adjusted for 1% signal change at TE=30 for 3T
        T2p1 = 1 / (1 / T2s1 - 1 / T2);    % Adjusted T2 prime for 3T
        T2p2 = 1 / (1 / T2s2 - 1 / T2);    % Adjusted T2 prime for 3T

        % Conditional setup for TR sweep
        if trSweep
            trList = 10:50; % Adjusted for 3T
            nc = 10;
            sweepperlist = trList * nc;
            label_x = 'TR';
            sweeplist = trList;
            savingDir = 'figures/3T/TR/';
        else
            ncList = 1:35; % Adjusted for 3T
            tr = 15;
            sweepperlist = ncList * tr;
            label_x = 'nc';
            sweeplist = ncList;
            savingDir = 'figures/3T/nc/';
        end

    case 0.55
        trList = 10:80;
        ncList = 5:35; 

        if graymatter
            T1 = 717;    % T1 in ms for 0.55T
            T2 = 112;    % T2 in ms for 0.55T
            T2s1 = 86;   
        else
            T1 = 493;    % T1 in ms for 0.55T
            T2 = 89;    % T2 in ms for 0.55T
            T2s1 = 72;   
        end
        T2s2 = 1 / (1 / T2s1 - 0.01 / TE); % Adjusted for 1% signal change at TE=30 for 0.55T
        T2p1 = (1 / (1 / T2s1 - 1 / T2)) * (3 / fieldStrength); % Scaled T2 prime for 0.55T
        T2p2 = (1 / (1 / T2s2 - 1 / T2)) * (3 / fieldStrength); % Scaled T2 prime for 0.55T
        %  setup for TR sweep
        if (trSweep && graymatter)
            nc = 6;
            sweepperlist = trList * nc;
            label_x = 'TR';
            sweeplist = trList;
            savingDir = 'figures/055T/TR/graymatter/';
        elseif (trSweep && ~graymatter)
            nc = 6;
            sweepperlist = trList * nc;
            label_x = 'TR';
            sweeplist = trList;
            savingDir = 'figures/055T/TR/whitematter/';
        elseif (~trSweep && graymatter)
            tr = 50;
            sweepperlist = ncList * tr;
            label_x = 'nc';
            sweeplist = ncList;
            savingDir = 'figures/055T/nc/graymatter/';
        elseif (~trSweep && ~graymatter)
            tr = 50;
            sweepperlist = ncList * tr;
            label_x = 'nc';
            sweeplist = ncList;
            savingDir = 'figures/055T/nc/whitematter/';
        end

    otherwise
        error('Invalid field strength selected');
end

% Initialize arrays for storing results
mnsig = zeros(length(sweeplist), length(flipAngles)); % Mean signal strength
mnsigdiff = mnsig;                  % Mean signal difference
mnsigvar = mnsig;                   % Mean signal variation
pdsigdiff = mnsig;                  % Percent signal difference
pdsigvar = mnsig;                   % Percent signal variation
pdminmaxfsig = mnsig;               % Percent min-max signal difference

% Main simulation loop for flip angles
for flipIndex = 1:length(flipAngles)
    flipAngle = flipAngles(flipIndex);
    faRad = flipAngle * pi / 180;
    arf = faRad / gamma / excitationPulseLength;

    % Nested loop for sweep period
    for sweepIndex = 1:length(sweepperlist)
        if trSweep
            tr = sweepperlist(sweepIndex) / nc;
        else
            nc = sweepperlist(sweepIndex) / tr;
        end

        % Time settings
        a = [excitationPulseLength (tr - excitationPulseLength)/2 (tr - excitationPulseLength)/2];
        obtr = ceil(5 * T1 / tr);
        ntr = obtr + 3 * nc;
        dt = repmat(a(:), ntr, 1);
        tt = cumsum(dt); % Cumulative sum to get actual time points

        % Frequency settings
        frrange = 1/(tr); % Frequency range in kHz (example value)
        nomstep = frrange / 2000; % Nominal step size
        stn = round(1 / (nc * tr) / nomstep); % Steps for frequency offset
        df = 1 / (nc * tr) / stn; % Frequency increment
        fff1 = (0:df:frrange / 2); % Frequency array
        freq = [-fff1(end:-1:2) fff1]; % Symmetric frequency array
        bz1 = freq / gammaBar; % Frequency offsets for Bloch simulation
        bz = repmat(bz1, 3 * ntr, 1); % Repeat for each time interval

        % Determine the sweep period for the current iteration
        sweepper = sweepperlist(sweepIndex);

        %         % Calculate additional frequency component for phase progression
        %         addfr = (1:ntr)' * (tr / sweepper) * (1 / (2 * tr));
        %
        %         % Calculate B1 phase progression
        %         b1ph = 2 * pi * (1:ntr)' .* tr .* addfr;

        % Spoiling phase to dephase spins across TRs
        % This helps in destroying residual transverse magnetization
        spoilseed = 360 / (sweepper / tr); % Spoiling increment
        addph = (1:ntr)' * (spoilseed / 360) * 2 * pi; % Additional phase per TR
        b1ph = cumsum(addph); % Cumulative phase for spoiling

        % Construct B1 field with amplitude and phase for each TR
        b1t = arf * exp(1i * b1ph);
        b1tt = [b1t zeros(size(b1t)) zeros(size(b1t))]; % B1 field in three dimensions
        b1ttt = b1tt';
        b1 = b1ttt(:) * ones(size(bz1)); % Replicate B1 field for all frequency offsets

        % Extract real and imaginary parts of B1 field
        bx = real(b1); % B1x component
        by = imag(b1); % B1y component


        % Initial magnetization
        m0 = [zeros(size(bz1)); zeros(size(bz1)); ones(size(bz1))];
        T1obj = T1 * ones(size(bz1));
        T2obj = T2 * ones(size(bz1));

        % Bloch Simulation
        [mx, my, mz] = blochsim4(m0, bx, by, bz, T1obj, T2obj, dt);
        toff = 2-spiralin;  % where to look at the signal 1 is right before the RF, 2 is just after, 3 is at the spin-echo time
        time = tt(toff:3:end);
        mxy = (mx(toff:3:end,:) + 1i*my(toff:3:end,:));
        % Signal processing
        % Initialize variables for signal processing
        noffs = 20;                             % Number of frequency offsets
        msig = zeros(2 * nc, noffs);            % Matrix to store signal for condition 1
        msig2 = zeros(2 * nc, noffs);           % Matrix to store signal for condition 2

        % Loop over frequency offsets
        for offsetIndex = 1:noffs
            offset = (offsetIndex - 1) / tr / nc / noffs; % Calculate current offset

            % Weighting factors for each condition
            f1 = ones([2 * nc, 1]) * cauch(T2p1, freq + offset) ./ sum(cauch(T2p1, freq + offset) * df);
            f2 = ones([2 * nc, 1]) * cauch(T2p2, freq + offset) ./ sum(cauch(T2p2, freq + offset) * df);

            % Calculate weighted average signal for each offset
            msig(:, offsetIndex) = abs(mean(f1 .* mxy(obtr:obtr + 2 * nc - 1, :), 2));
            msig2(:, offsetIndex) = abs(mean(f2 .* mxy(obtr:obtr + 2 * nc - 1, :), 2));
        end

        % Calculate combined signals and their differences
        comsig = sqrt(mean(msig(1:nc, :).^2, 1));   % RMS signal for condition 1
        comsig2 = sqrt(mean(msig2(1:nc, :).^2, 1)); % RMS signal for condition 2
        diffall = comsig2 - comsig;                 % Difference between signals
        diffsig = mean(diffall);                    % Mean difference
        pdsig = diffsig / mean(comsig) * 100;       % Percent difference in signal
        varsig = max(comsig) - min(comsig);         % Signal variation
        pdvar = varsig / mean(comsig) * 100;        % Percent variation

        % Store results in matrices
        mnsig(sweepIndex, flipIndex) = mean(comsig);
        mnsigdiff(sweepIndex, flipIndex) = diffsig;
        mnsigvar(sweepIndex, flipIndex) = varsig;
        pdsigdiff(sweepIndex, flipIndex) = pdsig;
        pdsigvar(sweepIndex, flipIndex) = pdvar;
        pdminmaxfsig(sweepIndex, flipIndex) = (max(diffall) - min(diffall)) / diffsig;


        %         % Update results
        %         mnsig(sweepIndex, flipIndex) = mean(comsig);
        %         mnsigdiff(sweepIndex, flipIndex) = diffsig;
        %         mnsigvar(sweepIndex, flipIndex) = varsig;
        %         pdsigdiff(sweepIndex, flipIndex) = pdsig;
        %         pdsigvar(sweepIndex, flipIndex) = pdvar;
        %         pdminmaxfsig(sweepIndex, flipIndex) = (max(diffall) - min(diffall)) / diffsig;

        %         % Conditional plotting for specific parameters
        %         if (flipAngle == 10) && (nc == 10)
        %             figure(1)
        %             plot((0:noffs-1)/tr/nc/noffs, comsig, (0:noffs-1)/tr/nc/noffs, comsig2);
        %             xlabel('Frequency Offset');
        %             ylabel('Signal Strength');
        %             legend('baseline','activation')
        %             title('Comparative Signal Strength for Flip Angle 10 and nc=10');
        %         end

    end % End of sweep loop
end % End of flip angle loop

%% Visualization of results
% Plotting and saving results
figure(3);
imagesc(sweeplist, flipAngles, mnsig');
colormap(hot);
colorbar;
title('Average Signal Strength');
xlabel(label_x);
ylabel('Flip angle (deg)');
axis square
saveas(gcf, [savingDir 'AverageSignalStrength.png']);
saveas(gcf, [savingDir 'AverageSignalStrength.pdf']);
saveas(gcf, [savingDir 'AverageSignalStrength.fig']);


figure(4);
imagesc(sweeplist, flipAngles, mnsigdiff');
colormap(hot);
colorbar;
title('Average fMRI Signal Difference');
xlabel(label_x);
ylabel('Flip angle (deg)');
axis square
saveas(gcf, [savingDir 'AveragefMRISignalDifference.png']);
saveas(gcf, [savingDir 'AveragefMRISignalDifference.pdf']);
saveas(gcf, [savingDir 'AveragefMRISignalDifference.fig']);




figure(5);
imagesc(sweeplist, flipAngles, pdsigdiff');
colorbar;
title('Percent fMRI Signal Difference');
xlabel(label_x);
ylabel('Flip angle (deg)');
axis square
saveas(gcf, [savingDir 'PercentfMRISignalDifference.png']);
saveas(gcf, [savingDir 'PercentfMRISignalDifference.pdf']);
saveas(gcf, [savingDir 'PercentfMRISignalDifference.fig']);




figure(6);
imagesc(sweeplist, flipAngles, mnsigvar');
colormap(hot);
colorbar;
title('Average Signal Variation');
xlabel(label_x);
ylabel('Flip angle (deg)');
axis square
saveas(gcf, [savingDir 'AverageSignalVariation.png']);
saveas(gcf, [savingDir 'AverageSignalVariation.pdf']);
saveas(gcf, [savingDir 'AverageSignalVariation.fig']);




figure(7);
imagesc(sweeplist, flipAngles, pdsigvar');
colorbar;
title('Percent Signal Variation');
xlabel(label_x);
ylabel('Flip angle (deg)');
axis square
saveas(gcf, [savingDir 'PercentSignalVariation.png']);
saveas(gcf, [savingDir 'PercentSignalVariation.pdf']);
saveas(gcf, [savingDir 'PercentSignalVariation.fig']);



% Conditional plot for TR sweep
if (trSweep == 1)
    normFactor = ones([length(flipAngles), 1]) * (sqrt(trList) ./ sqrt(trList - 5));
    figure(8);
    imagesc(sweeplist, flipAngles, mnsigdiff' ./ normFactor);
    colormap(hot);
    colorbar;
    title('Average fMRI Signal Difference (Normalized)');
    xlabel(label_x);
    ylabel('Flip angle (deg)');
    axis square
    saveas(gcf, [savingDir 'AveragefMRISignalDifferenceNormalized.png']);
    saveas(gcf, [savingDir 'AveragefMRISignalDifferenceNormalized.pdf']);
    saveas(gcf, [savingDir 'AveragefMRISignalDifferenceNormalized.fig']);


end

