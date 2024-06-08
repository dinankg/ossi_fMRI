% MATLAB code for Bloch Simulation in Balanced Steady-State Free Precession Imaging (bSSFP)
clc; clear; close all;

% Choose field strength (in Tesla)
fieldStrength = 0.55; % Options: 3, 0.55

% TR and NC List
trList = 10:80; % TR list in ms
ncList = 5:35; % Number of cycles list
TE = 55;
graymatter = 1;

% Define physical constants
gammaBar = 42570; % gamma/2pi in kHz/T
gamma = gammaBar * 2 * pi; % gamma in radians/s/T
excitationPulseLength = 3.2; % RF pulse duration in ms

% Adjust T1, T2, T2*, T2', flipAngleDegrees, and savingDir based on field strength
switch fieldStrength
    case 3
        T1 = 1433.2; % T1 in ms for 3T
        T2 = 92.6;   % T2 in ms for 3T
        T2s1 = 55;   % Example T2* value for 3T
        T2s2 = 1 / (1 / T2s1 - 0.01 / TE); % Adjusted for 1% signal change at TE=30 for 3T
        T2p1 = 1 / (1 / T2s1 - 1 / T2);    % T2 prime for 3T
        T2p2 = 1 / (1 / T2s2 - 1 / T2);    % Adjusted T2 prime for 3T
        flipAngleDegrees = 10;             % Flip angle in degrees for 3T
        savingDir = 'figures/3T/FA/';

    case 0.55
        flipAngleDegrees = 20;             % Flip angle in degrees for 0.55
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
        if graymatter
            savingDir = 'figures/055T/FA/graymatter/';
        else
            savingDir = 'figures/055T/FA/whitematter/';
        end

    otherwise
        error('Invalid field strength selected');
end

% Convert flip angle to radians
flipAngleRadians = flipAngleDegrees * pi / 180; % Flip angle in radians
arf = flipAngleRadians / gamma / excitationPulseLength; % Amplitude of RF field

%initialize outputs
mnsig = zeros(length(trList),length(ncList));
mnsigdiff = zeros(length(trList),length(ncList));
mnsigvar = zeros(length(trList),length(ncList));
pdsigdiff = zeros(length(trList),length(ncList));
pdsigvar = zeros(length(trList),length(ncList));

% Main simulation loop
for lpnc = 1:length(ncList)
    nc = ncList(lpnc); % Current number of cycles

    for lpsw = 1:length(trList)
        tr = trList(lpsw); % Current TR

        % Time settings
        timeIntervals = [excitationPulseLength (tr - excitationPulseLength)/2 (tr - excitationPulseLength)/2];
        obtr = ceil(5 * T1 / tr);
        ntr = obtr + 3 * nc;
        dt = repmat(timeIntervals(:), ntr, 1);
        tt = cumsum(dt); % Cumulative sum to get actual time points


        frrange = 1/tr; % Frequency range in kHz (example value)
        nomstep = frrange / 2000; % Nominal step size
        stn = round(1 / (nc * tr) / nomstep); % Steps for frequency offset
        df = 1 / (nc * tr) / stn; % Frequency increment
        fff1 = (0:df:frrange / 2); % Frequency array
        freq = [-fff1(end:-1:2) fff1]; % Symmetric frequency array
        bz1 = freq / gammaBar; % Frequency offsets for Bloch simulation
        bz = repmat(bz1, 3 * ntr, 1); % Repeat for each time interval


        % Calculate B1 field
        sweeper = tr*nc;
        spoilSeed = 360 / (sweeper/tr);
        addPhase = (1:ntr)' * (spoilSeed / 360) * 2 * pi;
        b1Phase = cumsum(addPhase);
        b1t = arf * exp(1i * b1Phase);
        b1Field = [b1t zeros(size(b1t)) zeros(size(b1t))];
        b1FieldTranspose = b1Field';
        b1 = b1FieldTranspose(:) * ones(size(bz1));
        bx = real(b1);
        by = imag(b1);


        % Initial magnetization
        m0 = [zeros(size(bz1)); zeros(size(bz1)); ones(size(bz1))];
        T1obj = T1 * ones(size(bz1));
        T2obj = T2 * ones(size(bz1));

        % Bloch Simulation
        [mx, my, mz] = blochsim4(m0, bx, by, bz, T1obj, T2obj, dt);
        toff = 2;  % Time offset for observing the signal
        time = tt(toff:3:end);
        mxy = (mx(toff:3:end,:) + 1i * my(toff:3:end,:));

        %Cauchy weighing
        noffs = 20; % Number of frequency offsets
        msig = zeros(2 * nc, noffs);
        msig2 = zeros(2 * nc, noffs);
        for lp = 1:noffs
            offset = (lp - 1) / tr / nc / noffs;
            f1 = ones([2 * nc, 1]) * cauch(T2p1, freq + offset) ./ sum(cauch(T2p1, freq + offset) * df);
            f2 = ones([2 * nc, 1]) * cauch(T2p2, freq + offset) ./ sum(cauch(T2p2, freq + offset) * df);
            msig(:, lp) = abs(mean(f1 .* mxy(obtr:obtr + 2 * nc - 1, :), 2));
            msig2(:, lp) = abs(mean(f2 .* mxy(obtr:obtr + 2 * nc - 1, :), 2));
        end

        % Calculate combined signals and their differences
        comsig = sqrt(mean(msig(1:nc, :).^2, 1));
        comsig2 = sqrt(mean(msig2(1:nc, :).^2, 1));
        diffsig = mean(comsig2 - comsig);
        pdsig = diffsig / mean(comsig) * 100;
        varsig = max(comsig) - min(comsig);
        pdvar = varsig / mean(comsig) * 100;

        % Store results
        mnsig(lpsw, lpnc) = mean(comsig);
        mnsigdiff(lpsw, lpnc) = diffsig;
        mnsigvar(lpsw, lpnc) = varsig;
        pdsigdiff(lpsw, lpnc) = pdsig;
        pdsigvar(lpsw, lpnc) = pdvar;
    end
end

%% Plotting and Visualization

% Define colormap
colormapType = hot;

% Average Signal Strength
figure;
imagesc(trList, ncList, mnsig');
colormap(colormapType);
colorbar;
title('Average Signal Strength');
xlabel('TR (ms)');
ylabel('nc');
axis square;
saveas(gcf, [savingDir 'AverageSignalStrength.png']);
saveas(gcf, [savingDir 'AverageSignalStrength.pdf']);
saveas(gcf, [savingDir 'AverageSignalStrength.fig']);


% Average fMRI Signal Difference
figure;
imagesc(trList, ncList, mnsigdiff');
colormap(colormapType);
colorbar;
title('Average fMRI Signal Difference');
xlabel('TR (ms)');
ylabel('nc');
axis square;
saveas(gcf, [savingDir 'AveragefMRISignalDifference.png']);
saveas(gcf, [savingDir 'AveragefMRISignalDifference.pdf']);
saveas(gcf, [savingDir 'AveragefMRISignalDifference.fig']);


% Percent fMRI Signal Difference
figure;
imagesc(trList, ncList, pdsigdiff');
colormap(colormapType);
colorbar;
title('Percent fMRI Signal Difference');
xlabel('TR (ms)');
ylabel('nc');
axis square;
saveas(gcf, [savingDir 'PercentfMRISignalDifference.png']);
saveas(gcf, [savingDir 'PercentfMRISignalDifference.pdf']);
saveas(gcf, [savingDir 'PercentfMRISignalDifference.fig']);



% Average Signal Variation
figure;
imagesc(trList, ncList, mnsigvar');
colormap(colormapType);
colorbar;
title('Average Signal Variation');
xlabel('TR (ms)');
ylabel('nc');
axis square;
saveas(gcf, [savingDir 'AverageSignalVariation.png']);
saveas(gcf, [savingDir 'AverageSignalVariation.pdf']);
saveas(gcf, [savingDir 'AverageSignalVariation.fig']);



% Percent Signal Variation
figure;
imagesc(trList, ncList, pdsigvar');
% colormap(colormapType);
colorbar;
title('Percent Signal Variation');
xlabel('TR (ms)');
ylabel('nc');
axis square;
saveas(gcf, [savingDir 'PercentSignalVariation.png']);
saveas(gcf, [savingDir 'PercentSignalVariation.pdf']);
saveas(gcf, [savingDir 'PercentSignalVariation.fig']);
%%

% normalize signal difference
figure;
normFactor = ones([length(ncList), 1]) * (sqrt(trList) ./ sqrt(trList - 5));
data_sig = mnsigdiff' ./ normFactor;
[maxVal,linearInd] = max(data_sig(:));
[row,col] = ind2sub(size(data_sig), linearInd);
imagesc(trList, ncList,data_sig );
colormap(colormapType);
colorbar;
hold on 
% Plot the marker at the location of the maximum value
plot(trList(col), ncList(row), 'rs', 'MarkerSize', 5, 'LineWidth', 2);hold off 

% Display the coordinates at the location of the maximum signal
text(trList(col), ncList(row), sprintf('nc = %d, TR = %d', round(ncList(row)), round(trList(col))), ...
    'Color', 'red', 'FontSize', 8, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

axis square
title('Average fMRI Signal Difference (Normalized)');
xlabel('TR (ms)');
ylabel('nc');

saveas(gcf, [savingDir 'AveragefMRISignalDifferenceNormalized.png']);
saveas(gcf, [savingDir 'AveragefMRISignalDifferenceNormalized.pdf']);
saveas(gcf, [savingDir 'AveragefMRISignalDifferenceNormalized.fig']);





