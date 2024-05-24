%% 3D OSSI pulse sequence coded with PulSeq
%{
Dinank Gupta and Mariama Salifu
University of Michigan, May 2024
TODO: Change RF to SLR/SPSP
%}

clc
clear

%% Constants:
gamT = 4.2576e7;   % Hz/Tesla
gamG = gamT/1e4;    % Hz/Gauss

%% Basic Scan Parameters
ACTUAL_SCAN = 1;
seq_params.scanner = 'inside';
seq_params.trig = 0;
seq_params.nc = 10;
seq_params.TR = 40e-3;
seq_params.ntp = 10;
seq_params.nDummyLoops = 2;
seq_params.nshot_spiral = 4;
seq_params = getparams(seq_params);
seq_params.oversamp=100;%Fully sample these many spiral samples.

sys=mr.opts('maxGrad',seq_params.maxgrad,'gradUnit', ...
    'mT/m','maxSlew',seq_params.maxslew, 'slewUnit', 'T/m/s');
sys.adcDeadTime=1e-5;
% sys.rfDeadTime = 1e-4;
%% RF pulse segment
seq = mr.Sequence(sys);

% RF pulse
[rf, gz_rf,gzReph] = mr.makeSincPulse(seq_params.alpha*pi/180, 'Duration', 3e-3, ...
    'SliceThickness', seq_params.slicethickness, 'apodization', 0.42, ...
    'timeBwProduct', 4, 'system', sys);
gzPreph = gzReph;
rf.delay = rf.deadTime;
gz_rf.delay = ceil((mr.calcDuration(gzReph))/sys.gradRasterTime)*sys.gradRasterTime;
rf.delay = gz_rf.delay+gz_rf.riseTime; %-5e-7
gzReph.delay =  mr.calcDuration(gz_rf);
gzcomb = mr.addGradients({gzPreph,gz_rf,gzReph},'system',sys);
%% Spiral segment
% Parameters for spiral readouts in toppe compatible units
ge_sys = [];
maxgrad_gcm = grad_convertion(seq_params.maxgrad);
maxslew_gcms = slew_convertion(seq_params.maxslew,'s');
ge_sys.fov = seq_params.fov.*100;
ge_sys.N = seq_params.N;
ge_sys.raster = sys.gradRasterTime;
ge_sys.maxSlew = slew_convertion(seq_params.maxslew,'ms');  % assumes G/cm/ms
ge_sys.maxGrad = maxgrad_gcm;  % assumes G/cm

% create spiral readout with multiple shot
[gx1,gy1,t,spiral_readout_length] = makevdspiral2(seq_params.fov(1)*100,seq_params.N(1), ...
    seq_params.nshot_spiral,seq_params.oversamp,ge_sys.maxGrad, ...
    0.9*seq_params.maxslew ,sys.gradRasterTime);
%making it divisible by 10
gx1 = [gx1;zeros(10-rem(length(gx1),10),size(gx1,2))];
gy1 = [gy1;zeros(10-rem(length(gy1),10),size(gy1,2))];

gx1_Hzcm = gx1*gamG*100;gy1_Hzcm = gy1*gamG*100; %Converting to Hz/m
% Making struct with all spiral shots.
gx = cell(1,seq_params.nshot_spiral);
gy = cell(1,seq_params.nshot_spiral);
% last_sample_time = sys.gradRasterTime*(size(gx1,1)-1);
for i = 1:seq_params.nshot_spiral
    gx{i} = mr.makeArbitraryGrad('x',gx1_Hzcm(:,i),sys,'delay',.1e-3);
    gy{i}  = mr.makeArbitraryGrad('y',gy1_Hzcm(:,i),sys,'delay',.1e-3);
    sum(isnan(gx1_Hzcm(:,i)))
    sum(isnan(gy1_Hzcm(:,i)))


end

%% kz-encoding segment
dkz = 1/seq_params.fov(3); %1/m
res = seq_params.fov(3)/seq_params.N(3);%m
seq_params.kzmax = 1/res;%1/m
z_enc_trap = mr.makeTrapezoid('z','Area',seq_params.kzmax,'system',sys);
%^CHECK. From the comments in the function it seems to want 1/kzmax.

% Grads will go from -Gz_area/2 to Gz_area/2
zp = (-(seq_params.N(3)-1)/2:(seq_params.N(3)-1)/2)/((seq_params.N(3)-1)/2)/2; %zpe scaling for each TR
%List of gradient area:
seq_params.zp_scale = zp;%/zp(1); %Scaling of the kz traps

%% ADC
adc = mr.makeAdc(spiral_readout_length,sys,'Dwell',sys.gradRasterTime,'delay',.1e-3);
%% Define delays
TRmin = mr.calcDuration(gzcomb)+mr.calcDuration(gx{1}) + 2*mr.calcDuration(z_enc_trap);
delayTR = ceil((seq_params.TR - TRmin)/seq.gradRasterTime)*seq.gradRasterTime;

%% Putting it all together
b1ph = (pi/seq_params.nc).*((0:seq_params.nc-1).^2).';  % Phase cycling - Only works for even values of Nc
seq_params.cal_tp = 1; % Time point to use for calibration (no kz)

ndisdaq = seq_params.nDummyLoops;
adc_flag=0;
for itp = 0:seq_params.ntp % 0th timepoint will serve as disdaq.
    for nd=1:ndisdaq % n_disdaq will be set to 1 after adc is on.
        if(adc_flag==0) %No gy blip, no adc
            segmentID = 0;
        elseif(adc_flag==1)
            segmentID=seq_params.nshot_spiral;
        end
        for ishots = 1:seq_params.nshot_spiral

            for ikz = 1:seq_params.N(3)

                segmentID = ishots+1 + (adc_flag*seq_params.nshot_spiral);
                for iNc = 1:seq_params.nc
                    % transmit and receive phase (OSS phase cycling)
                    rf.phaseOffset  = b1ph(iNc);
                    adc.phaseOffset = b1ph(iNc);
                    %RF
                    seq.addBlock(rf,gzcomb, mr.makeLabel('SET', 'TRID', segmentID))%
                    %kz-encode
                    seq.addBlock(mr.scaleGrad(z_enc_trap,seq_params.zp_scale(ikz)));
                    if adc_flag
                        seq.addBlock(gx{ishots},gy{ishots},adc)
                    else
                        seq.addBlock(gx{ishots},gy{ishots})
                    end
                    %kz-encode balance
                    seq.addBlock(mr.scaleGrad(z_enc_trap,-1*seq_params.zp_scale(ikz)));
                    seq.addBlock(mr.makeDelay(delayTR));

                end
            end
        end
        if(nd==ndisdaq && adc_flag==0)% After the disdaq loop ends once, set adc to 1.
            adc_flag=1;
            ndisdaq=1;%After disdaqs are done, set the ndisdaq to 1
            break
        end
    end

end

%%
% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

seq.plot('timerange',[0 seq_params.TR*seq_params.nc]); %*seq_params.nc
%
% k-space trajectory calculation
% [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
%
% % plot k-spaces
% figure; plot(t_ktraj, ktraj'); title('k-space components as functions of time'); % plot the entire k-space trajectory
% figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
% hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); title('2D k-space');
%% Output for execution
% seq.setDefinition('gradRasterTime', seq_params.gradRasterTime);
% seq.setDefinition('adcRasterTime', seq_params.adcRasterTime);
% seq.setDefinition('blockDurationRaster', seq_params.gradRasterTime);
seq.setDefinition('FOV', seq_params.fov);
seq.setDefinition('Name', 'OSSI');
seq.write('external.seq');

%% Seeing in GE format to make sure units look ok
sysGE = toppe.systemspecs('maxGrad', 10, ... % G/cm
    'maxSlew', 20, ... % G/cm/ms
    'maxRF', 0.2, ... % Gauss. Must be >= peak RF in sequence
    'rfDeadTime',100,...
    'adcDeadTime', 0, ... % us. Half of 40us since applied both before + after ADC window.
    'psd_rf_wait', 148, ... % RF/gradient delay (us)
    'psd_grd_wait', 156); % ADC/gradient delay (us)

seq2ge('external.seq',sysGE,'ossi.tar');

system('tar -xvf ossi.tar');

figure; toppe.plotseq(sysGE,'timeRange',[0 seq_params.TR*seq_params.nc]); % Plot 1 frame

%% Plotting multiple TRs:
%{
for k=1:5:300
    toppe.plotseq(sysGE,'timeRange',[k*seq_params.TR, (k+1)*seq_params.TR]);
    drawnow
end
%}
%%
if ACTUAL_SCAN
    exp_number = 0;
    while true
        exp_number = exp_number + 1;
        exp_name = ['experiments/',datestr(now,'yyyy-mm-dd_'),seq_params.scanner,'_exp',num2str(exp_number)];
        if ~isfolder(exp_name)
            mkdir(exp_name)
            break
        end
    end
    [status] = copyfile('ossi.tar', ['./',exp_name,'/ossi.tar']);
    [status] = copyfile('external.seq', ['./',exp_name,'/ossi.seq']);
    system(append('tar -xvf ','./',exp_name,'/','ossi.tar',...
        ' -C ',exp_name,'/'));

    % save('seq.mat','seq')
    % [status] = copyfile('seq.mat', ['./',exp_name,'/seq.mat']);
    save('variables')
    [status] = copyfile('variables.mat', ['./',exp_name,'/variables.mat']);


    toppe.utils.coppe('target',seq_params.scanner,'cv',788,'use_pw',1);
end

