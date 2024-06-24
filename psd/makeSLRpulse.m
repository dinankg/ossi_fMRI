function [rf, gz, gzr, delay] = makeSLRpulse(flip,varargin)
%makeSLRpulse make an SLR pulse
%     a wrapper to a python function(see below). See supported params below
%     in the 'parser' section. Currently it will probably only work on
%     Linux. On my system I could install the required Python library by
%     executing "pip3 install sigpy"
%
% sigpy.mri.rf.dzrf = dzrf(n=64, tb=4, ptype='st', ftype='ls', d1=0.01, d2=0.01, cancel_alpha_phs=False)
%     Primary function for design of pulses using the SLR algorithm.
%
%     Args:
%         n (int): number of time points.
%         tb (int): pulse time bandwidth product.
%         ptype (string): pulse type, 'st' (small-tip excitation), 'ex' (pi/2
%             excitation pulse), 'se' (spin-echo pulse), 'inv' (inversion), or
%             'sat' (pi/2 saturation pulse).
%         ftype (string): type of filter to use: 'ms' (sinc), 'pm'
%             (Parks-McClellan equal-ripple), 'min' (minphase using factored pm),
%             'max' (maxphase using factored pm), 'ls' (least squares).
%         d1 (float): passband ripple level in :math:'M_0^{-1}'.
%         d2 (float): stopband ripple level in :math:'M_0^{-1}'.
%         filterType (str): filter type to use, e.g. sinc (ms),
%         least-squares (ls), etc. Refer to sigpy.rf documentation.
%
%     Returns:
%         rf (array): designed RF pulse.
%
%     References:
%         Pauly, J., Le Roux, Patrick., Nishimura, D., and Macovski, A.(1991).
%         Parameter Relations for the Shinnar-LeRoux Selective Excitation
%         Pulse Design Algorithm.
%         IEEE Transactions on Medical Imaging, Vol 10, No 1, 53-65.

validPulseUses = mr.getSupportedRfUse();

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeSLRpulse';

    % RF params
    addRequired(parser, 'flipAngle', @isnumeric);
    addOptional(parser, 'system', mr.opts(), @isstruct);
    addParameter(parser, 'duration', 1e-3, @isnumeric);
    addParameter(parser, 'freqOffset', 0, @isnumeric);
    addParameter(parser, 'phaseOffset', 0, @isnumeric);
    addParameter(parser, 'timeBwProduct', 4, @isnumeric);
    addParameter(parser, 'passbandRipple', 0.01, @isnumeric);
    addParameter(parser, 'stopbandRipple', 0.01, @isnumeric);
    addParameter(parser, 'filterType', 'ls', @isstr);
    %addParamValue(parser, 'apodization', 0, @isnumeric);
    %addParamValue(parser, 'centerpos', 0.5, @isnumeric);
    % Slice params
    addParameter(parser, 'maxGrad', 0, @isnumeric);
    addParameter(parser, 'maxSlew', 0, @isnumeric);
    addParameter(parser, 'sliceThickness', 0, @isnumeric);
    addParameter(parser, 'delay', 0, @isnumeric);
    addParameter(parser, 'dwell', 0, @isnumeric); % dummy default value
    % whether it is a refocusing pulse (for k-space calculation)
    addParameter(parser, 'use', 'excitation', @(x) any(validatestring(x,validPulseUses)));
    % optional Python command
    addParameter(parser, 'pythonCmd', '', @(x)isstring(x)||ischar(x));
end
parse(parser, flip, varargin{:});
opt = parser.Results;
opt.centerpos=0.5; % fixme

if opt.dwell==0
    opt.dwell=opt.system.rfRasterTime;
end

% find/check python
if ~isempty(opt.pythonCmd)
    [status, result]=system([opt.pythonCmd ' --version']);
    if status~=0
        error(['provided python executable ''' opt.pythonCmd ''' returns an error on the version check']);
    end
    python=opt.pythonCmd;
elseif ispc()
    % on Windows we rely on the PATH settings
    [status, result]=system('python --version');
    if status==0
        python='python';
    else
        [status, result]=system('py --version');
        if status~=0
            error('python executable not found, please check your system PATH settings');
        end
        python='py';
    end
else
    % this probably only works on linux and maybe also on mac
    [status, result]=system('which python3');
    if status==0
        python=strip(result);
    else
        [status, result]=system('which python');
        if status==0
            python=strip(result);
        else
            error('python executable not found');
        end
    end
end

add_opt='';

switch opt.use
    case 'excitation'
        if opt.flipAngle <= pi/6
            ptype='st';
        else
            ptype='ex';
        end
          %add_opt=',cancel_alpha_phs=True';
    case 'refocusing'
        ptype='se';
    case 'inversion'
        ptype='inv';
    case 'saturation'
        ptype='sat';
    otherwise
        ptype='st';
end

N = round(opt.duration/opt.dwell);
% on Windows it looks like the $ and '' are not needed and ; can be used in place of \n
if ispc()
    cmd=[python ' -c "import sigpy.mri.rf;pulse=sigpy.mri.rf.dzrf(' num2str(N) ...
        ',' num2str(opt.timeBwProduct) ',ptype=''' ptype '''' ...
        ',d1=' num2str(opt.passbandRipple) ',d2=' num2str(opt.stopbandRipple) ...
        ',ftype=''' opt.filterType '''' add_opt ');print(*pulse)"'];
else
    cmd=[python ' -c $''import sigpy.mri.rf\npulse=sigpy.mri.rf.dzrf(' num2str(N) ...
        ',' num2str(opt.timeBwProduct) ',ptype=\''' ptype '\''' ...
        ',d1=' num2str(opt.passbandRipple) ',d2=' num2str(opt.stopbandRipple) ...
        ',ftype=\''' opt.filterType '\''' add_opt ')\nprint(*pulse)'''];
end

% Display the command being executed
disp('Executing Python command:');
disp(cmd);

% Execute the command and capture the output
[status, result] = system(cmd);

% Display the result of the command execution
disp('Python command output:');
disp(result);

if status~=0
    error('Executing Python command failed: %s', result);
end

lines = regexp(result,'\n','split'); % the response from the python call contains some garbage
% look for a usable result vector
for i=1:length(lines)
    try
        signal = str2num(lines{i});
        if length(signal)==N
            break;
        end
    catch
        continue;
    end
end
if length(signal)~=N
    error('could not find usable data in the response of the Python command');
end

BW = opt.timeBwProduct/opt.duration;
t = ((1:N)-0.5)*opt.dwell;
flip = abs(sum(signal))*opt.dwell*2*pi;
signal = signal*opt.flipAngle/flip;

rf.type = 'rf';
rf.signal = signal;
rf.t = t;
rf.shape_dur=N*opt.dwell;
rf.freqOffset = opt.freqOffset;
rf.phaseOffset = opt.phaseOffset;
rf.deadTime = opt.system.rfDeadTime;
rf.ringdownTime = opt.system.rfRingdownTime;
rf.delay = opt.delay;
if ~isempty(opt.use)
    rf.use=opt.use;
end
if rf.deadTime > rf.delay
    rf.delay = rf.deadTime;
end

if nargout > 1
    assert(opt.sliceThickness > 0,'SliceThickness must be provided');
    if opt.maxGrad > 0
        opt.system.maxGrad = opt.maxGrad;
    end
    if opt.maxSlew > 0
        opt.system.maxSlew = opt.maxSlew;
    end

    amplitude = BW/opt.sliceThickness;
    area = amplitude*opt.duration;
    gz = mr.makeTrapezoid('z', opt.system, 'flatTime', opt.duration, ...
        'flatArea', area);
    gzr= mr.makeTrapezoid('z', opt.system, 'Area', -area*(1-opt.centerpos)-0.5*(gz.area-area));
    if rf.delay > gz.riseTime
        gz.delay = ceil((rf.delay - gz.riseTime)/opt.system.gradRasterTime)*opt.system.gradRasterTime; % round-up to gradient raster
    end
    if rf.delay < (gz.riseTime+gz.delay)
        rf.delay = gz.riseTime+gz.delay; % these are on the grad raster already which is coarser
    end
end

% v1.4 finally eliminates RF zerofilling
% if rf.ringdownTime > 0
%     tFill = (1:round(rf.ringdownTime/1e-6))*1e-6;  % Round to microsecond
%     rf.t = [rf.t rf.t(end)+tFill];
%     rf.signal = [rf.signal, zeros(size(tFill))];
% end


if rf.ringdownTime > 0 && nargout > 3
    delay=mr.makeDelay(mr.calcDuration(rf)+rf.ringdownTime);
end

end
