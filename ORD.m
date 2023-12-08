%------------------------------------------------------------------------------%
% Copyright (c) 2023 Zoltan Gabos
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software, to deal in the Software without restriction, including without 
% limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% The watermarks in the figures must be included in all copies or substantial 
% portions of the Software and in the exported figures.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%------------------------------------------------------------------------------%


function mu = ORD(y,t,cut,comb,comb_n,fl,p_min,f_low,f_high,f_d,hardenning)

% cut initial transient: --------------------------------------------------
if cut == 1
    % pick the cutting time instance:
    figure(2)
    plot(y); ylabel('$x$ (m)'); box on; % plotting the signal
    xlabel('$t$ (s)');
    [n,p] = ginput(1);                  % hand picking end of transient
    close(figure(2));                   % close the figure
    
    % cut the signal:
    y = y(ceil(n):end);
    t = t(1:end-ceil(n)+1);
end
%--------------------------------------------------------------------------

% comb filter: ------------------------------------------------------------
y_comb = y;                             % relabel the data for filtering
var_1 = smoothdata(y);                  % smoothing for zero crossing
t_1 = [];                               % crossing time instances

% cut the signal at the first crossing:
for i=1:length(var_1)
    
    if i+1 > length(var_1)              % abort iteration
        break
    end
    
    if sign(var_1(i))<sign(var_1(i+1))  % crossing from - to +
        y_comb = y_comb(i:end);         % signal cutting
        t = t(1:end-i+1);
        break
    end
end

% plot original signal (onli if comb == 1)
if comb == 1
    figure(2)
    hold on
    box on;
    xlabel('$t$ (s)');
    ylabel('$x$ (m)');
    plot(t,y_comb)
end

% comb filter:
var_2 = envelope(y_comb,fl,'peak');     % envelope of signal
var_1 = smoothdata(y_comb);             % smoothing for zero crossing
t_indices = [];                         % crossing indices
for i=1:length(var_1)
    
    if i+1 > length(var_1)              % abort iteration
        break
    end
    
    if sign(var_1(i))<sign(var_1(i+1))  % detect zero crossings
        t_1(end+1) = t(i);              % store zero crossings (ZC)
        t_indices(end+1) = i;           % store zero crossing indices
    end
end

T_1 = zeros(1,length(t_1));             % time periods from ZC

for i = 1:length(t_1)-1                 % time period calculation
    T_1(i) = t_1(i+1)-t_1(i);           % storing time periods
end

f_1 = 1./T_1(1:end-1);                  % instantaneous frequencies (Hz)
f_sim = mean(f_1(1:comb_n));            % comb filter frequency (Hz)
f_amp = mean(var_2(1:t_indices(comb_n))); % comb filter amplitude (4243)
sim_sin = f_amp*sin(f_sim*2*pi.*t);     % simulated signal for c. filter
y_comb = y_comb-sim_sin;                % applying comb filter

% plot filtered signal (only if comb == 1)
if comb == 1
    plot(t,y_comb)
end
legend('raw signal','filtered signal')
%--------------------------------------------------------------------------

% signal sectioning for transients: ---------------------------------------
y_peak = var_2;                         % envelope of the signal

% plot signal envelope to pick the maximum amplitudes:
figure(3)
plot(y_peak)                            % plot signal amplitudes

[n,p] = ginput(1);                      % pick amplitude treshold
findpeaks(y_peak,'MinPeakHeight',p)     % find peaks above the treshold
% store the peak values and locations:
[pks,locs] = findpeaks(y_peak,'MinPeakHeight',p_min); % min Peak value: 5e-4
% make temporary data:
locs_temp = locs;
pks_temp = pks;

% !!! 4 transients are considered !!! - this is not optional.
[n_t,p_p] = ginput(4);                  % pick initial transient amplitudes
locs = zeros(1,4);
pks = zeros(1,4);
% store initial transient amplitudes:
for i = 1:length(n_t)
    % find closes point for the initial transients:
    n_temp = find(min(abs(locs_temp-n_t(i)))==abs(locs_temp-n_t(i)));
    locs(i) = locs_temp(n_temp);
    pks(i) = pks_temp(n_temp);
end

% cut the transients ends before the next perturbation: (again, 4 transients
% are possible, no more no less)
[n,p] = ginput(4);
y1 = y_comb(locs(1):floor(n(1)));
y2 = y_comb(locs(2):floor(n(2)));
y3 = y_comb(locs(3):floor(n(3)));
y4 = y_comb(locs(4):floor(n(4)));
t1 = t(locs(1):floor(n(1)));
t2 = t(locs(2):floor(n(2)));
t3 = t(locs(3):floor(n(3)));
t4 = t(locs(4):floor(n(4)));
n_end = ceil(n);

% store the transient signals:
y_a = {y1,y2,y3,y4};
t_a = {t1,t2,t3,t4};

%--------------------------------------------------------------------------

% calculate the BBCs: -----------------------------------------------------
for j = 1:4
    % temporary signal data:
    y_temp = y_a{j};
    t_temp = t_a{j};
    
    A_inst = envelope(y_temp,fl,'peak');    % instantaneous amplitude
    t_1 = [];                               % crossing time instances
    t_2 = [];                               % crossing indices
    y_temp = smoothdata(y_a{j});            % smoothing for zero crossing
    
    % calculate instantaneous frequencies:
    for i=1:length(y_temp)
        
        if i+1 > length(y_temp)             % abort iteration
            break
        end
        
        if sign(y_temp(i))<sign(y_temp(i+1))% detect zero crossings
            t_1(end+1) = t_temp(i);         % store zero crossings
            t_2(end+1) = i;                 % store zero crossing indices
        end
    end
    
    T_1 = zeros(1,length(t_1));             % time periods from ZC
    
    for i = 1:length(t_1)-1                 % time period calculation
        T_1(i) = t_1(i+1)-t_1(i);
    end
    
    f_1 = 1./T_1(1:end-1);                  % instantaneous frequencies (Hz)

    for i = 1:length(f_1)                   % cut the signal if the instantaneous 
                                            % frequency is too far from the forcing frequeny
        if abs(f_1(i)-f_sim) > 0.3*f_sim
            f_1 = f_1(1:i-1);
            A_inst = A_inst(1:t_2(i-1));
            break
        end
    end
%--------------------------------------------------------------------------

% interpolating BBCs: -----------------------------------------------------
    
    % cut the instantaneous frequencies accoriding to A_inst:
    f2 = zeros(1,length(A_inst));
    f2(1,1:t_2(1)-1) = f_1(1);
    
    % linearly interpolate f2 to get a continuous instantaneous frequency
    % data:
    for i = 1:length(f_1)-1
        f2(1,t_2(i):t_2(i+1)-1) = linspace(f_1(i),f_1(i+1),t_2(i+1)-t_2(i));
    end
    
    f2(1,t_2(end-1):end) = f_1(end-1);      % end correction of f_inst
    figure(4); hold on; plot(f2,A_inst)     % plot the measured BBCs
    xlabel('frequency (Hz)');
    ylabel('$A$ (m)');

    % curve fitting:
    figure(4)
    xlim([f_low f_high])

    % pick the instant where the transient goes back to the operational
    % frequency:
    [n,p] = ginput(1);
    
    % cut the data if the instantaneous frequency is less than the picked
    % value:
    if hardenning == 1
        for i = 1:length(f2)
            if f2(i) <= n
                f2 = f2(1:i-2);
                A_inst = A_inst(1:i-2);
                break
            end
        end
    else
        for i = 1:length(f2)
            if f2(i) >= n
                f2 = f2(1:i-2);
                A_inst = A_inst(1:i-2);
                break
            end
        end
    end
    
    cut_f = n;                          % cutting frequency
    if hardenning == 1
        span = linspace(cut_f,f_high,f_d);  % mesh for BBC interp
    else
        span = linspace(f_low,cut_f,f_d);
    end
    
    % shift data for interpolation: (shift back occurs after
    % interpolation), the curve fitting toolbox work better if we shift the
    % data for small amplitude oscillations.
    for i = -10:10-1
        if max(y)>10^i && max(y)<10^(i+1)   % find the power of 10
            n_10 = i;
            break
        end
    end
    
    % shift occurs here:
    if n_10 < 0
        n_lab = 10^(n_10+13);
    else
        n_lab = 1;
    end
    
    % interpolation:
    [xData, yData] = prepareCurveData( f2', A_inst );
    if hardenning == 1
        ft = fittype( ['2*sqrt(((x*2*pi)^2-(',num2str(cut_f),'*2*pi)^2)/abs(3*mu*',num2str(n_lab),'))'], 'independent', 'x', 'dependent', 'y' );
    else
        ft = fittype( ['2*sqrt(((x*2*pi)^2-(',num2str(cut_f),'*2*pi)^2)/-abs(3*mu*',num2str(n_lab),'))'], 'independent', 'x', 'dependent', 'y' );
    end
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = 0.141886338627215; % optional.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    mu = coeffvalues(fitresult);
    if hardenning == 1
        if j ~= 1 % average fitting from transients
            b = (b + 2*sqrt(((span*2*pi).^2-(cut_f*2*pi)^2)./abs(3*mu*n_lab)))/2;
        else
            b = 2*sqrt(((span*2*pi).^2-(cut_f*2*pi)^2)./abs(3*mu*n_lab));
        end
    else
        if j ~= 1
            b = (b + 2*sqrt(((span*2*pi).^2-(cut_f*2*pi)^2)./-abs(3*mu*n_lab)))/2;
        else
            b = 2*sqrt(((span*2*pi).^2-(cut_f*2*pi)^2)./-abs(3*mu*n_lab));
        end
    end % here n_lab is to shift back the data.
    clc                     % delete size watinings
end
%--------------------------------------------------------------------------

% plot BBC estimation: ----------------------------------------------------
figure(4)
hold on
box on
plot(span,b,'--','color',"#77AC30")
if hardenning == 1
    legend('sample 1','sample 2','sample 3','sample 4','BBC fit','Location','northwest')
else
    legend('sample 1','sample 2','sample 3','sample 4','BBC fit','Location','northeast')
end

% figure 3 labels:
figure(3)
plot(t,y_peak);
hold on
box on
plot(locs*(t(2)-t(1)),pks,'g+')
plot(n_end*(t(2)-t(1)),y_peak(n_end),'r+')
xlabel('$t$ (s)');
ylabel('$A$ (m)');
%--------------------------------------------------------------------------
if hardenning == 1
    mu = abs(3*mu*n_lab); % output: hardenning nonlin param
else
    mu = -abs(3*mu*n_lab);% output: softenning nonlin param
end

if comb ==1
    c_start = 2;
else
    c_start = 3;
end

for i = c_start:4
figure(i)
annotation('textbox',[0.13 0.1 0.1 0.1], 'Interpreter', 'Latex', ...
    'String','\copyright BME GPK MM, Gabos 2023','EdgeColor','none',...
    'Color', [0.5 0.5 0.5])
end

end % end of function.
