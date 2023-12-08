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

% clear wp, cw, and close plots:
clear
clc
close all

% plot settings:
set(0,'defaultAxesFontSize',15)
set(0, 'DefaultLineLineWidth', 1.5)
set(0,'DefaultAxesXGrid','off')
set(0,'DefaultAxesYGrid','off')
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');

% open file:
%
% data expected in variable y (m): plots and labels are fitted to this data
% type.
[file,path] = uigetfile;
load([path,'\',file]);

% example run: (change according to test setups)
cut = 1;            % cut the transient signal from the beggining
comb = 1;           % plot the filtered results: 1, do not plot: 0
comb_n = 50;        % number of averaged periods for the comb filter freq. and amplitude
fl = 700;           % Hilbert FIR filter length (matlab envelope documentation)
p_min = 5e-4;       % minimum peak value
f_low = 58;         % lower frequency bound for plots
f_high = 66;        % higher frequency bound for plots
f_d = 1000;         % frequency resolution for BBC
hardenning = 1;     % hardenning: 1, softenning: 0

% open ORD function: plots the raw and filtered signal, the envelope of the
% filtered signal and the measured and fitted (averaged) BBCs. Output: \mu
% ~ nonlinear parameter. See further information in the documentation filte
% and in the article: https://dx.doi.org/10.2139/ssrn.4573777
mu = ORD(y,t,cut,comb,comb_n,fl,p_min,f_low,f_high,f_d,hardenning)









