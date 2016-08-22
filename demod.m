% demod: demodulates double precision TDO oscillations
% To play with resolution: change L, npad, and/or nstep

dt_fast = 5e-8; % timebase for TDO signal
dt_slow = 1e-6; % timebase for everything else
[file, folder] = uigetfile('./LoganBVH/PFL/*.txt');
fullFile = fullfile(folder, file);
% data = [Bdot Idot Proj1 Proj2 PDO Field BProj1 BProj2 Field_fixed]
data = dlmread(fullFile,'\t', 5 ,1);
npnts = size(data,1); % number of points in TDO signal
npnts_slow = npnts*dt_fast/dt_slow; % number of points in everything else
time_fast = 0:dt_fast:dt_fast*(npnts-1);
time_slow = 0:dt_slow:dt_slow*(npnts_slow-1);
% Take an FFT of length L+npad and then advance by nstep
% Zero-pad the end of the time-domain signal to help resolution
L = 512; % length of oscillation data to be FFT'd
npad = 2048-L; % zero-padding (for speed should be 2^N-L, or not)
nstep = 64; % length of step
Fs = 1/dt_fast; % sampling frequency
freq = zeros(npnts/nstep,1);
time = zeros(npnts/nstep,1);
pad = zeros(npad,1);
% Perform sliding FFT
for i = 1 : npnts/nstep-nstep
    fftdata = [data(i*nstep:i*nstep+L-1,5); pad];
    tdo_fft = fft(fftdata);
    %tdo_fft = fft(data(i*nstep:i*nstep+L-1,5)); % if no zero-padding
    P2 = abs(tdo_fft/L);
    P1 = P2(1:L/2+1);
    %f = Fs*(0:(L/2))/L; % if no zero-padding
    f = Fs*(0:(L/2))/(L+npad);
    [~, ind] = max(P1);
    freq(i) = f(ind);
    time(i) = dt_fast*nstep*i;
end
freq = freq(1:end-nstep);
time = time(1:end-nstep);
% find point in time_slow that matches the last value in time 
[~, t] = find(ismembertol(time_slow,time(end),1e-5)==1);
timeB = time_slow(1:t)';
field = data(1:t,9); 
% plot data
figure
yyaxis left
plot(time,freq)
yyaxis right
plot(timeB, field)
% save data
dlmwrite([file(1:end-7) '_demod.txt'], [time freq]);
dlmwrite([file(1:end-7) '_field.txt'], [timeB field]);


