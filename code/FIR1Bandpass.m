clear all;
f_samp = 330; % in kHz

%Unnormalised Specifications for Band-Pass Fitler
fs1 = 44.4;
fp1 = 48.4;
fp2 = 68.4;
fs2 = 72.4;
ft = 4;

%Calculating A
A = -20*log10(0.15);

%Getting appropriate Alpha using A
if(A < 21)
    alpha = 0;
elseif(A <51)
    alpha = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    alpha = 0.1102*(A-8.7);
end

N_min = ceil((A-8)/(2*2.285*(ft/f_samp)*2*pi));       
%empirical formula for N_min

%Window length for Kaiser Window
n = (2*N_min + 1) + 18;
% By trail and error method, we increase n by 2, the least n satisfying all the conditions was 69, i.e. (2*N_min+ 1) + 18


%Ideal bandpass impulse response of length "n"
IdealBPF = IdealLPF(((fs2+fp2)/f_samp)*pi,n) - IdealLPF(((fs1+fp1)/f_samp)*pi,n);

%Kaiser Window of length "n" with shape paramter Alpha calculated above
kaiser_win = (kaiser(n,alpha))';

FIR_BandPass = IdealBPF .* kaiser_win;

%Magnitude and Phase response in normalized domain
fvtool(FIR_BandPass); 

%Magnitude response in Unnormalised frequencies
[H,f] = freqz(FIR_BandPass,1,10000, f_samp);

figure(2)
plot(f,abs(H))
title("Magnitude plot |H(w)| for FIR Band Pass filter")
xlabel("Frequency in kHz")
axis on
grid

%Impulse response of the Band PassFilter
figure(3)
plot(FIR_BandPass)
title("Impulse response of the FIR Band Pass Filter")
xlabel('samples')
axis on
grid