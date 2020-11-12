clear all;
f_samp = 260; % in kHz

%Unnormalised Specifications for Band-Stop Fitler
fp1 = 39;
fs1 = 43;
fs2 = 63;
fp2 = 67;
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

N_min = ceil((A-8)/(2*2.285*(ft/f_samp)*2*pi));       %empirical formula for N_min

%Window length for Kaiser Window 
n=(2*N_min + 1) + 14; % By trail and error method, we increase n
              %least n satisfying all the conditions was 55, i.e. (2*Nmin +1)+14

%Ideal bandstop impulse response of length "n"
IdealBSF =  IdealLPF(pi,n) -IdealLPF(((fs2+fp2)/f_samp)*pi,n) + IdealLPF(((fs1+fp1)/f_samp)*pi,n);

%Kaiser Window of length "n" with shape paramter alpha calculated above
kaiser_win = (kaiser(n,alpha))';

FIR_BandStop = IdealBSF .* kaiser_win
fvtool(FIR_BandStop);         %Magnitude and Phase response in normalized domain

%Magnitude response in Unnormalised frequencies
[H,f] = freqz(FIR_BandStop,1,10000, f_samp);

figure(2)
plot(f,abs(H))
title("Magnitude plot |H(w)| for FIR Band Stop filter")
xlabel("Frequency in kHz")
axis on
grid

%Impulse response of the Band Stop Filter
figure(3)
plot(FIR_BandStop)
title("Impulse response of the FIR Band Stop Filter")
xlabel('samples')
axis on
grid