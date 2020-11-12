clear all;

%% Unnormalised Specifications for Band-Pass Fitler
Wp1u = 48.4;%  in kHz
Wp2u = 68.4;
Ws1u = 44.4;
Ws2u = 72.4;

Sampling_Frequency = 330;

%% Normalised Specification for Band-Pass filter
Wp1n = 2*pi*Wp1u/Sampling_Frequency;
Wp2n = 2*pi*Wp2u/Sampling_Frequency;
Ws1n = 2*pi*Ws1u/Sampling_Frequency;
Ws2n = 2*pi*Ws2u/Sampling_Frequency;

%% Using Bilinear Transformation
wp1 = tan(Wp1n/2);
wp2 = tan(Wp2n/2);
ws1 = tan(Ws1n/2);
ws2 = tan(Ws2n/2);

%% Using Band-pass Transformation
B = wp2 - wp1;
w0 = sqrt(wp1*wp2);
wlp1 = (wp1^2 - w0^2)/(B*wp1);
wlp2 = (wp2^2 - w0^2)/(B*wp2);
wls1 = (ws1^2 - w0^2)/(B*ws1);
wls2 = (ws2^2 - w0^2)/(B*ws2);

%% Design of Low-Pass Butterworth filter
wlp = 1;
wls = min(-wls1,wls2);
D1 = (1/0.85)^2 - 1;
D2 = (1/0.15)^2 - 1;

Nmin = ceil(log(sqrt(D2/D1))/ log(wls/wlp));
Wcmin = wlp/((D1)^(1/(2*Nmin)));
Wcmax = wls/((D2)^(1/(2*Nmin)));

Wc = 1.07; %Chosen in between Wcmin and Wcmax

p = zeros([Nmin,1]);
for k = 1:1:Nmin
    theta = (2 * k + 7)*pi;
    p(k) = Wc* exp(1i*theta/16);
    msg = "p" + k + " = " + p(k);
    %disp(msg);
end

%TF with poles p(1...8) and numerator Wc^Nmin and no zeroes
[num,den] = zp2tf([], p, Wc^Nmin); 
%numerator is chosen to make the DC Gain 1
                                    
%% Tranforming back to Band-Pass transfer function and then discrete domain
syms s z;
%analog LPF Transfer function
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);
%bandpass transformation to get analog BPF Transfer function
analog_bpf(s) = analog_lpf((s*s + w0*w0)/(B*s)); 
%bilinear transformation to get discrete BPF system function 
discrete_bpf(z) = analog_bpf((z-1)/(z+1));         

%% coefficients of analog low-pass filter
[nls, dls] = numden(analog_lpf(s));                  
nls = sym2poly(expand(nls));                          
dls = sym2poly(expand(dls));                         
k = dls(1);    
dls = dls/k;
nls = nls/k;

%% coefficients of analog band-pass filter
[ns, ds] = numden(analog_bpf(s));                  
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                         
k = ds(1);    
ds = ds/k;
ns = ns/k;

%% coeffs of discrete band-pass filter
[nz, dz] = numden(discrete_bpf(z));                                         
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              
k = dz(1);                                              
dz = dz/k;
nz = nz/k;

%% Magnitude response of the digital Band Pass filter in log scale
fvtool(nz,dz)                                           

%% Magnitude response of the digital Band Pass filter
[H,f] = freqz(nz,dz,10000, Sampling_Frequency);
figure(2)
plot(f,abs(H))
title("Magnitude plot |H(w)| for Butterworth Band Pass filter")
xlabel("Frequency in kHz")
axis on

%% Magnitude response of the Analog Band pass filter
figure(3)
f = linspace(-1.5, 1.5, 30000);
h= freqs(ns,ds,f);
plot(f,abs(h))
title("Magnitude plot |H(w)| for Butterworth Band Pass filter")
xlabel("Analog Frequency")
axis on

%% Magnitude response of the Analog Low pass filter
figure(4)
f = linspace(-4, 4, 80000);
h= freqs(nls,dls,f);
plot(f,abs(h))
ylim([0 1.2])
title("Magnitude plot |H(w)| for Butterworth Low Pass filter")
xlabel("Analog Frequency")
axis on