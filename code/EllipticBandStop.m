clear all;

%% Unnormalised Specifications for Band-Stop Fitler
Ws1u = 43;%  in kHz
Ws2u = 63;
Wp1u = 39;
Wp2u = 67;

Sampling_Frequency = 260;

%% Normalised Specification for Band-Stop filter
Ws1n = 2*pi*Ws1u/Sampling_Frequency;
Ws2n = 2*pi*Ws2u/Sampling_Frequency;
Wp1n = 2*pi*Wp1u/Sampling_Frequency;
Wp2n = 2*pi*Wp2u/Sampling_Frequency;

%% Using Bilinear Transformation
ws1 = tan(Ws1n/2);
ws2 = tan(Ws2n/2);
wp1 = tan(Wp1n/2);
wp2 = tan(Wp2n/2);

%% Using Band-stop Transformation
B = wp2 - wp1;
w0 = sqrt(wp1*wp2);
wls1 = (B*ws1)/(w0^2 - ws1^2);
wls2 = (B*ws2)/(w0^2 - ws2^2);
wlp1 = (B*wp1)/(w0^2 - wp1^2);
wlp2 = (B*wp2)/(w0^2 - wp2^2);

%% Design of Low-Pass Elliptic Filter
wlp = 1;
wls = min(-wls2, wls1);
D1 = (1/0.85)^2 - 1;
D2 = (1/0.15)^2 - 1;

k1 = sqrt(D1/D2);
k = wlp/wls;
%% Calculating the complete elliptic integral(and its complement) of first kind for k and k_1
[K,Kc] = ellipk(k);
[K1,K1c] = ellipk(k1);

%% Calcultating the minimum degree of the filter and recalculating k and k_1 for exact(and more stringent) passband and stopband characteristics
N_min = ceil((K1c*K)/(K1*Kc));
k = ellipdeg(N_min,k1);
wls = wlp/k; % new Omega_Ls(more stringent)
L = floor(N_min/2);
r = (N_min-(2*L));
i = (1:L)';
u = (2*i-1)/N_min;
%% Finding zeroes and poles of the LPF
zeta = cde(u,k);
zeroes_lpf = (1j)./(k*zeta);
v0 = (-1j)*asne(1j/sqrt(D1),k1)/N_min;
poles_lpf = 1j*cde(u-1j*v0,k);
pole_0 = 1j*sne(1j*v0,k);
%% Finding the Transfer function
Constant_coeff = 1;
for i=1:L
    Constant_coeff = Constant_coeff*((abs(poles_lpf(i)/abs(zeroes_lpf(i))))^2);
end
if r==1
    Constant_coeff = Constant_coeff*(-pole_0);
else
    Constant_coeff = Constant_coeff/sqrt(1+D1);
end
zeroes_lpf = [zeroes_lpf,conj(zeroes_lpf)]';
poles_lpf = [poles_lpf,conj(poles_lpf),pole_0]';
[num, den] = zp2tf(zeroes_lpf, poles_lpf, Constant_coeff);           

%% Tranforming back to Band-Stop transfer function and then discrete domain
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = analog_lpf((B*s)/(s*s + w0*w0));        %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              %bilinear transformation

%% coefficients of analog low-pass filter
[nls, dls] = numden(analog_lpf(s));                  
nls = sym2poly(expand(nls));                          
dls = sym2poly(expand(dls));                         
k = dls(1);    
dls = dls/k;
nls = nls/k;

%% coeffs of analog band stop filter
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%% coeffs of digital band stop filter
[nz, dz] = numden(discrete_bsf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;

%% Magnitude response of the digital Band Stop filter in log scale
fvtool(nz,dz)                                           

%% Magnitude response of the digital Band Stop filter
[H,f] = freqz(nz,dz,10000, 260);
figure(2)
plot(f,abs(H))
title("Magnitude plot |H(w)| of Elliptic Band stop filter")
xlabel("Frequency in kHz")
axis on
grid 


%% Magnitude response of the Analog Band pass filter
figure(3)
f = linspace(-2.5, 2.5, 1000);
h= freqs(ns,ds,f);
plot(f,abs(h))
title("Magnitude plot |H(w)| for Elliptic Band Stop filter")
xlabel("Analog Frequency")
axis on

%% Magnitude response of the Analog Low pass filter
figure(4)
f = linspace(-4, 4, 2000);
h= freqs(nls,dls,f);
plot(f,abs(h))
ylim([0 1.2])
title("Magnitude plot |H(w)| for Elliptic Low Pass filter")
xlabel("Analog Frequency")
axis on