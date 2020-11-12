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

%% Design of Low-Pass Chebyshev Filter
wlp = 1;
wls = min(-wls2, wls1);
D1 = (1/0.85)^2 - 1;
D2 = (1/0.15)^2 - 1;

Nmin = ceil(acosh(sqrt(D2/D1))/ acosh(wls/wlp));

Bk = asinh(1/sqrt(D1))/Nmin;
A0 = pi/(2*Nmin);

p1 = -sin(A0)*sinh(Bk)+1i*cos(A0)*cosh(asinh(Bk));
p2 = -sin(A0)*sinh(Bk)-1i*cos(A0)*cosh(asinh(Bk));
p3 = -sin(3*A0)*sinh(Bk)+1i*cos(3*A0)*cosh(Bk);
p4 = -sin(3*A0)*sinh(Bk)-1i*cos(3*A0)*cosh(Bk);

p = [p1 p2 p3 p4];

%for k = 1:1:Nmin
%   msg = "p" + k + " = " + p(k);
%    disp(msg);
%end

[num,den] = zp2tf([], p, (p1*p2*p3*p4)/sqrt(1+D1));
%TF with poles p(1...4) and numerator (p1*p2*p3*p4)/sqrt(1+D1) and no zeroes
% Numerator is chosen such that the DC gain in 1/sqrt(1+D1) as Nmin is even

%% Tranforming back to Band-Stop transfer function and then discrete domain
syms s z;
%analog LPF Transfer Function
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s); 
%bandstop transformation to analog BSF transfer function
analog_bsf(s) = analog_lpf((B*s)/(s*s + w0*w0));        
%bilinear transformation into dicrete BSF transfer function
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              

%% coefficients of analog low-pass filter
[nls, dls] = numden(analog_lpf(s));                  
nls = sym2poly(expand(nls));                          
dls = sym2poly(expand(dls));                         
k = dls(1);    
dls = dls/k;
nls = nls/k;

%% coeffs of analog band stop filter
[ns, ds] = numden(analog_bsf(s));               
%numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          
%collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%% coeffs of digital band stop filter
[nz, dz] = numden(discrete_bsf(z));                     
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              
k = dz(1);                                              
dz = dz/k;
nz = nz/k;

%% Magnitude response of the digital Band Stop filter in log scale
fvtool(nz,dz)                                           

%% Magnitude response of the digital Band Stop filter
[H,f] = freqz(nz,dz,10000, 260);
figure(2)
plot(f,abs(H))
title("Magnitude plot |H(w)| of Chebyshev Band stop filter")
xlabel("Frequency in kHz")
axis on
grid 


%% Magnitude response of the Analog Band stop filter
figure(3)
f = linspace(-2.5, 2.5, 50000);
h= freqs(ns,ds,f);
plot(f,abs(h))
title("Magnitude plot |H(w)| for Chebyshev Band Stop filter")
xlabel("Analog Frequency")
axis on

%% Magnitude response of the Analog Low pass filter
figure(4)
f = linspace(-4, 4, 80000);
h= freqs(nls,dls,f);
plot(f,abs(h))
ylim([0 1.2])
title("Magnitude plot |H(w)| for Chebyshev Low Pass filter")
xlabel("Analog Frequency")
axis on