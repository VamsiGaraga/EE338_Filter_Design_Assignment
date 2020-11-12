function hn = IdealLPF(wc,M) %Length M, Passband = wc
mid = (M-1)/2;  
n = [0:1:(M-1)];
m = n - mid + eps; % here a small value is being added so that h[0] will become h[eps] and will not diverge
hn = sin(wc*m)./(pi*m);

