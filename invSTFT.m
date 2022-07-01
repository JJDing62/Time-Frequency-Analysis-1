% for the inverse STFT
clear
[x1, fs] = audioread('D:\新文件\TFW\Chord.wav');
x1=x1(:,1).';
L=length(x1);
sgm=500;
tau=[0:length(x1)-1]/fs;
dt=0.01; df=1; dtau=1/fs;
te=floor(tau(end)/dt)*dt;
t=[0:dt:te];
B=1.9143/sgm^.5; Q=floor(B/dtau); Nq=2*Q+1;
N=1/dtau/df;
if N < Nq
    'N is too small'
    return
end
N1=round(N/2);
f1=[N1-N:N1-1]*df;
f= 0:df:1000;
fa=(f(1)-f1(1))/df+1;  fb=(f(end)-f1(1))/df+1;
y=Gabor51(x1,tau,t,f1,sgm);  
% inverse STFT
tn=round((t-tau(1))/dtau)+1;  Lt=length(t);
m=round(f1/df);   
m1=j*2*pi/N*m;
gs=exp(-pi*dtau^2*[-Q:Q].^2*sgm)*sgm^0.25*dtau;
x1u=zeros(1,L);
wiv=zeros(1,L);
for a1=1:Lt
    y2=y(:,a1).'.*exp(-m1*(Q-tn(a1)+1));
    y2=y2([N-N1+1:N,1:N-N1]);
    x2=ifft(y2);
    ta1=Q+1-min(tn(a1)-1,Q);
    ta2=Q+1+min(L-tn(a1),Q);
    tiv=[ta1:ta2]-Q-1+tn(a1);
    wv=gs(ta1:ta2);
    x3=x2(ta1:ta2)./wv;
    x1u(tiv)=x1u(tiv)+x3.*wv;
    wiv(tiv)=wiv(tiv)+wv; 
end
x1v=x1u./wiv;
norm(x1-x1v)/norm(x1)


figure(1)
image(t,f,abs(y(fa:fb,:))/max(max(abs(y)))*600)
colormap(gray(256))
set(gca,'Ydir','normal')
figure(2)
subplot(211)
plot(tau,x1)
xlim([min(tau),max(tau)])
subplot(212)
plot(tau,real(x1v))
xlim([min(tau),max(tau)])
