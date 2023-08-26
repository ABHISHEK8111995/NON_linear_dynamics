function F2=modal_freq_step(d,di,D,l,i_len,g,den,ip)
clc
clear
syms w
d=150;
di=15;
D=160;
l=1000;
i_len=160;
g=76923; %mpa
den=7.850*10^-6; %kg/mm^3
ip=2250000;
% l=50;
L1=(l-i_len)/5;
L2=(i_len)/5;

j1=pi*(d^4-di^4)/32;
j2=pi*(D^4-di^4)/32;
% g=76923;
ma1=(den*j1*L1/6);
ma2=(den*j2*L2/6);
% disp('............FINITE ELEMENT METHOD.........')

m=cell(1,10);
k=cell(1,10);

%element 1
m{1}=ma1*[2 1;1 2];
k1=(g*j1/L1)*[1 -1;-1 1];
k{1}=k1;

%element 2

m{2}=ma1*[2 1;1 2];
k2=(g*j1/L1)*[1 -1;-1 1];
k{2}=k2;

%element 3
m{3}=ma1*[2 1;1 2];
k3=(g*j1/L1)*[1 -1;-1 1];
k{3}=k3;


%element 4
m{4}=ma1*[2 1;1 2];
k4=(g*j1/L1)*[1 -1;-1 1];
k{4}=k4;

%element 5
m{5}=ma1*[2 1;1 2];
k5=(g*j1/L1)*[1 -1;-1 1];
k{5}=k5;


%element 6
m{6}=ma2*[2 1;1 2];
k6=(g*j2/L2)*[1 -1;-1 1];
k{6}=k6;

%element 7
m{7}=ma2*[2 1;1 2];
k7=(g*j2/L2)*[1 -1;-1 1];
k{7}=k7;


%element 8
m{8}=ma2*[2 1;1 2];
k8=(g*j2/L2)*[1 -1;-1 1];
k{8}=k8;


%element 9
m{9}=ma2*[2 1;1 2];
k9=(g*j2/L2)*[1 -1;-1 1];
k{9}=k9;


%element 10
m{10}=ma2*[2 1;1 2]+[0 0;0 ip];
k10=(g*j2/L2)*[1 -1;-1 1];
k{10}=k10;



%global Mass matrix
MASS_MATRIX=zeros(11,11);
for i=1:10
    MASS_MATRIX(i:i+1,i:i+1)=MASS_MATRIX(i:i+1,i:i+1)+m{i};
end
  
%global stiffness matrix
STIFFNESS_MATRIX=zeros(11,11);
for i=1:10
    STIFFNESS_MATRIX(i:i+1,i:i+1)=STIFFNESS_MATRIX(i:i+1,i:i+1)+k{i};
end

%boundary condition
MASS_MATRIX=MASS_MATRIX(2:11,2:11)
STIFFNESS_MATRIX=STIFFNESS_MATRIX(2:11,2:11)


MASS_MATRIX=vpa(-w^2*MASS_MATRIX);
A2=(MASS_MATRIX+STIFFNESS_MATRIX);


%solving for frequency
A2=det(A2);
frequency2=solve(A2,w);
freq=double(frequency2)
F2=min(freq(freq>0))

end



