function recon=SC(w, T, Tc, d0, elp)
%depending on required output the function names can be [recon,imgcon] or conduct or recon or imgcon
%Input variables are the frequency (THz), Temp (K), critical temp Tc (K)
% gap@zero_kelvin[d0](meV), elp = (epsilon_infinity -1)
%to run the code just put a table or list of frequencies in THz , give
%correct temperature and estimated d0. 
k= 0.086173303;
h = 4.135667662;
del=d0*tanh(1.74*sqrt(Tc/T - 1));
T=k*T;
fermi = @(E)1./(exp(E/T)+1);
g=size(w);
%as w is a vector we need to pre assign the output vectors. 
imgcon=zeros(size(w));
recon=zeros(size(w));
%conduct=zeros(size(w));

kew=@(w, E)abs(E.^2 + h*w*E + del.^2)./sqrt((E + h*w).^2 - del.^2);
%defining integrand for 1st expression of sigma1
A=@(w, E)(fermi(E)-fermi(E+h*w)).*kew(w,E)./(sqrt(E.^2 - del.^2));
%integrand for 2nd expression of sigma1 (real part of conductivity)
B=@(w, E)(1-2*fermi(E+h*w)).*kew(w,E)./(sqrt(E.^2 - del.^2));
%Integrand for imaginary part of conductivity (sigma2)
B2=@(w, E)(1-2*fermi(E+h*w)).*kew(w,E)./(sqrt(del.^2 - E.^2));
for j=1:g
%calculating integral for 1st expression of sigma1 (A)
intA=integral(@(E)A(w(j),E),del,(10^5).*del);
%calculating integral for 2nd expression of sigma1 (B)
intB=heaviside(h*w(j) -2*del).*integral(@(E)B(w(j),E),del-h*w(j),-del);
%Calculating integral for imaginary part of conductivity
%imgcon=(heaviside(h*w -2*del).*integral(@(E)B2(w,E),del-h*w,del));
imgcon(j)=real(heaviside(h*w(j) -2*del).*integral(@(E)B2(w(j),E),del-h*w(j),del) + heaviside(2*del-h*w(j)).*integral(@(E)B2(w(j),E),-del,del))./(h*w(j))-w(j)*elp;
%imgcon(j)= piecewise(h*w(j) <= 2*del,(1/(h*w(j)))*integral(@(E)B2(w(j),E),del-h*w(j),del), h*w(j) > 2*del, integral(@(E)B2(w(j),E),-del,del)./(h*w(j)));
%real part of conductivity
recon(j)=(2*intA + intB)./(h*w(j));%
%total complex conductivity
%conduct(j)=complex(recon(j),imgcon(j));
end
conduct = [recon,imgcon];