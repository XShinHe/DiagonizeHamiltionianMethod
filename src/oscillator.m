% Harmonic oscillaor system

clear;clf;

% define the task
workpath='../test/doublewell/'
jobtitle='doublewell'

%-- parameters of 1-D infinite well system
%----- particle mass
m=1
%----- width of well
L=20;
x=-L/2:L/1000:L/2;
%----- size of Hilbert space
nbase=50;
% wavefunction and eigenvalues
ps=@(x,n,L)sqrt(2/L)*sin(n*pi*(x./L+1/2));
En=@(n,L,m)(n*pi/L)^2*1/(2*m);

%-- External potential of the fitting system
%----- for oscillator
Vp=@(x)0.5*x.^2;
%----- for double-well
%Vp=@(x)0.5*(x.^2-2).^2;

hmt=zeros(nbase);
for i=1:nbase
    for j=i:nbase
        hmt(i,j)=En(i,L,m)*(i==j)+trapz(x,Vp(x).*ps(x,i,L).*ps(x,j,L));
        hmt(j,i)=hmt(i,j);
    end
end
[V,E]=eig(hmt)

save ../test/oscillator/hmt_HM.dat hmt -ASCII
save ../test/oscillator/hmt_HM_V.dat V -ASCII
save ../test/oscillator/hmt_HM_E.dat E -ASCII

mx=20;
hold on
for p=1:mx
    f=0*x
        for q=1:nbase
            f=f+ps(x,q,L).*V(q,p)
        end
    plot(x,f.^2+E(p,p),'r')
    text(-L/2+L/50,E(p,p)+0.1,strcat('n=',num2str(p),' E=',num2str(E(p,p))))
end
title(strcat('Density distribution and energys'))
plot(x,Vp(x),'b')
ylim([0,E(mx+1,mx+1)])
hold off
F=getframe(gcf)
imwrite(F.cdata,strcat(workpath,'density_',jobtitle,'.png')) 

clf
for p=1:mx
    f=0*x;
        for q=1:nbase
            f=f+ps(x,q,L).*V(q,p)
        end
    hold on
    plot(x,f,'r')
    title(strcat('n=',num2str(p),' E=',num2str(E(p,p))))
    hold off
     
    F=getframe(gcf)
    imwrite(F.cdata,strcat(workpath,'ps_',jobtitle,'_',num2str(p),'.png'))
    clf
end

