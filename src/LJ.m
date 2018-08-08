%Lennard-Jones system
clear;clf;

% define the task
workpath='../test/LJ/'
jobtitle='LJ'

m=20.1797*1822.39/2
s=274.9/52.92%(lenard jones bottom position, 5.195)
bs=2^(1/6)*s%(lenard jones balance position, bs=1.122s)
a=0.5*s;   % left section
b=10*s;% right section
L=a+b
x=-a:L/1000:b;
nbase=100;

ps=@(x,n,L)sqrt(2/L)*sin(n*pi*(x./L+a/L));

%-- parameters of L-J potential
eg=35.6/315775%(0.00307eV=296J/mol)
VLJ=@(x)4*eg*((s./(bs+x)).^12-(s./(bs+x)).^6);


hmt=zeros(nbase);
for i=1:nbase
    for j=i:nbase
        hmt(i,j)=(i*pi/L)^2/(2*m)*(i==j)+integral(@(x)VLJ(x).*ps(x,i,L).*ps(x,j,L),-a,b);
        hmt(j,i)=hmt(i,j);
    end
end
[V,E]=eig(hmt)

save ../test/LJ/hmt_HM.dat hmt -ASCII
save ../test/LJ/hmt_HM_V.dat V -ASCII
save ../test/LJ/hmt_HM_E.dat E -ASCII

scale=0.00001
deltaE=0.000001
hold on
for p=1:10
    f=0*x;
    for q=1:nbase
        f=f+ps(x,q,L).*V(q,p)
    end
    plot(x,scale*f.^2+E(p,p),'r')
    %text(-0.9*s,E(p,p)+deltaE,strcat('n=',num2str(p),' E=',num2str(E(p,p))))
end
title(strcat('Density distribution'))
plot(x,VLJ(x),'b')
xlim([-0.5*a,0.5*b])
ylim([-0.0001,0.00002])
hold off
F=getframe(gcf)
imwrite(F.cdata,strcat(workpath,'psof_LJ_tot','.png')) 

clf
scale=0.00005
for p=1:5
    f=0*x;
    for q=1:nbase
        f=f+ps(x,q,L).*V(q,p)
    end
    hold on
    plot(x,VLJ(x),'b')
    plot(x,scale*f.^2+E(p,p),'r')
    title(strcat('n=',num2str(p),' E=',num2str(E(p,p))))
    xlim([-a,3*s])
    ylim([-0.0001,0.00002])
    hold off
    F=getframe(gcf)
    imwrite(F.cdata,strcat(workpath,'PS2_',jobtitle,'_',num2str(p),'.png'))
    clf
end

