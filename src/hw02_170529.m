cd d:\\codefile\\mldatafile
m=20.1797*1822.39/2%(ʹ��Լ��������ԭ�ӵ�λ�����Ե�������Ϊ1)
eg=35.6/315775%(0.00307eV=296J/mol,ʹ��ԭ�ӵ�λ��)
s=274.9/52.92%(lenard jones ��ƽ����룬5.195)
bs=2^(1/6)*s%(lenard jones ����׵ļ��Ϊ2^(1/6)s=1.122s)
%������LJ�Ƶ���ز���

a=10*s;%a�ǣ������������
b=0.5*s;%b�ǲ����������
% x=-a:a/1000:b;%���ֳ����λ���
% y=b:a/1000:a;
% z=-a:a/1000:a;
L=2*a;%(�߽�ѡ���е���֣�)

nbase=500;
ps=@(x,n,L)sqrt(2/L)*sin(n*pi*(x./L+1/2));
VLJ=@(x)4*eg*((s./(bs-x)).^12-(s./(bs-x)).^6);%(������ײ�Ϊ�������)


hmt=zeros(nbase);
for i=1:nbase
    for j=i:nbase
        hmt(i,j)=(i*pi/L)^2/(2*m)*(i==j)+integral(@(x)VLJ(x).*ps(x,i,L).*ps(x,j,L),-a,b)+integral(@(x)VLJ(b).*ps(x,i,L).*ps(x,j,L),b,a);
        hmt(j,i)=hmt(i,j);
    end
end
[V,E]=eig(hmt)

cd d:\\codefile\\mldatafile
save hmt_LJ.dat hmt -ASCII
save hmt_LJ_V.dat V -ASCII
save hmt_LJ_V.mat V -ASCII
save hmt_LJ_E.dat E -ASCII


hold on
for p=1:20
     f=@(x)0;
     for q=1:nbase
         f=@(x)(f(x)+ps(x,q,a).*V(q,p))
     end
     plot(x,0.00001*f.^2+E(p,p),'r')%��������
     text(-0.9*s,E(p,p)+0.000001,strcat('n=',num2str(p),' E=',num2str(E(p,p))))
end
title(strcat('LJ�ơ���ǰ10����ϲ�����ʾ��ͼ������������Ϊ׼��'))
plot(x,VLJ(x),'b')
xlim([-3*s,3*s])
ylim([0,0.00020])
hold off

F=getframe(gcf)
imwrite(F.cdata,strcat('psof_LJ_tot','.png')) 

clf
for p=1:5
   f=0*x;
     for q=1:nbase
         f=f+ps(x,q,a).*V(q,p)
     end
     hold on
     plot(x,f.^2,'r')%��������
     title(strcat('n=',num2str(p),' E=',num2str(E(p,p))))
     xlim([-6*s,6*s])
     hold off
     %text(-0.9*s,E(p,p)+0.000001,strcat('n=',num2str(p),' E=',num2str(E(p,p))))
     F=getframe(gcf)
     imwrite(F.cdata,strcat('Rof_LJ_',num2str(p),'.png'))
     clf
end