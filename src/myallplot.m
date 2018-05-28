%myplot.m
function []=myallplot(mx,nbase,V,E,ps,Vp,L,x1,x2)
x=x1:(x2-x1)/1000:x2
hold on
for p=1:mx
     f=@(x)0;
     for q=1:nbase
         f=@(x)(f(x)+ps(x,q,L).*V(q,p))
     end
     plot(x,0.00001*f(x).^2+E(p,p),'r')%��������
     text(0.99*x1+0.01*x2,E(p,p)+0.000001,strcat('n=',num2str(p),' E=',num2str(E(p,p))))
end
title(strcat('LJ�ơ���ǰ10����ϲ�����ʾ��ͼ������������Ϊ׼��'))
plot(x,Vp(x),'b')
xlim([x1,x2])
ylim([-0.0001,0.00005])
hold off
F=getframe(gcf)
imwrite(F.cdata,strcat('psof_LJ_tot','.png')) 