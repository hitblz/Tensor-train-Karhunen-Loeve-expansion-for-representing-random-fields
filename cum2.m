function  c=cum2(id,x,y)
switch id
    case 1
        sigma=1;
        b=1;
        L=1;
        c=sigma*sigma*exp(-sqrt(sum((x-y).^2,2))./(b*L));
        
    case 2
        sigma=1;
        b=1;
        L=1;
%         c=sigma*sigma*exp(-norm(x-y,2)^2/((b*L)^2));
        c=sigma*sigma*exp(-sum((x-y).^2,2)./((b*L)^2));
        
    case 3
        sigma=1;
        b=1;
        L=1;
        u=sqrt(sum((x-y).^2,2))./(b*L);
        u1=u~=0;
        u2=u(u1);
        c1=sigma*sigma*u2.*besselk(1,u2);
        u3=find(u1);
        c=sigma*sigma*ones(size(u));
        c(u3)=c1;
      
    case 4
        c=min(x,y);
        
    case 5
%         b=0.5;
%         L=1;
%         lamdaf=@(n) exp(-n/(b*L));
%         c=0;
%         N=40;
%         for k=1:N
%             c=c+lambdaf(k)*sin(x)*sin(y);    
%         end
       T=1;
       numeig=80;
       eigvals=arrayfun(@(k) 4*T*T/(pi*pi*(2*k-1)*(2*k-1)), (1:numeig)');
       eigfuns=arrayfun(@(k) @(x) sqrt(2)*sin(x/sqrt(eigvals(k))), (1:numeig)', 'UniformOutput',0);        
       c=0;
       for k=1:numeig
           c=c+eigvals(k)*eigfuns{k}(x).*eigfuns{k}(y);
       end
        
end
end