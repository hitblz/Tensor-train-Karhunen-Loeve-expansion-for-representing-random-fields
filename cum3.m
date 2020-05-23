function c=cum3(id,x,y,z)
switch id
    case 1
        c=min([x,y,z],[],2);
    case 2
        sigma=1;
        b=1;
        L=1;
%         c=sigma*sigma*exp(-norm(x-y,2)^2/((b*L)^2));
        c=sigma*sigma*sigma*exp(-sum((x-y).^2+(x-z).^2+(y-z).^2,2)./((b*L)^2));
    case 3
        sigma=1;
        b=1;
        L=1;
        cov=[2 1 1;1 2 1;1 1 2];
        c=sigma*sigma*sigma*exp(exp(-sum(([x y z]*cov).*[x y z],2)./((b*L)^2)));
        
    case 4
       T=1;
       numeig=80;
       eigvals=arrayfun(@(k) 4*T*T/(pi*pi*(2*k-1)*(2*k-1)), (1:numeig)');
       eigfuns=arrayfun(@(k) @(x) sqrt(2)*sin(x/sqrt(eigvals(k))), (1:numeig)', 'UniformOutput',0);        
       c=0;
       for k=1:numeig
           c=c+eigvals(k)*eigfuns{k}(x).*eigfuns{k}(y).*eigfuns{k}(z);
       end
       
end
end