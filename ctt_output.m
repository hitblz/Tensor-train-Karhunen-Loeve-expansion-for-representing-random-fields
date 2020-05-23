function out=ctt_output(d,fun, idsets)
if d==2
    out=cell(2,2);
    nI1=size(idsets{1,2},1);
    Ax=arrayfun(@(c) chebfun(@(x) fun([x,idsets{1,2}(c,:)]), [0 1], 'vectorize','splitting','on'), 1:nI1, ...
        'UniformOutput',0);
    
    Amidid=zeros(nI1*nI1,2);
    Amidid(:,1)=reshape(repmat((1:nI1)',[1 nI1]),nI1*nI1,1);
    Amidid(:,2)=reshape(repmat(1:nI1,[nI1 1]),nI1*nI1,1);
    Amid=reshape(fun([idsets{1,1}(Amidid(:,1),:) idsets{1,2}(Amidid(:,2),:)]),nI1,nI1);
    Amidinv=pinv(Amid);
    
    out{1,1}=chebmatrix(Ax)*Amidinv;
    out{1,2}=Amidinv;
    
    Ay=arrayfun(@(c) chebfun(@(y) fun([idsets{1,1}(c,:),y]), [0 1], 'vectorize','splitting','on'), (1:nI1)', ...
        'UniformOutput',0);
    out{2,1}=chebmatrix(Ay);
    out{2,2}=1;
elseif d>2
    out=cell(d,2);
    %for output
    for k=1:d
        if k==1
            nI1=size(idsets{1,2},1);
            Ax=arrayfun(@(c) chebfun(@(x) fun([x,idsets{1,2}(c,:)]), [0 1], 'vectorize','splitting','on'), 1:nI1, ...
                'UniformOutput',0);
            
            Amidid=zeros(nI1*nI1,2);
            Amidid(:,1)=reshape(repmat((1:nI1)',[1 nI1]),nI1*nI1,1);
            Amidid(:,2)=reshape(repmat(1:nI1,[nI1 1]),nI1*nI1,1);
            Amid=reshape(fun([idsets{k,1}(Amidid(:,1),:) idsets{k,2}(Amidid(:,2),:)]),nI1,nI1);
            Amidinv=pinv(Amid);
            
            out{k,1}=chebmatrix(Ax)*Amidinv;
            out{k,2}=Amidinv;
            
        elseif k<d
            nIkm1=size(idsets{k-1,1},1);
            nIk=size(idsets{k,2},1);
            id=zeros(nIkm1,nIk,2);
            id(:,:,1)=repmat((1:nIkm1)',[1 nIk]);
            id(:,:,2)=repmat(1:nIk,[nIkm1 1]);
            if all([nIkm1 nIk]<=100)
            Ax=arrayfun(@(c1,c2) chebfun(@(x) fun([idsets{k-1,1}(c1,:) x idsets{k,2}(c2,:)]), [0 1], 'vectorize','splitting','on'),id(:,:,1), id(:,:,2), ...
                'UniformOutput',0);
            else
            Ax=cell(nIkm1,nIk);
            numele=nIkm1*nIk;
            id1=id(:,:,1);
            id2=id(:,:,2);
            parfor kk=1:numele
                Ax{kk}=chebfun(@(x) fun([idsets{k-1,1}(id1(kk),:) x idsets{k,2}(id2(kk),:)]), [0 1], 'vectorize','splitting','on');
            end
            end
            
            Amidid=zeros(nIk*nIk,2);
            Amidid(:,1)=reshape(repmat((1:nIk)',[1 nIk]),nIk*nIk,1);
            Amidid(:,2)=reshape(repmat(1:nIk,[nIk 1]),nIk*nIk,1);
            Amid=reshape(fun([idsets{k,1}(Amidid(:,1),:) idsets{k,2}(Amidid(:,2),:)]),nIk,nIk);
            Amidinv=pinv(Amid);
            
            out{k,1}=chebmatrix(Ax)*Amidinv;
            out{k,2}=Amidinv;
        else
            nIdm1=size(idsets{d-1,1},1);
            Ax=arrayfun(@(c) chebfun(@(x) fun([idsets{d-1,1}(c,:) x]), [0 1], 'vectorize','splitting','on'), (1:nIdm1)', ...
                'UniformOutput',0);
            
            out{k,1}=chebmatrix(Ax);
            out{k,2}=1;
        end
    end
else
    error('invalid d');
end
end