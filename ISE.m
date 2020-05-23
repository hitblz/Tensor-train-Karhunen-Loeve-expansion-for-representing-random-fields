function out=ISE(k,idsets,fun,d,option,ndoetemp,tol)
%Interpolation set expansion algorithm in
%tensor train decomposition
%fun=fun(ind)

%opt: 1 add one candidate each time
%       2 add maximum number of candidates


% tol=1e-6;
% ndoetemp=100;
maxrank=100;%for option=2
p = haltonset(2,'Skip',1e3,'Leap',4);

if k==1
    errmax=1;
    %     cardinterp=0;
    flag=0;
    %     while errmax>tol
    for cardinterp=1:maxrank
        %         cardinterp=cardinterp+1;
        nI2rig=size(idsets{2,2},1);%cardinality of I>2
        nI1=size(idsets{1,2},1);%cardinality of I>1
        
        % construct vector-valued function of x
        %         Ax=arrayfun(@(c) @(x) fun([x, idsets{1,2}(c,:)]), 1:nI1, 'UniformOutput',0);
        M1=@(x) reshape(fun([reshape(x*ones(1,nI1),length(x)*nI1,1),...
            idsets{1,2}(reshape(ones(length(x),1)*(1:nI1),length(x)*nI1,1),:)]),...
            length(x),nI1);
        
        Amidid=zeros(nI1*nI1,2);
        Amidid(:,1)=reshape(repmat([1:nI1]',[1 nI1]),nI1*nI1,1);
        Amidid(:,2)=reshape(repmat(1:nI1,[nI1 1]),nI1*nI1,1);
        Amid=reshape(fun([idsets{1,1}(Amidid(:,1),:) idsets{1,2}(Amidid(:,2),:)]),nI1,nI1);
        Ainv=pinv(Amid);
        
        % construct vector-valued function of y
        %         id=zeros(nI1,nI2rig,2);
        %         id(:,:,1)=repmat([1:nI1]',[1 nI2rig]);
        %         id(:,:,2)=repmat(1:nI2rig,[nI1 1]);
        %         Ay=arrayfun(@(c1,c2) @(y) fun([idsets{1,1}(c1,:),y,idsets{2,2}(c2,:)]), id(:,:,1),id(:,:,2), 'UniformOutput',0);
        
        M2=@(y) reshape(fun([idsets{1,1}(reshape(ones(length(y),1)*(1:nI1),length(y)*nI1,1),:), ...
            reshape((y-floor(y))*ones(1,nI1),length(y)*nI1,1),idsets{2,2}(reshape(ceil(y)*ones(1,nI1),length(y)*nI1,1),:)]),...
            length(y),nI1)';
        
        f1=@(x,y) fun([x,y-floor(y),idsets{2,2}(ceil(y),:)]);% extend y
        %         f2=@(x,y) cellfun(@(c) feval(c,x), Ax)*Ainv*cellfun(@(c) feval(c,y-floor(y)), Ay(:,ceil(y)));
        
        ri=randi(10000);
        doetemp=p(ri+1:ri+ndoetemp,:);
        doetemp(:,2)=doetemp(:,2)*nI2rig;
        valdoetemp1=f1(doetemp(:,1),doetemp(:,2));
        valdoetemp2=sum(M1(doetemp(:,1)).*(Ainv*M2(doetemp(:,2)))',2);
        %         doetemp=num2cell(doetemp,2);
        %         valdoetemp2=cellfun(@(c) f2(c(1),c(2)), doetemp);
        [qmax1,~]=max(abs(valdoetemp1));
        [errmax,id]=max(abs(valdoetemp1-valdoetemp2));
        errmax=errmax/qmax1;
        if cardinterp==1
            errmaxv=errmax;%if option==1, errmaxv is equivalent to errmax
        else
            errmaxv=[errmaxv;errmax];
        end
        if errmax<=tol
            break
        end
        xopt=doetemp(id,1);
        yopt=doetemp(id,2);
        cancoo=[xopt,yopt-floor(yopt)];
        idopt=ceil(yopt);
        idsetsL=union(idsets{1,1},cancoo(1),'rows');
        idsetsR=union(idsets{1,2},[cancoo(2) idsets{2,2}(idopt,:)],'rows');
        if size(idsetsL,1)~=size(idsetsR,1)
            flag=1;
            break
        end
        idsets{1,1}=idsetsL;
        idsets{1,2}=idsetsR;
        if option==1
            break
        end
    end
    
elseif k<d-1
    errmax=1;
%     cardinterp=0;
    flag=0;
%     while errmax>tol
    for cardinterp=1:maxrank
%         cardinterp=cardinterp+1;
        nIkm1lef=size(idsets{k-1,1},1);
        nIkp1rig=size(idsets{k+1,2},1);
        nIk=size(idsets{k,1},1);
        
        %         id=zeros(nIkm1lef,nIk,2);
        %         id(:,:,1)=repmat([1:nIkm1lef]', [1 nIk]);
        %         id(:,:,2)=repmat(1:nIk, [nIkm1lef 1]);
        %         Ax=arrayfun(@(c1,c2) @(x) fun([idsets{k-1,1}(c1,:),x,idsets{k,2}(c2,:)]), id(:,:,1),id(:,:,2), 'UniformOutput',0);
        
        M1=@(x) reshape(fun([idsets{k-1,1}(reshape(ceil(x)*ones(1,nIk),length(x)*nIk,1),:), ...
            reshape((x-floor(x))*ones(1,nIk),length(x)*nIk,1),idsets{k,2}(reshape(ones(length(x),1)*(1:nIk),length(x)*nIk,1),:)]),...
            length(x),nIk);
        
        Amidid=zeros(nIk*nIk,2);
        Amidid(:,1)=reshape(repmat([1:nIk]',[1 nIk]),nIk*nIk,1);
        Amidid(:,2)=reshape(repmat(1:nIk,[nIk 1]),nIk*nIk,1);
        Amid=reshape(fun([idsets{k,1}(Amidid(:,1),:) idsets{k,2}(Amidid(:,2),:)]),nIk,nIk);
        Ainv=pinv(Amid);
        
        %         id=zeros(nIk,nIkp1rig,2);
        %         id(:,:,1)=repmat([1:nIk]', [1 nIkp1rig]);
        %         id(:,:,2)=repmat(1:nIkp1rig, [nIk 1]);
        %         Ay=arrayfun(@(c1,c2) @(y) fun([idsets{k,1}(c1,:),y,idsets{k+1,2}(c2,:)]),  id(:,:,1),id(:,:,2), 'UniformOutput',0);
        
        M2=@(y) reshape(fun([idsets{k,1}(reshape(ones(length(y),1)*(1:nIk),length(y)*nIk,1),:), ...
            reshape((y-floor(y))*ones(1,nIk),length(y)*nIk,1),idsets{k+1,2}(reshape(ceil(y)*ones(1,nIk),length(y)*nIk,1),:)]),...
            length(y),nIk)';
        
        f1=@(x,y) fun([idsets{k-1,1}(ceil(x),:), x-floor(x),y-floor(y),idsets{k+1,2}(ceil(y),:)]);
        %         f2=@(x,y) cellfun(@(c) feval(c,x-floor(x)), Ax(ceil(x),:))*Ainv*cellfun(@(c) feval(c,y-floor(y)), Ay(:,ceil(y)));
        
        ri=randi(10000);
        doetemp=p(ri+1:ri+ndoetemp,:);
        doetemp(:,1)=doetemp(:,1)*nIkm1lef;
        doetemp(:,2)=doetemp(:,2)*nIkp1rig;
        valdoetemp1=f1(doetemp(:,1),doetemp(:,2));
        ma1=M1(doetemp(:,1));
        ma2=M2(doetemp(:,2));
        valdoetemp2=sum(ma1.*(Ainv*ma2)',2);
        [qmax1,~]=max(abs(valdoetemp1));
        [errmax,id]=max(abs(valdoetemp1-valdoetemp2));
        errmax=errmax/qmax1;
        if cardinterp==1
            errmaxv=errmax;
        else
            errmaxv=[errmaxv;errmax];
        end
        if errmax<=tol
            break
        end
        %         xopt=doetemp{id}(1);
        %         yopt=doetemp{id}(2);
        xopt=doetemp(id,1);
        yopt=doetemp(id,2);
        cancoo=[xopt-floor(xopt),yopt-floor(yopt)];
        idopt=ceil([xopt,yopt]);
        idsetsL=union(idsets{k,1},[idsets{k-1,1}(idopt(1),:) cancoo(1)],'rows');
        idsetsR=union(idsets{k,2},[cancoo(2) idsets{k+1,2}(idopt(2),:)],'rows');
        if size(idsetsL,1)~=size(idsetsR,1)
            flag=1;
            break
        end
        idsets{k,1}=idsetsL;
        idsets{k,2}=idsetsR;
        if option==1
            break
        end
    end
    
else%k==d-1
    errmax=1;
%     cardinterp=0;
    flag=0;
%     while errmax>tol
    for cardinterp=1:maxrank
%         cardinterp=cardinterp+1;
        nIdm2lef=size(idsets{d-2,1},1);
        nIdm1=size(idsets{d-1,1},1);
        
        %         id=zeros(nIdm2lef,nIdm1,2);
        %         id(:,:,1)=repmat([1:nIdm2lef]', [1 nIdm1]);
        %         id(:,:,2)=repmat(1:nIdm1, [nIdm2lef 1]);
        %         Ax=arrayfun(@(c1,c2) @(x) fun([idsets{d-2,1}(c1,:),x,idsets{d-1,2}(c2,:)]), id(:,:,1),id(:,:,2), 'UniformOutput',0);
        M1=@(x) reshape(fun([idsets{d-2,1}(reshape(ceil(x)*ones(1,nIdm1),length(x)*nIdm1,1),:), ...
            reshape((x-floor(x))*ones(1,nIdm1),length(x)*nIdm1,1),idsets{d-1,2}(reshape(ones(length(x),1)*(1:nIdm1),length(x)*nIdm1,1),:)]),...
            length(x),nIdm1);
        
        Amidid(:,1)=reshape(repmat([1:nIdm1]',[1 nIdm1]),nIdm1*nIdm1,1);
        Amidid(:,2)=reshape(repmat(1:nIdm1,[nIdm1 1]),nIdm1*nIdm1,1);
        Amid=reshape(fun([idsets{d-1,1}(Amidid(:,1),:) idsets{d-1,2}(Amidid(:,2),:)]),nIdm1,nIdm1);
        Ainv=pinv(Amid);
        
        %         Ay=arrayfun(@(c) @(y) fun([idsets{d-1,1}(c,:), y]), [1:nIdm1]', 'UniformOutput',0);
        M2=@(y) reshape(fun([idsets{d-1,1}(reshape(ones(length(y),1)*(1:nIdm1),length(y)*nIdm1,1),:), ...
            reshape((y-floor(y))*ones(1,nIdm1),length(y)*nIdm1,1)]),...
            length(y),nIdm1)';
        
        f1=@(x,y)  fun([idsets{d-2,1}(ceil(x),:), x-floor(x),y]);
        %         f2=@(x,y) cellfun(@(c) feval(c,x-floor(x)), Ax(ceil(x),:))*Ainv*cellfun(@(c) feval(c,y), Ay);
        
        ri=randi(10000);
        doetemp=p(ri+1:ri+ndoetemp,:);
        doetemp(:,1)=doetemp(:,1)*nIdm2lef;
        valdoetemp1=f1(doetemp(:,1),doetemp(:,2));
        valdoetemp2=sum(M1(doetemp(:,1)).*(Ainv*M2(doetemp(:,2)))',2);
        %         doetemp=num2cell(doetemp,2);
        %         valdoetemp2=cellfun(@(c) f2(c(1),c(2)), doetemp);
        [qmax1,~]=max(abs(valdoetemp1));
        [errmax,id]=max(abs(valdoetemp1-valdoetemp2));
        errmax=errmax/qmax1;
        if cardinterp==1
            errmaxv=errmax;
        else
            errmaxv=[errmaxv;errmax];
        end
        if errmax<=tol
            break
        end
        %         xopt=doetemp{id}(1);
        %         yopt=doetemp{id}(2);
        xopt=doetemp(id,1);
        yopt=doetemp(id,2);
        cancoo=[xopt-floor(xopt),yopt];
        idopt=ceil(xopt);
        idsetsL=union(idsets{d-1,1},[idsets{d-2,1}(idopt,:) cancoo(1)],'rows');
        idsetsR=union(idsets{d-1,2},cancoo(2),'rows');
        if size(idsetsL,1)~=size(idsetsR,1)
            flag=1;
            break
        end
        idsets{d-1,1}=idsetsL;
        idsets{d-1,2}=idsetsR;
        if option==1
            break
        end
    end
end
if errmax<=tol
    fprintf('Dimension (%i, %i) converged at rank %i.\n',k,k+1,size(idsets{k,1},1));
end
out={idsets errmaxv flag};
end