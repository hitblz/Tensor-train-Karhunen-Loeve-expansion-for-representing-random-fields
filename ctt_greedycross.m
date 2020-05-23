function out=ctt_greedycross(d, fun, option, ndoetemp,tol2)
% continuous tensor train decomposition based on a greedy cross
% interpolation algorithm

%d--an integer no less than two, dimension of function fun
%fun=@(ind) fun(ind)

%opt: 1 add one candidate each time
%       2 add maximum number of candidates

%tol2=1e-6;
%initial pivot set: one point selected from a QMC set
qmp=sobolset(d,'Skip',randi(10000));
nqmc=10000;
doeini=qmp(1:nqmc,:);
cand=fun(doeini);
[qmax,idini]=max(abs(cand));
inipivot=doeini(idini,:);
idsets=cell(d-1,2);
for k=1:d-1
    idsets{k,1}=inipivot(1:k);%I<=k
    idsets{k,2}=inipivot(k+1:d);%I>k
end
% figure
%special case: d=2
if d==2
    %ndoetemp=100;
    maxrank=1000;
    errmax=1;
    f1=@(x,y) fun([x,y]);% for selecting candidates
    cardinterp=0;
    flag=0;
    figure
    h=animatedline;
    addpoints(h,cardinterp,errmax);
    drawnow;
    set(gca,'YScale','log');
%     while errmax>tol2
    for cardinterp=1:maxrank
%         cardinterp=cardinterp+1;
        nI1=size(idsets{1,2},1);%cardinality of I>1
        %(backup)
        %errs=zeros(n,1);
        %cancoos=zeros(n,2);
        
        % construct vector-valued function of x
        Axid=num2cell(1:nI1);
        Ax=cellfun(@(c) @(x) fun([x,repmat(idsets{1,2}(c,:),length(x),1)]), Axid, 'UniformOutput',0);
        Amidid=zeros(nI1*nI1,2);
        Amidid(:,1)=reshape(repmat([1:nI1]',[1 nI1]),nI1*nI1,1);
        Amidid(:,2)=reshape(repmat(1:nI1,[nI1 1]),nI1*nI1,1);
        Amid=reshape(fun([idsets{1,1}(Amidid(:,1),:) idsets{1,2}(Amidid(:,2),:)]),nI1,nI1);
        Ainv=pinv(Amid);
        % construct vector-valued function of y
        Ay=arrayfun(@(c) @(y) fun([repmat(idsets{1,1}(c,:),length(y),1),y]), [1:nI1]', 'UniformOutput',0);
        
        %f2=@(x,y) cellfun(@(c) feval(c,x), Ax)*Ainv*cellfun(@(c)
        %feval(c,y), Ay); %original form
        f2=@(x,y) sum((cell2mat(cellfun(@(c) feval(c,x), Ax, 'UniformOutput', false))*Ainv).*...
            cell2mat(cellfun(@(c) feval(c,y), Ay', 'UniformOutput', false)),2);
        
        ri=randi(10000);
        doetemp=qmp(ri+1:ri+ndoetemp,:);
        %         doetemp=num2cell(doetemp,2);
        %         valdoetemp=cellfun(@(c) f1(c(1),c(2))-f2(c(1),c(2)), doetemp);
        valdoetemp=f1(doetemp(:,1),doetemp(:,2))-f2(doetemp(:,1),doetemp(:,2));
        
        [errmax,id]=max(abs(valdoetemp));
        errmax=errmax/qmax;        
        if cardinterp==1
            errdm=errmax;
        else
            errdm=[errdm,errmax];
        end
        if errmax<=tol2
            break
        end
        addpoints(h,cardinterp,errmax);
        drawnow;
        xopt=doetemp(id,1);
        yopt=doetemp(id,2);
        cancoo=[xopt,yopt];
        idsetsL=union(idsets{1,1},cancoo(1),'rows');
        idsetsR=union(idsets{1,2},cancoo(2),'rows');
        if size(idsetsL,1)~=size(idsetsR,1)
            flag=1;
            break
        end
        idsets{1,1}=idsetsL;
        idsets{1,2}=idsetsR;
%         if cardinterp==maxrank
        
    end
    if errmax<=tol2
        fprintf('Dimension (1,2) converged at rank %i.\n',size(idsets{1,1},1));
    end
    fprintf('sweeping completed');
    close(gcf)
%     out1=cell(2,2);
%     nI1=size(idsets{1,2},1);
%     Ax=arrayfun(@(c) chebfun(@(x) fun([x,idsets{1,2}(c,:)]), [0 1], 'vectorize','splitting','on'), 1:nI1, ...
%         'UniformOutput',0);
%     
%     
%     Amidid=zeros(nI1*nI1,2);
%     Amidid(:,1)=reshape(repmat([1:nI1]',[1 nI1]),nI1*nI1,1);
%     Amidid(:,2)=reshape(repmat(1:nI1,[nI1 1]),nI1*nI1,1);
%     Amid=reshape(fun([idsets{1,1}(Amidid(:,1),:) idsets{1,2}(Amidid(:,2),:)]),nI1,nI1);
%     Amidinv=pinv(Amid);
%     
%     out1{1,1}=chebmatrix(Ax);
%     out1{1,2}=Amidinv;
%     
%     Ay=arrayfun(@(c) chebfun(@(y) fun([idsets{1,1}(c,:),y]), [0 1], 'vectorize','splitting','on'), [1:nI1]', ...
%         'UniformOutput',0);
%     out1{2,1}=chebmatrix(Ay);
%     out1{2,2}=1;
    
%     out={out1,idsets,errdm,flag};
out={idsets,errdm,flag};
    
elseif d>2
    
    switch option
        case 1
            %***************************************************************************
            maxhswps=1000;
            %***************************************************************************
            nhswps=0;
            errd=ones(d-1,1);
            flagd=zeros(d-1,1);
            figure
            H=arrayfun(@(x) animatedline, 1:d-1, 'UniformOutput',false);
            for k=1:d-1
                    addpoints(H{k},nhswps,errd(k));
            end
            drawnow;
            set(gca,'YScale','log');
%             while max(errd)>tol2
            while nhswps<maxhswps
                nhswps=nhswps+1;
                %left-to-right sweep:
                for k=1:d-1
                    if errd(k)<=tol2 || flagd(k)==1
                        continue
                    else
                        outISE=ISE(k,idsets,fun,d,1,ndoetemp,tol2);
                        idsets=outISE{1};
                        errd(k)=outISE{2};
                        flagd(k)=outISE{3};
                    end
                end
                if nhswps==1
                    errdm=errd;
                else
                    errdm=[errdm,errd];
                end                
                for k=1:d-1
                    addpoints(H{k},nhswps,errd(k));
                end
                drawnow;
                set(gca,'YScale','log');
                if max(errd)<=tol2 || min(flagd)==1 || nhswps==maxhswps
                    break
                end
                %right-to-left sweep:
                nhswps=nhswps+1;
                for k=d:-1:2
                    if errd(k-1)<=tol2 || flagd(k-1)==1
                        continue
                    else
                        outISE=ISE(k-1,idsets,fun,d,1,ndoetemp,tol2);
                        idsets=outISE{1};
                        errd(k-1)=outISE{2};
                        flagd(k-1)=outISE{3};
                    end
                end
                errdm=[errdm,errd];
                for k=1:d-1
                    addpoints(H{k},nhswps,errd(k));
                end
                drawnow;
                set(gca,'YScale','log');
                if max(errd)<=tol2 || min(flagd)==1 || nhswps==maxhswps
                    break
                end
            end
            
        case 2
            errdm=cell(d-1,1);
            for k=1:d-1
                outISE=ISE(k,idsets,fun,d,2);
                idsets=outISE{1};
                errdm{k}=outISE{2};
            end
    end
    fprintf('sweeping completed');
%     close(gcf)
    out={idsets,errdm,flagd};    
    %out=struct('funmat',out1,'interpsets',idsets,'abserrors',errd);
else
    error('invalid d');
end
end