
add='E:/my papers/Random field representation/v3/figures/';
dimpar=1;% 1 2 3 ***
dimphy=1; %1 2 3 and >=dimpar
cum2id=4;
cum3id=4;
switch dimpar
    case 1
        switch dimphy
            case 1
                igtid=4;% ID of isogeometric transform
                nrb=igf1d(igtid); % construct isogeometric transform
                x0=@(xi1) nrbeval(nrb,xi1);
                x1=@(u) u(1,:)';% a row vector
                x=@(xi1) x1(x0(xi1));
            case 2
                igtid=3;% ID of isogeometric transform
                nrb=igf1d(igtid); % construct isogeometric transform
                x0=@(xi1) nrbeval(nrb,xi1);
                x1=@(u) u(1:2,:)';% points listed along the second dimension
                x=@(xi1) x1(x0(xi1));
            case 3
                igtid=1;% ID of isogeometric transform
                nrb=igf1d(igtid); % construct isogeometric transform
                x=@(xi1) nrbeval(nrb,xi1)';
        end
         T=1;
         numeig=80;
         eigvals=arrayfun(@(k) 4*T*T/(pi*pi*(2*k-1)*(2*k-1)), (1:numeig)');
         eigfuns=arrayfun(@(k) @(x) sqrt(2)*sin(x/sqrt(eigvals(k))), (1:numeig)', 'UniformOutput',0);
       %FEM 
        C2t=@(ind) cum2(cum2id, x(ind(:,1)), x(ind(:,2)));
        %eigfem(d,n,fun,type,neig,m)
%         out=eigfem(1,100,C2t,1,20,3);
%         modefem=out{1};
%         eigvalfem=out{2};


        %method 2 (for comparison)
          
        %task 1: assess the accuracy of the reconstructed covariance
        %function with an experiment
        %out=ctt_greedycross(d, fun, option, ndoetemp,tol2)
        tic
        cttout=ctt_greedycross(2, C2t, 1, 400, 1e-6);
        toc
        interpsets1=cttout{1};
        cttout2=ctt_output(2,C2t,interpsets1);
        toc
        C2tcell=cttout2(:,1);
        xgrid=0:0.01:1;
        ygrid=0:0.01:1;
        Z1=permute(fevalcbm(C2tcell{1},xgrid'),[3 2 1])*permute(fevalcbm(C2tcell{2},ygrid'),[1 3 2]);
        [X,Y]=meshgrid(xgrid',ygrid);
        Z0=reshape(C2t([X(:),Y(:)]),length(xgrid),length(xgrid));
        
        figure;
        scatter(interpsets1{1},interpsets1{2});box on
        xlabel('\it \xi');
        ylabel('\it \xi`');
        set(gca,'FontName', 'Times New Roman', 'Fontsize', 14,'SortMethod','childorder');
        set(gcf,'color','w');
        set(gcf,'Units','centimeters');
        set(gcf,'Position',[1.5 1.5 12 10]);
        set(gcf,'PaperUnits','centimeters');
        saveas(gcf,'interpsets1.fig');
        export_fig(gcf, [add 'interpsets1.pdf'], '-r2000', '-transparent');
        
        figure;
        surf(X,Y,Z1, 'LineWidth',0.05);
        xlabel('\it \xi');
        ylabel('\it \eta');
        zlabel('covariance');
        set(gca,'FontName', 'Times New Roman', 'Fontsize', 14,'SortMethod','childorder');
        set(gcf,'color','w');
        set(gcf,'Units','centimeters');
        set(gcf,'Position',[1.5 1.5 12 10]);
        set(gcf,'PaperUnits','centimeters');
        saveas(gcf,'reconcov.fig');
        export_fig(gcf, [add 'reconcov.pdf'], '-r2000', '-transparent');
        
        figure;
        surf(X,Y,Z0-Z1,'LineWidth',0.05);
        set(gca,'FontName', 'Times New Roman', 'Fontsize', 14,'SortMethod','childorder');
        xlabel('\it \xi');
        ylabel('\it \eta');
        zlabel('errors');
        set(gcf,'color','w');
        set(gcf,'Units','centimeters');
        set(gcf,'Position',[1.5 1.5 12 10]);
        set(gcf,'PaperUnits','centimeters');
        saveas(gcf,'reconcoverr.fig');
        export_fig(gcf, [add 'reconcoverr.pdf'], '-r2000', '-transparent');
        
        %task2: global random test of errors, iteration error
        %investigations, rank variability
        ndoetemps=[50 100 200 400];
        ntest=100;
        ndoe=1000;
        tol=1e-6;
        testerrv=zeros(ntest,4);
        ranks=zeros(ntest,4);
        ctt1errc=cell(ntest,4);
        for nn=1:4
            for k=1:ntest
%                 ctt_greedycross(d, fun, option, ndoetemp,tol2)
                cttout=ctt_greedycross(2, C2t, 1, ndoetemps(nn),tol);
                interpsets=cttout{1};
                cttout2=ctt_output(2,C2t,interpsets);
                C2tcell=cttout2(:,1);                              
                doe=rand(ndoe,2);
                refv=C2t(doe);
                testerrv(k,nn)=sqrt(mean((sum(permute(fevalcbm(C2tcell{1},doe(:,1)),[3 2 1]).*permute...
                    (fevalcbm(C2tcell{2},doe(:,2)),[3 1 2]),2)-refv).^2)/mean(refv.^2));
                ranks(k,nn)=size(C2tcell{1},2);
                toc
            end
        end
        
        figure
        boxplot(testerrv,'Labels',{'50','100','200','400'});
        set(gca,'FontName', 'Times New Roman', 'Fontsize', 14,'YScale', 'log');
        xlabel('\it m_{\rm 1}');
        ylabel('\it \epsilon_{\rm g}');
        set(gcf,'color','w','Units','centimeters','Position',[1.5 1.5 12 10],'PaperUnits','centimeters');
        saveas(gcf,'epg_1.fig');
        export_fig(gcf, [add 'epg_1.pdf'], '-r2000', '-transparent');
        
        figure
        boxplot(ranks,'Labels',{'50','100','200','400'});
        set(gca,'FontName', 'Times New Roman', 'Fontsize', 14,'SortMethod','childorder');
        xlabel('\it m_1');
        ylabel('rank');
        set(gcf,'color','w');
        set(gcf,'Units','centimeters');
        set(gcf,'Position',[1.5 1.5 12 10]);
        set(gcf,'PaperUnits','centimeters');
        saveas(gcf,'rank_1.fig');
        export_fig(gcf, [add 'rank_1.pdf'], '-r2000', '-transparent');
        
        tols= [1e-4 1e-5];
        testerrv2=zeros(ntest,3);
        ranks2_=zeros(ntest,3);
        testerrv2(:,3)=testerrv(:,2);
        ranks2_(:,3)=ranks(:,2);
        for nn=1:2
            for k=1:ntest
                cttout=ctt_greedycross(2, C2t, 1, 100,tols(nn));
                interpsets=cttout{1};
                cttout2=ctt_output(2,C2t,interpsets);
                C2tcell=cttout2(:,1);                
                doe=rand(ndoe,2);
                refv=C2t(doe);
                testerrv2(k,nn)=sqrt(mean((sum(permute(fevalcbm(C2tcell{1},doe(:,1)),[3 2 1]).*permute...
                    (fevalcbm(C2tcell{2},doe(:,2)),[3 1 2]),2)-refv).^2)/mean(refv.^2));
                ranks2_(k,nn)=size(C2tcell{1},2);
                toc
            end
        end
        
        figure
        boxplot(testerrv2,'Labels',{10^{-4},10^{-5},10^{-6}});
        set(gca,'FontName', 'Times New Roman', 'Fontsize', 14,'SortMethod','childorder','YScale', 'log');
        xlabel('\it tol');
        ylabel('\it \epsilon_{\rm g}');
        set(gcf,'color','w','Units','centimeters','Position',[1.5 1.5 12 10],'PaperUnits','centimeters');
        saveas(gcf,'epg_2.fig');
        export_fig(gcf, [add 'epg_2.pdf'], '-r2000', '-transparent');
        
        figure
        boxplot(ranks2_,'Labels',{'10^{-4}','10^{-5}','10^{-6}'});
        set(gca,'FontName', 'Times New Roman', 'Fontsize', 14,'SortMethod','childorder');
        xlabel('\it tol');
        ylabel('rank');
        set(gcf,'color','w','Units','centimeters','Position',[1.5 1.5 12 10],'PaperUnits','centimeters');
        saveas(gcf,'rank_2.fig');
        export_fig(gcf, [add 'rank_2.pdf'], '-r2000', '-transparent');
        
        
        %task3: modes computation
        n3=3;
        S2ac=cell(1,n3);
        mode1ac=cell(1,n3);
        ctterr3c=cell(1,n3);
        for k=1:n3
            cttout=ctt_greedycross(2, C2t, 1, 400, 1e-6);
            interpsets=cttout{1};
            cttout2=ctt_output(2,C2t,interpsets);
            C2tcell=cttout2(:,1);
            clear cttout
            [mode1a,S2a]=svdre3(C2tcell{1},C2tcell{2});
            mode1ac{k}=mode1a;%*****************x-mode obtained**********************
            S2ac{k}=diag(S2a);
            toc
        end
        
        refids=[1 10 80];
        lrefid=length(refids);
        for k=1:lrefid
            eigfunref=chebfun(eigfuns{refids(k)},[0 1],'vectorize','splitting','on');
            figure
            plot(eigfunref);hold on
            t=mode1ac{1}(:,refids(k));
            t=t{1};
            u=[sum((eigfunref-t).^2), sum((eigfunref+t).^2)];
            [~,s]=min(u);
            if s==1
                sign=1;
            else
                sign=-1;
            end
            plot(sign*t,'--');hold on
            u=[sum((eigfunref-mode1ac{2}(:,refids(k))).^2), sum((eigfunref+mode1ac{2}(:,refids(k))).^2)];
            [~,s]=min(u);
            if s==1
                sign=1;
            else
                sign=-1;
            end
            t=mode1ac{2}(:,refids(k));
            t=t{1};
            plot(sign*t,':');hold on
            t=mode1ac{3}(:,refids(k));
            t=t{1};
            u=[sum((eigfunref-t).^2), sum((eigfunref+t).^2)];
            [~,s]=min(u);
            if s==1
                sign=1;
            else
                sign=-1;
            end
            plot(sign*t,'-.');hold on
            legend('exact','Location','Southeast');
            set(gca,'FontName', 'Times New Roman', 'Fontsize', 14,'SortMethod','childorder');
            xlabel('\it \xi');
            ylabel(['\it f_{\rm ' num2str(refids(k)) '}(\xi)']);
            set(gcf,'color','w');
            set(gcf,'Units','centimeters');
            set(gcf,'Position',[1.5 1.5 12 10]);
            set(gcf,'PaperUnits','centimeters');
            saveas(gcf,['m1xmode_{' num2str(refids(k)) '}.fig']);
            export_fig(gcf, [add 'm1xmode_{' num2str(refids(k)) '}.pdf'], '-r2000', '-transparent');
        end
        figure
        numeig2=80;
        scatter(1:numeig2, eigvals(1:numeig2)); hold on
        scatter(1:numeig2, S2ac{1}(1:numeig2),'*'); hold on
        scatter(1:numeig2, S2ac{2}(1:numeig2),'x'); hold on
        scatter(1:numeig2, S2ac{3}(1:numeig2),'^'); hold on
        box on
        legend('exact','TT-K-L','Location','Northeast');
        set(gca,'FontName', 'Times New Roman', 'Fontsize', 14,'SortMethod','childorder');
        xlabel('order');
        ylabel('eigenvalues');
        set(gcf,'color','w');
        set(gcf,'Units','centimeters');
        set(gcf,'Position',[1.5 1.5 12 10]);
        set(gcf,'PaperUnits','centimeters');
        saveas(gcf,'eigvalx.fig');
        export_fig(gcf, [add 'eigvalx.pdf'], '-r2000', '-transparent');
        
    case 2
        switch dimphy
            case 2
                igtid=1;% ID of isogeometric transform ***
                nrb=igf2d(igtid); % construct isogeometric transform
                x0=@(xi1,xi2) nrbeval(nrb,[xi1,xi2]);% xi1 is a 2d
                x1=@(u) u(1:2,:)';
                x=@(xi1,xi2) x1(x0(xi1,xi2));
            case 3
                igtid=2;% ID of isogeometric transform ***
                nrb=igf2d(igtid); % construct isogeometric transform
                x=@(xi1,xi2) nrbeval(nrb,[xi1,xi2])';% xi1 is a 2d
        end
        xr=@(xi1,xi2) x(xi2,xi1);
        %low-rank approximation of cumulant function
        C2t=@(ind) cum2(cum2id, x(ind(:,1),ind(:,2)), xr(ind(:,3),ind(:,4)));
        C2tref=@(ind) cum2(cum2id, x(ind(:,1),ind(:,2)), x(ind(:,3),ind(:,4)));
        %         C2tref2=@(ind) exp(-abs(ind(:,1)-ind(:,3))-abs(ind(:,2)-ind(:,4)));
        %B2t=@(x,y,y',x') C2t(x,y,x',y')
        %B2t=@(ind) C2t([ind(:,1),ind(:,2),ind(:,4),ind(:,3)]);
        tic
        cttout1=ctt_greedycross(4, C2t, 1,800,1e-6);
        interpsets2=cttout1{1};
        cttout2=ctt_output(4,C2t,interpsets2);
        cum2t1=toc;
        C2tcell=cttout2(:,1);
        clear cttout2 
        ndoe=1000;
        doe=rand(ndoe,4);
        refv=C2t(doe);
        refv2=C2tref(doe(:,[1 2 4 3]));
        diffv=zeros(ndoe,1);
        diffv2=zeros(ndoe,1);
        cdata={fevalcbm(C2tcell{1},doe(:,1)),fevalcbm(C2tcell{2},doe(:,2)),fevalcbm(C2tcell{3},doe(:,3)), fevalcbm(C2tcell{4},doe(:,4))};
        for k=1:ndoe
            diffv(k)=cdata{1}(:,:,k)*cdata{2}(:,:,k)*cdata{3}(:,:,k)*cdata{4}(:,:,k)-refv(k);
            diffv2(k)=cdata{1}(:,:,k)*cdata{2}(:,:,k)*cdata{3}(:,:,k)*cdata{4}(:,:,k)-refv2(k);
        end
        clear cdata
        epg2m2=sqrt(mean(diffv.^2)/mean(refv.^2));
        epg2m2v2=sqrt(mean(diffv2.^2)/mean(refv2.^2));
        ranks2=[size(C2tcell{2},1),size(C2tcell{3},1),size(C2tcell{4},1)];
        fprintf('Global random test completed');
        tim2=toc;
        cum2t2=tim2-cum2t1;
        
        % for eigen decomposition
        xm1=C2tcell{1};
        xm2=sum(C2tcell{2}*C2tcell{3})*C2tcell{4};%convert to Inf*n chebfun
        [mode1a,S2a]=svdre3(xm1,xm2);
        fprintf('1-mode obtained');
                        
        ym1=(mode1a'*C2tcell{1})*C2tcell{2};
        ym2=C2tcell{3}*sum(C2tcell{4}*mode1a);
        [mode2a,S2b]=svdre3(ym1,ym2);%*****************y-mode obtained**********************
        fprintf('2-mode obtained');
        tim3=toc;
        cum2t3=tim3-cum2t2;
        %****************************************************************************
        %task 1: compare with FEM
        %out=eigfem(n,fun,type,neig,m)
        neig=length(diag(S2b));
        ngrid=40;
        out=eigfem(ngrid,C2tref,1,neig,3);
        modefem=out{1};
        eigvalfem=out{2};
        
        figure
        scatter(1:neig, diag(S2b(1:neig,1:neig))); hold on
        scatter(1:neig, diag(eigvalfem),'*'); hold on
        box on
        legend('TT-K-L','FEM','Location','Northeast');
        set(gca,'FontName', 'Times New Roman', 'Fontsize', 14,'SortMethod','childorder','YScale','log');
        xlabel('number');
        ylabel('eigenvalues');
        set(gcf,'color','w');
        set(gcf,'Units','centimeters');
        set(gcf,'Position',[1.5 1.5 13.5 10]);
        set(gcf,'PaperUnits','centimeters');
        saveas(gcf,'eigval2.fig');
        export_fig(gcf, [add 'eigval2.pdf'], '-r2000', '-transparent');
        
        modesel=[1 2 4];
        nmsel=length(modesel);
        ut=0:1/ngrid:1;
        ut=ut(:);
        cdata=cell(1,nmsel+1);
        cdata{1}=fevalcbm(mode1a,ut);
        for k=1:nmsel
            cdata{k+1}=fevalcbm(mode2a(:,modesel(k)),ut);
        end
        [id1,id2]=meshgrid(1:ngrid+1);
        fungrid=cell(1,nmsel);
        for k=1:nmsel
            fungrid{k}=arrayfun(@(c1,c2) cdata{1}(:,:,c1)*cdata{k+1}(:,:,c2),id1,id2);% f(ut(id1))*g(ut(id2))
        end
        utr=1:-1/ngrid:0;
        utr=utr(:);
        [Xr,Yr]=meshgrid(ut,utr);
        for k=1:nmsel
            figure
            surf(Xr,Yr,modefem(:,:,modesel(k)),'LineStyle','none');hold on
            switch k
                case 1
                    signk=1;
                case 2
                    signk=-1;
                case 3
                    signk=-1;
            end
            mesh(ut(id1),ut(id2),signk*fungrid{k},'FaceAlpha',0.5);
            legend('FEM','TT-K-L');
            set(gca,'FontName', 'Times New Roman', 'Fontsize', 14,'SortMethod','childorder');
            xlabel('\it \xi');
            ylabel('\it \eta');
            set(gcf,'color','w');
            set(gcf,'Units','centimeters');
            set(gcf,'Position',[1.5 1.5 14 10]);
            set(gcf,'PaperUnits','centimeters','PaperSize',[13.5 10]);
            saveas(gcf,['m2xymode_' num2str(k) '.fig']);
            print(gcf, '-dpdf',[add 'm2xymode_' num2str(k) '.pdf' ], '-r2000');
            
        end
        
        
    case 3%in this case, dimphy==3
        igtid=2;% ID of isogeometric transform
        nrb=igf3d(igtid); % construct isogeometric transform
        x=@(xi1, xi2, xi3) nrbeval(nrb,[xi1, xi2, xi3])';
        xr=@(xi1,xi2,xi3) x(xi3,xi2,xi1);
        %C2=@(xi1, xi2, xi3, xi1_, xi2_, xi3_) cum2_2(x(xi1, xi2, xi3), x(xi1_, xi2_, xi3_));%***
        
        %low-rank approximation of cumulant function
        C2t=@(ind) cum2(cum2id, x(ind(:,1),ind(:,2),ind(:,3)), xr(ind(:,4),ind(:,5),ind(:,6)));
        C2tref=@(ind) cum2(cum2id, x(ind(:,1),ind(:,2),ind(:,3)), x(ind(:,4),ind(:,5),ind(:,6)));
        %B2t=@(x,y,z,z',y',x') C2t(x,y,z,x',y',z')
        %         B2t=@(ind) C2t([ind(:,1),ind(:,2),ind(:,3),ind(:,6),ind(:,5),ind(:,4)]);
        tic;
        cttout1=ctt_greedycross(6, C2t, 1,800,1e-4);
        interpsets2=cttout1{1};
        cttout2=ctt_output(6,C2t,interpsets2);
        cum2t1=toc;
        C2tcell=cttout2(:,1);
        clear cttout2
        ndoe=1000;
        doe=rand(ndoe,6);
        refv=C2t(doe);
        refv2=C2tref(doe(:,[1 2 3 6 5 4]));
        cdata={fevalcbm(C2tcell{1},doe(:,1)),fevalcbm(C2tcell{2},doe(:,2)),fevalcbm(C2tcell{3},doe(:,3)), fevalcbm(C2tcell{4},doe(:,4)),fevalcbm(C2tcell{5},doe(:,5)),fevalcbm(C2tcell{6},doe(:,6))};
        diffv=arrayfun(@(k) cdata{1}(:,:,k)*cdata{2}(:,:,k)*cdata{3}(:,:,k)*cdata{4}(:,:,k)*cdata{5}(:,:,k)*cdata{6}(:,:,k)-refv(k),(1:ndoe)');
        diffv2=arrayfun(@(k) cdata{1}(:,:,k)*cdata{2}(:,:,k)*cdata{3}(:,:,k)*cdata{4}(:,:,k)*cdata{5}(:,:,k)*cdata{6}(:,:,k)-refv2(k),(1:ndoe)');
        clear cdata
        epg2m3=sqrt(mean(diffv.^2)/mean(refv.^2));
        epg2m3v2=sqrt(mean(diffv2.^2)/mean(refv2.^2));
        if epg2m3>0.05
            pause
        end
        ranks2=[size(C2tcell{2},1),size(C2tcell{3},1),size(C2tcell{4},1),size(C2tcell{5},1),size(C2tcell{6},1)];
        fprintf('Global random test completed');
        tim2=toc;
        cum2t2=tim2-cum2t1;
        
        % for eigen decomposition
        zint=sum(mul(C2tcell{3},C2tcell{4}));
        xm1=C2tcell{1}*sum(C2tcell{2}*(zint*C2tcell{5}));
        xm2=C2tcell{6};
        %reconstruct the two chebmatrices for fast computations
        %         xm1=cellfun(@(f) chebfun(@(x) feval(f,x),[0,1],'vectorize','splitting','on'),xm1);
        %         xm2=cellfun(@(f) chebfun(@(x) feval(f,x),[0,1],'vectorize','splitting','on'),xm2);
        [mode1a,S2a]=svdre3(xm1,xm2);%*****************x-mode obtained**********************
        fprintf('x-mode obtained');
        
        m1a=mode1a'*C2tcell{1};
        ym1=m1a*C2tcell{2}*zint;
        m2a=sum(C2tcell{6}*mode1a);
        ym2=C2tcell{5}*m2a;
        ym1=cellfun(@(f) chebfun(@(x) feval(f,x),[0,1],'vectorize','splitting','on'),ym1);
        ym2=cellfun(@(f) chebfun(@(x) feval(f,x),[0,1],'vectorize','splitting','on'),ym2);
        
        [mode2a,S2b]=svdre3(ym1,ym2);%*****************y-mode obtained**********************
        fprintf('y-mode obtained');
        
        zm1=(mode2a'*(m1a*C2tcell{2}))*C2tcell{3};
        zm2=C2tcell{4}*sum(C2tcell{5}*m2a*mode2a);
        zm1r=cellfun(@(f) chebfun(@(x) feval(f,x),[0,1],'vectorize','splitting','on'),zm1);
        zm2r=cellfun(@(f) chebfun(@(x) feval(f,x),[0,1],'vectorize','splitting','on'),zm2);
        [mode3a,S2c]=svdre3(zm1,zm2);%*****************z-mode obtained**********************
        fprintf('z-mode obtained');
        tim3=toc;
        cum2t3=tim3-cum2t2;
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% third order cumulants of second-order factors
switch dimpar
    case 1        
        % third order cumulant of random factors
        C3t=@(ind) cum3(cum3id, x(ind(:,1)), x(ind(:,2)), x(ind(:,3)));
        tic;
        cttout1=ctt_greedycross(3, C3t, 1,800,1e-5);        
        interpsets3=cttout1{1};
        cttout2=ctt_output(3,C3t,interpsets3);
        cum3t1=toc;
        C3tcell=cttout2(:,1);
        clear cttout2
        ndoe=1000;
        doe=rand(ndoe,3);
        refv=C3t(doe);
        cdata={fevalcbm(C3tcell{1},doe(:,1)),fevalcbm(C3tcell{2},doe(:,2)),fevalcbm(C3tcell{3},doe(:,3))};
        diffv=arrayfun(@(k) cdata{1}(:,:,k)*cdata{2}(:,:,k)*cdata{3}(:,:,k)-refv(k), (1:ndoe)');
        epg3m1=sqrt(mean(diffv.^2)/mean(refv.^2));
        ranks3=size(C3tcell{2});
        tim2=toc;
        cum3t2=tim2-cum3t1;
        
        tic;
        cttout1=ctt_greedycross(2, C2t, 1, 400,1e-6);
        interpsets2=cttout1{1};
        cttout2=ctt_output(2,C2t,interpsets2);
        cum2t1=toc;
        C2tcell=cttout2(:,1);
        clear cttout2
        ndoe=1000;
        doe=rand(ndoe,2);
        refv=C2t(doe);
        cdata={fevalcbm(C2tcell{1},doe(:,1)),fevalcbm(C2tcell{2},doe(:,2))};
        diffv=arrayfun(@(k) cdata{1}(:,:,k)*cdata{2}(:,:,k)-refv(k), (1:ndoe)');
        epg2m1=sqrt(mean(diffv.^2)/mean(refv.^2));
        ranks2=size(C2tcell{1},2);
        tim3=toc;
        cum2t2=tim3-cum2t1;
        [mode1a,S2a]=svdre3(C2tcell{1},C2tcell{2});
        tim4=toc;
        cum2t3=tim4-tim3;
        
        tic;    
        ni=size(mode1a,2);
        cum3cell=cell(1,3);
        for k=1:3
            sz1=size(C3tcell{k},1);
            sz2=size(C3tcell{k},2);
            ct=C3tcell{k};%chebmatrix
            mt=zeros(ni,sz2,sz1);
            for k2=1:sz1
                mt(:,:,k2) = mul(mode1a',ct(k2,:),'par');
                fprintf('%d  %d', k,k2);
                toc
            end
            cum3cell{k}=permute(mt,[3 1 2]);% in accordance with literature
        end
        cum3t3=toc;
%         for k=1:3
%             bk=cum3cell{k};
%             for i=1:ni
%                 bk(:,i,:)=bk(:,i,:)*sqrt(S2a(i,i));
%             end
%             cum3cell{k}=bk;
%         end
    case 2
        % third order cumulant of random factors
        C3t=@(ind) cum3(cum3id, x(ind(:,1),ind(:,2)), x(ind(:,3),ind(:,4)), x(ind(:,5),ind(:,6)));
        %out=ctt_greedycross(d, fun, option, ndoetemp,tol2)
        tic;
        cttout1 = ctt_greedycross(6, C3t, 1,800,1e-5);
        interpsets3=cttout1{1};
        cttout2=ctt_output(6,C3t,interpsets3);
        cum3t1=toc;
        C3tcell=cttout2(:,1);
        clear cttout2
        ndoe=10000;
        doe=rand(ndoe,6);
        refv=C3t(doe);
        cdata={fevalcbm(C3tcell{1},doe(:,1)),fevalcbm(C3tcell{2},doe(:,2)),fevalcbm(C3tcell{3},doe(:,3)),...
            fevalcbm(C3tcell{4},doe(:,4)),fevalcbm(C3tcell{5},doe(:,5)),fevalcbm(C3tcell{6},doe(:,6))};
        diffv=arrayfun(@(k) cdata{1}(:,:,k)*cdata{2}(:,:,k)*cdata{3}(:,:,k)*cdata{4}(:,:,k)*cdata{5}(:,:,k)*cdata{6}(:,:,k)-refv(k),(1:ndoe)');
        clear cdata
        epg3m2=sqrt(mean(diffv.^2)/mean(refv.^2));
        ranks3=[size(C3tcell{2},1),size(C3tcell{3},1),size(C3tcell{4},1),size(C3tcell{5},1),size(C3tcell{6},1)];
        tim2=toc;
        cum3t2=tim2-cum3t1;
        
        ni=size(mode2a,2);
        cum3cell=cell(1,3);
        for k=1:3
            sz1=size(C3tcell{2*k-1},1);
            sz2=size(C3tcell{2*k},2);
            ct=C3tcell{2*k-1};%chebmatrix
            mt=zeros(ni,sz2,sz1);
            for k2=1:sz1
                mt1=mul(mode1a',ct(k2,:),'par')*C3tcell{2*k};
                mt(:,:,k2) = mul(mode2a',mt1,'par');
                fprintf('%d  %d', k,k2);
                toc
            end
            cum3cell{k}=permute(mt,[3 1 2]);
        end
        tim3=toc;
        cum3t3=tim3-cum3t2;
    case 3
        % third order cumulant of random factors
        C3t=@(ind) cum3(cum3id, x(ind(:,1),ind(:,2),ind(:,3)), x(ind(:,4),ind(:,5),ind(:,6)), x(ind(:,7),ind(:,8),ind(:,9)));
        tic;
        cttout1=ctt_greedycross(9, C3t, 1,800,1e-4);
        interpsets3=cttout1{1};
%         save ex3_2.mat interpsets3
        cttout2=ctt_output(9,C3t,interpsets3);
        cum3t1=toc;
        C3tcell=cttout2(:,1);
        clear cttout2
        % clear cttout1
        ndoe=1000;
        doe=rand(ndoe,9);
        refv=C3t(doe);
        cdata={fevalcbm(C3tcell{1},doe(:,1)),fevalcbm(C3tcell{2},doe(:,2)),fevalcbm(C3tcell{3},doe(:,3)),...
            fevalcbm(C3tcell{4},doe(:,4)),fevalcbm(C3tcell{5},doe(:,5)),fevalcbm(C3tcell{6},doe(:,6)),...
            fevalcbm(C3tcell{7},doe(:,7)),fevalcbm(C3tcell{8},doe(:,8)),fevalcbm(C3tcell{9},doe(:,9))};
        diffv=arrayfun(@(k) cdata{1}(:,:,k)*cdata{2}(:,:,k)*cdata{3}(:,:,k)*cdata{4}(:,:,k)*cdata{5}(:,:,k)...
            *cdata{6}(:,:,k)*cdata{7}(:,:,k)*cdata{8}(:,:,k)*cdata{9}(:,:,k)-refv(k),(1:ndoe)');
        clear cdata
        epg3m3=sqrt(mean(diffv.^2)/mean(refv.^2));
        ranks3=[size(C3tcell{2},1),size(C3tcell{3},1),size(C3tcell{4},1),size(C3tcell{5},1),size(C3tcell{6},1),...
            size(C3tcell{7},1),size(C3tcell{8},1),size(C3tcell{9},1)];
        tim2=toc;
        cum3t2=tim2-cum3t1;
        
        ni=size(mode3a,2);
        cum3cell=cell(1,3);
        for k=1:3
            sz1=size(C3tcell{3*k-2},1);
            sz2=size(C3tcell{3*k},2);
            ct=C3tcell{3*k-2};%chebmatrix
            mt=zeros(ni,sz2,sz1);
            for k2=1:sz1
                mt1=mul(mode1a',ct(k2,:),'par')*C3tcell{3*k-1};
                mt1 = mul(mode2a',mt1,'par')*C3tcell{3*k};
                mt(:,:,k2) = mul(mode3a',mt1,'par');
                fprintf('%d  %d', k,k2);
                toc
            end
            cum3cell{k}=permute(mt,[3 1 2]);
        end
        tim3=toc;
        cum3t3=tim3-cum3t2;
end

szsC3=zeros(3,3);
for k=1:3
    if k==3
        szsC3(k,[1 2])=size(cum3cell{k});
        szsC3(k,3)=1;
    else
        szsC3(k,:)=size(cum3cell{k});
    end
end
t=cum3cell{3}*cum3cell{3}';
% L3=chol(t,'lower');
A=cum3cell{1}(1,:,:);
t2=0;
for k=1:szsC3(2,2)
    p=permute(cum3cell{2}(:,k,:),[1 3 2]);
    t2=t2+p*t*p';
end
t2m=0.5*(t2+t2');
p=permute(A,[2 3 1]);
t3=p*t2m*p';
t3m=0.5*(t3+t3');
[U3,S3]=eig(t3m);
[ss3,ind]=sort(abs(diag(S3)),'descend');
S3=diag(ss3);
U3=U3(:,ind);

nU3=size(U3,2);
tolcum=1-0.0001;
norm3=norm(diag(S3),2);
dS3=diag(S3);
for k=1:nU3
    r3k=norm(dS3(1:k),2)/norm3;
    U3f=U3(:,1:k);
    switch dimpar
        case 1
           norm2=norm(S2a,'fro');
           r2k=norm(U3f'*S2a*U3f,'fro')/norm2;
        case 2
           norm2=norm(S2b,'fro');
           r2k=norm(U3f'*S2b*U3f,'fro')/norm2; 
        case 3
           norm2=norm(S2c,'fro');
           r2k=norm(U3f'*S2c*U3f,'fro')/norm2; 
    end
    if min(r2k,r3k)>=tolcum
        break
    end
end
U3f=U3(:,1:k);
cum3core=cell(1,3);%third-order cumulants of third-order factors
for i=1:3
    cum3core{i}=cell2mat(cellfun(@(c) c*U3f, mat2cell(cum3cell{i},szsC3(i,1),szsC3(i,2),ones(1,szsC3(i,3))), 'UniformOutput', false));
end

%% Final assessment using global random test
% alpha(\xi,\eta,\zeta;\theta) =
% mode1a(\xi)*mode2a(\eta)*mode3a(\zeta)*U3f*\gamma(\theta)
switch dimpar
    case 1
        covgamma=U3f'*S2a*U3f;
        mode1b=mode1a*U3f;
        ndoe=10000;
        doe=rand(ndoe,2*dimpar);
        refv=C2t(doe);
        cdata={fevalcbm(mode1b,doe(:,1)),fevalcbm(mode1b,doe(:,2))};
        diffv=arrayfun(@(k) cdata{1}(:,:,k)*covgamma*cdata{2}(:,:,k)'-refv(k),(1:ndoe)');
        clear cdata
        epg2m3f=sqrt(mean(diffv.^2)/mean(refv.^2));
        
        doe=rand(ndoe,3*dimpar);
        refv=C3t(doe);
        nk=size(U3f,2);
        cdata2={fevalcbm(mode1b,doe(:,1)),fevalcbm(mode1b,doe(:,2)), fevalcbm(mode1b,doe(:,3))};
        
    case 2
        covgamma=U3f'*S2b*U3f;
        mode2b=mode2a*U3f;
        ndoe=10000;
        doe=rand(ndoe,2*dimpar);
        refv=C2t(doe);
        cdata={fevalcbm(mode1a,doe(:,1)),fevalcbm(mode2b,doe(:,2)),fevalcbm(mode2b,doe(:,3)),fevalcbm(mode1a,doe(:,4))};
        diffv=arrayfun(@(k) cdata{1}(:,:,k)*cdata{2}(:,:,k)*covgamma*cdata{3}(:,:,k)'*cdata{4}(:,:,k)'-refv(k),(1:ndoe)');
        clear cdata
        epg2m3f=sqrt(mean(diffv.^2)/mean(refv.^2));
        doe=rand(ndoe,3*dimpar);
        refv=C3t(doe);
        cdata={fevalcbm(mode1a,doe(:,1)),fevalcbm(mode2b,doe(:,2)), ...
            fevalcbm(mode1a,doe(:,3)),fevalcbm(mode2b,doe(:,4)),...
            fevalcbm(mode1a,doe(:,5)),fevalcbm(mode2b,doe(:,6))};
        nk=size(U3f,2);
        cdata2=cell(1,3);
        for i=1:3
            cdata2{i}=cell2mat(arrayfun(@(k) cdata{2*i-1}(:,:,k)*cdata{2*i}(:,:,k),(1:ndoe)','UniformOutput',false));
        end
        clear cdata 
    case 3
        covgamma=U3f'*S2c*U3f;
        mode3b=mode3a*U3f;
        ndoe=10000;
        doe=rand(ndoe,2*dimpar);
        refv=C2t(doe);
        cdata={fevalcbm(mode1a,doe(:,1)),fevalcbm(mode2a,doe(:,2)),fevalcbm(mode3b,doe(:,3)), fevalcbm(mode3b,doe(:,4)),fevalcbm(mode2a,doe(:,5)),fevalcbm(mode1a,doe(:,6))};
        diffv=arrayfun(@(k) cdata{1}(:,:,k)*cdata{2}(:,:,k)*cdata{3}(:,:,k)*covgamma*cdata{4}(:,:,k)'*cdata{5}(:,:,k)'*cdata{6}(:,:,k)'-refv(k),(1:ndoe)');
        clear cdata
        epg2m3f=sqrt(mean(diffv.^2)/mean(refv.^2));
        doe=rand(ndoe,3*dimpar);
        refv=C3t(doe);
        cdata={fevalcbm(mode1a,doe(:,1)),fevalcbm(mode2a,doe(:,2)),fevalcbm(mode3b,doe(:,3)), ...
            fevalcbm(mode1a,doe(:,4)),fevalcbm(mode2a,doe(:,5)),fevalcbm(mode3b,doe(:,6)),...
            fevalcbm(mode1a,doe(:,7)),fevalcbm(mode2a,doe(:,8)),fevalcbm(mode3b,doe(:,9))};
        nk=size(U3f,2);
        cdata2=cell(1,3);
        for i=1:3
            cdata2{i}=cell2mat(arrayfun(@(k) cdata{3*i-2}(:,:,k)*cdata{3*i-1}(:,:,k)*cdata{3*i}(:,:,k),(1:ndoe)','UniformOutput',false));
        end
        clear cdata        
end
cdata2=arrayfun(@(k) permute(cdata2{k},[3 2 1]),1:3,'UniformOutput',false);
cum3test=cell(1,3);
for i=1:3
    cum3test{i}=cell2mat(cellfun(@(c) c*cdata2{i}', mat2cell(cum3core{i},szsC3(i,1),nk,ones(1,szsC3(i,3))), 'UniformOutput', false));
    cum3test{i}=permute(cum3test{i},[1 3 2]);
end
diffv=arrayfun(@(k) cum3test{1}(:,:,k)*cum3test{2}(:,:,k)*cum3test{3}(:,:,k)-refv(k),(1:ndoe)');
epg3m3f=sqrt(mean(diffv.^2)/mean(refv.^2));


% fourth order cumulants of second-order factors
% switch dimpar
%     case 1
%         % fourth order cumulant of random factors
%         C4t=@(ind) cum4_2(x(ind(1)), x(ind(2)), x(ind(3)), x(ind(4)));
%         cttout=ctt_greedycross(4, C4t, 1);
%         C4tcell=cttout{1};
%         C4tcell{1,1}=chebmatrix(C4tcell{1,1});
%         interpsets4=cttout{2};
%         ctt4err=cttout{3};
%         clear cttout
%
%         mode1a=chebmatrix(mode1a);
%         ni=size(mode1a,2);
%         cum4cell=cell(1,4);
%         for k=1:4
%             sz1=size(C4tcell{k,1},1);
%             sz2=size(C4tcell{k,1},2);
%             id=num2cell(1:sz1);
%             cum4cell{k}=reshape(cell2mat(cellfun(@(c) mode1a'*(C4tcell{k,1}(c,:)*C4tcell{k,2}),id,'UniformOutput',0)),[ni,sz2,sz1]);
%             cum4cell{k}=permute(cum4cell{k},[3 1 2]);% in accordance with literature
%         end
%     case 2
%         % third order cumulant of random factors
%         C4t=@(ind) cum4_2(x(ind(1),ind(2)), x(ind(3),ind(4)), x(ind(5),ind(6)), x(ind(7),ind(8)));
%         cttout=ctt_greedycross(8, C4t, 1);
%         C4tcell=cttout{1};
%         C4tcell{1,1}=chebmatrix(C4tcell{1,1});
%         interpsets4=cttout{2};
%         ctt4err=cttout{3};
%         clear cttout
%
%         mode1a=chebmatrix(mode1a);
%         ni=size(mode2a,2);
%         cum4cell=cell(1,4);
%         for k=1:4
%             sz1=size(C4tcell{2*k-1,1},1);
%             sz2=size(C4tcell{2*k,2},2);
%             id=num2cell(1:sz1);
%             cum4cell{k}=reshape(cell2mat(cellfun(@(c) mode2a'*(mode1a'*(C4tcell{2*k-1,1}(c,:)*C4tcell{2*k-1,2}))*(C4tcell{2*k,1}*C4tcell{2*k,2}),...
%                 id,'UniformOutput',0)),[ni,sz2,sz1]);
%             cum4cell{k}=permute(cum4cell{k},[3 1 2]);
%         end
%     case 3
%         % third order cumulant of random factors
%         C4t=@(ind) cum3_2(x(ind(1),ind(2),ind(3)), x(ind(4),ind(5),ind(6)), x(ind(7),ind(8),ind(9)), x(ind(10),ind(11),ind(12)));
%         cttout=ctt_greedycross(12, C4t, 1);
%         C4tcell=cttout{1};
%         C4tcell{1,1}=chebmatrix(C4tcell{1,1});
%         interpsets4=cttout{2};
%         ctt4err=cttout{3};
%         clear cttout
%
%         mode1a=chebmatrix(mode1a);
%         ni=size(mode3a,2);
%         cum4cell=cell(1,4);
%         for k=1:4
%             sz1=size(C4tcell{3*k-2,1},1);
%             sz2=size(C4tcell{3*k,1},2);
%             id=num2cell(1:sz1);
%             cum4cell{k}=reshape(cell2mat(cellfun(@(c) (mode3a'*(mode2a'*(mode1a'*(C4tcell{3*k-2,1}(c,:)*C4tcell{3*k-2,2}))*(C4tcell{3*k-1,1}*C4tcell{3*k-1,2}))*(C4tcell{3*k,1}*C4tcell{3*k,2})),...
%                 id,'UniformOutput',0)),[ni,sz2,sz1]);
%             cum4cell{k}=permute(cum4cell{k},[3 1 2]);
%         end
% end
%
% szsC4=zeros(4,3);
% for k=1:3
%     szsC4(k,:)=size(cum4cell{k});
% end
% cum4core=cell(1,4);%fourth order cumulants of third-order factors
% for i=1:4
%     cum4core{i}=cellfun(@(c) c*U3, mat2cell(cum4cell{i},szsC4(i,1),szsC4(i,2),ones(1,szsC4(i,3))));
% end