function nrb=igf2d(id)
switch id
    case 1%2D parametric space example: an annulus
        %inner and outer radius
        ro=2;
        ri=1;
        %knots
        knots=cell(1,2);
        knots{1}=[0 0 1 1];
        knots{2}=[0 0 0 1 1 2 2 3 3 4 4 4]./4;
        %unweighted control points
        coefs=zeros(4,2,9);
        coefs(1,1,:)=[ri ri 0 -ri -ri -ri 0 ri ri];
        coefs(1,2,:)=[ro ro 0 -ro -ro -ro 0 ro ro];
        coefs(2,1,:)=[0 ri ri ri 0 -ri -ri -ri 0];
        coefs(2,2,:)=[0 ro ro ro 0 -ro -ro -ro 0];
        %coefs(3,:,:) all zeros
        %weights
        coefs(4,1,:)=[1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1];
        coefs(4,2,:)=[1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1];
        %weighted control points
        coefs(1,:,:)=coefs(1,:,:).*coefs(4,:,:);
        coefs(2,:,:)=coefs(2,:,:).*coefs(4,:,:);
        %make an NURBS geometry object
        nrb = nrbmak(coefs, knots);
        
    case 2
        L=1;
        H=1;        
        %knots
        knots=cell(1,2);
        knots{1}=[0 0 1 1];
        knots{2}=[0 0 1 1];        
        %unweighted control points
        coefs=zeros(4,2,2);
        coefs(1,:,:)=[-L/2,L/2;-L/2,L/2];       
        coefs(2,:,:)=[-L/2,-L/2;L/2,L/2];
        coefs(3,:,:)=[0,H;H,0];
        %weights
        coefs(4,:,:)=[1 1;1 1];
        %weighted control points
        coefs(1,:,:)=coefs(1,:,:).*coefs(4,:,:);
        coefs(2,:,:)=coefs(2,:,:).*coefs(4,:,:);
        coefs(3,:,:)=coefs(3,:,:).*coefs(4,:,:);
        
        %make an NURBS geometry object
        nrb = nrbmak(coefs, knots);
end
end