function nrb=igf3d(id)
switch id 
    case 1
       %inner and outer radius
        ro=2;
        ri=1;
        R=4;
        %knots
        knots=cell(1,3);
        knots{1}=[0 0 1 1];
        knots{2}=[0 0 0 1 1 2 2 3 3 4 4 4]./4;       
        knots{3}=[0 0 0 1 1 1];
        %unweighted control points
        coefs=zeros(4,2,9,3);
        coefs(1,1,:,1)=[-R+ri -R+ri -R+0 -R+(-ri) -R+(-ri) -R+(-ri) -R+0 -R+ri -R+ri];
        coefs(1,2,:,1)=[-R+ro -R+ro -R+0 -R+(-ro) -R+(-ro) -R+(-ro) -R+0 -R+ro -R+ro];
        coefs(1,1,:,2)=[-R+ri -R+ri -R+0 -R+(-ri) -R+(-ri) -R+(-ri) -R+0 -R+ri -R+ri];
        coefs(1,2,:,2)=[-R+ro -R+ro -R+0 -R+(-ro) -R+(-ro) -R+(-ro) -R+0 -R+ro -R+ro];
        coefs(1,:,:,3)=zeros(2,9);
        coefs(2,1,:,1)=[0 ri ri ri 0 -ri -ri -ri 0];
        coefs(2,2,:,1)=[0 ro ro ro 0 -ro -ro -ro 0];
        coefs(2,1,:,2)=[0 ri ri ri 0 -ri -ri -ri 0];
        coefs(2,2,:,2)=[0 ro ro ro 0 -ro -ro -ro 0];
        coefs(2,1,:,3)=[0 ri ri ri 0 -ri -ri -ri 0];
        coefs(2,2,:,3)=[0 ro ro ro 0 -ro -ro -ro 0];
        coefs(3,:,:,1)=zeros(2,9);
        coefs(3,1,:,2)=[-R+ri -R+ri -R+0 -R+(-ri) -R+(-ri) -R+(-ri) -R+0 -R+ri -R+ri];
        coefs(3,2,:,2)=[-R+ro -R+ro -R+0 -R+(-ro) -R+(-ro) -R+(-ro) -R+0 -R+ro -R+ro];
        coefs(3,1,:,3)=[-R+ri -R+ri -R+0 -R+(-ri) -R+(-ri) -R+(-ri) -R+0 -R+ri -R+ri];
        coefs(3,2,:,3)=[-R+ro -R+ro -R+0 -R+(-ro) -R+(-ro) -R+(-ro) -R+0 -R+ro -R+ro];
        %weights
        coefs(4,1,:,1)=[1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1];
        coefs(4,2,:,1)=[1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1];
        coefs(4,1,:,2)=[1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1]*sqrt(0.5);
        coefs(4,2,:,2)=[1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1]*sqrt(0.5);
        coefs(4,1,:,3)=[1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1];
        coefs(4,2,:,3)=[1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1 sqrt(0.5) 1];

    case 2
        %inner and outer radius
        ro=1.2;
        ri=1;
        %knots
        knots=cell(1,3);
        knots{1}=[0 0 1 1];
        knots{2}=[0 0 0 1 1 1];       
        knots{3}=[0 0 0 1 1 1];
        %unweighted control points
        coefs=zeros(4,2,3,3);
        coefs(1,:,:,1)=[ri ro;ri ro;0 0]';
        coefs(1,:,:,2)=[ri ro;ri ro;0 0]';
        coefs(1,:,:,3)=zeros(2,3);
        coefs(2,:,:,1)=[0 0;ri ro;ri ro]';
        coefs(2,:,:,2)=[0 0;ri ro;ri ro]';
        coefs(2,:,:,3)=zeros(2,3);
        coefs(3,:,:,1)=zeros(2,3);
        coefs(3,:,:,2)=[ri ro;ri ro;ri ro]';
        coefs(3,:,:,3)=[ri ro;ri ro;ri ro]';
        %weights
        coefs(4,:,:,1)=[1 1; sqrt(0.5) sqrt(0.5); 1 1]';
        coefs(4,:,:,2)=[1 1; sqrt(0.5) sqrt(0.5); 1 1]'*sqrt(0.5);
        coefs(4,:,:,3)=[1 1; sqrt(0.5) sqrt(0.5); 1 1]';
end
        %weighted control points
        coefs(1,:,:,:)=coefs(1,:,:,:).*coefs(4,:,:,:);
        coefs(2,:,:,:)=coefs(2,:,:,:).*coefs(4,:,:,:);
        coefs(3,:,:,:)=coefs(3,:,:,:).*coefs(4,:,:,:);
        %make an NURBS geometry object
        nrb = nrbmak(coefs, knots); 
end