function [f,s]=svdre3(ym1,ym2)
%
tol=1e-15;
[QL,RL]=qr_h3(ym1,tol);
ym2=cellfun(@transpose,ym2)';
% ym2=chebmatrix(cheb2cell(ym2)');
[~,RR]=qr_h3(ym2,tol);
% QR=cellfun(@transpose,QR)';
[U,S,~]=svd(RL*RR');
[ss,ind]=sort(abs(diag(S)),'descend');
% out={QL*U(:,ind),diag(ss)};
% out={QL*U(:,ind),diag(ss),QR*V(ind,:)};
f=QL*U(:,ind);
s=diag(ss);
end