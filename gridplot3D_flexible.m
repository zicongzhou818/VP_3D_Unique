function [] = gridplot3D_flexible(X,Y,Z,hx,hy,hz)
[n1,n2,n3]=size(X);
% X=imresize3(X,[n1+1, n2+1, n3+1]);
% Y=imresize3(Y,[n1+1, n2+1, n3+1]);
% Z=imresize3(Z,[n1+1, n2+1, n3+1]);
for i=1:hx:n1
    for j=1:hy:n2
        % z-axes lines
        x=squeeze(X(i,j,:));
        y=squeeze(Y(i,j,:));
        z=squeeze(Z(i,j,:));
        plot3(x,y,z,'k'),hold on
    end
end
for i=1:hx:n1
    for k=1:hz:n3       
        x=squeeze(X(i,:,k));
        y=squeeze(Y(i,:,k));
        z=squeeze(Z(i,:,k));
        plot3(x,y,z,'k'),hold on
    end
end
for j=1:hy:n2
    for k=1:hz:n3
        x=squeeze(X(:,j,k));
        y=squeeze(Y(:,j,k));
        z=squeeze(Z(:,j,k));
        plot3(x,y,z,'k'),hold on
    end
end
