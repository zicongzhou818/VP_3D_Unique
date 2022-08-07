function [x_new, y_new, z_new]=cut_off_small_3D(x,y,z,N,theta,X_rat,Y_rat,Z_rat)
cut_off_bound=floor(0.25*N);
r=sqrt((x-(N+1)*X_rat)^2+(y-(N+1)*Y_rat)^2+(z-(N+1)*Z_rat)^2);
    if r>cut_off_bound
         x_new=x;
         y_new=y;
         z_new=z;
    else
        a=(cos((pi*r)/(cut_off_bound))+1)*0.5;
        u=[x-(N+1)*0.5; y-(N+1)*0.5; z-(N+1)*0.5];
        Alpha=a*theta;
        u=[cos(Alpha) -sin(Alpha) 0;...
           sin(Alpha) cos(Alpha)  0;...
           0           0          1]*u;
       
        u=[cos(Alpha) 0 -sin(Alpha);...
           0          1          0;...
           sin(Alpha) 0 cos(Alpha)]*u;
       
        u=[1          0          0; ...
           0      cos(Alpha) -sin(Alpha);...
           0      sin(Alpha) cos(Alpha)]*u;
        x_new=u(1)+(N+1)*0.5;
        y_new=u(2)+(N+1)*0.5;
        z_new=u(3)+(N+1)*0.5;
    end
end