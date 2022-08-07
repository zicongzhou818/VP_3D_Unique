function [x_new, y_new, z_new]=shift_var(x,y,z,N,b,var)
cut_off_bound=floor(0.5*N);
r=sqrt((x-(N+1)*0.5)^2+(y-(N+1)*0.5)^2+(z-(N+1)*0.5)^2);
    if r>cut_off_bound
         x_new=x;
         y_new=y;
         z_new=z;
    else
        if var==1
            a=(cos((pi*r)/(cut_off_bound))+1)*0.5*b;
            x_temp=(x)/N; 
            x_temp=(x_temp+a*0.5)*N; 
            x_new=x_temp; 
            y_new=y;
            z_new=z;
        end
        if var==2
            a=(cos((pi*r)/(cut_off_bound))+1)*0.5*b;
            x_new=x; 
            y_temp=(y)/N; 
            y_temp=(y_temp+a*0.5)*N; 
            y_new=y_temp; 
            z_new=z;
        end
        if var==3
            a=(cos((pi*r)/(cut_off_bound))+1)*0.5*b;
            x_new=x; 
            y_new=y;
            z_temp=(z)/N; 
            z_temp=(z_temp+a*0.5)*N;  
            z_new=z_temp; 
        end
    end
end