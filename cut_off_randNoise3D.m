function [x_new, y_new, z_new]=cut_off_randNoise3D(x, y, z, N, deg)
bound=floor(0.40*N);
r=sqrt((x-(N+1)*0.5)^2+(y-(N+1)*0.5)^2+(z-(N+1)*0.5)^2);
if r>bound
     x_new=x;
     y_new=y;
     z_new=z;
else
    x_new=x+deg*(sin(rand(1)-rand(1)));
    y_new=y+deg*(sin(rand(1)-rand(1)));
    z_new=z+deg*(sin(rand(1)-rand(1)));
end