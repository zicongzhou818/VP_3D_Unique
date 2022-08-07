function  [phi1,phi2,phi3,U1,U2,U3]=Check_Id_from_given3D(N,given_x,given_y,given_z)
%% initialization
better = 1;
ei = 0;
r = 1;
ti = 0;
ts = 1e-3;
Npts = N-2;
imax = 5.0e2;
ts_r = 1e-6;
rTol = 1.0e-8;
%% preset control functions corresponsive to the identity map
f1 = zeros(Npts,Npts,Npts);
f2 = zeros(Npts,Npts,Npts);
f3 = zeros(Npts,Npts,Npts);
%% preset new deformation based on control functions that are corresponsive to the identity map
[X,Y,Z] = ndgrid(1:N,1:N,1:N);
u_1 = given_x-X;
u_2 = given_y-Y;
u_3 = given_z-Z;
%% reading prescriptions
[u1y,u1x,u1z] = gradient(u_1,1);
[u2y,u2x,u2z] = gradient(u_2,1);
[u3y,u3x,u3z] = gradient(u_3,1);
tail_u = u2y.*u3z - u3y.*u2z + u1x.*u3z - u3x.*u1z + u1x.*u2y - u2x.*u1y;
div_u = u1x + u2y + u3z;
[jd_u, cv1_u, cv2_u, cv3_u]=compute_JD_and_Curl3D(u_1 , u_2, u_3, 1);
%% (ps: the final phi initially is the same as given phi here because the current phi is id)
%% compute ssd
Pn=(div_u + jd_u + tail_u);
ssd_initial=sum(sum(sum(Pn(2:Npts+1,2:Npts+1,2:Npts+1).^2+cv1_u(2:Npts+1,2:Npts+1,2:Npts+1).^2+cv2_u(2:Npts+1,2:Npts+1,2:Npts+1).^2+cv3_u(2:Npts+1,2:Npts+1,2:Npts+1).^2)));
ssd_old=ssd_initial;

%% main loop
while ti<imax && ts>ts_r && r>rTol
    ti = ti+1;
    if better % form gradients wrt F when ssd decreases
        ei = ei+1;
        [Pny,Pnx,Pnz] = gradient(Pn,1);

        [~, a1x, ~] = gradient(Pn.*(u2y.*u3z - u3y.*u2z),1);
        [a1y, ~, ~] = gradient(Pn.*(u3x.*u2z - u2x.*u3z),1);
        [~, ~, a1z] = gradient(Pn.*(u2x.*u3y - u2y.*u3x),1);

        [~, a2x, ~] = gradient(Pn.*(u3y.*u1z - u1y.*u3z),1);
        [a2y, ~, ~] = gradient(Pn.*(u1x.*u3z - u1y.*u3x),1);
        [~, ~, a2z] = gradient(Pn.*(u3x.*u1y - u1x.*u3y),1);

        [~, a3x, ~] = gradient(Pn.*(u1y.*u2z - u2y.*u1z),1);
        [a3y, ~, ~] = gradient(Pn.*(u2x.*u1z - u1x.*u2z),1);
        [~, ~, a3z] = gradient(Pn.*(u1x.*u2y - u2x.*u1y),1);

        del_div_u1 = a1x+a1y+a1z;
        del_div_u2 = a2x+a2y+a2z;
        del_div_u3 = a3x+a3y+a3z;
        
        [~, b1x, ~] = gradient(u3z + u2z,1);
        [b1y, ~, ~] = gradient(-u2x,1);
        [~, ~, b1z] = gradient( -u3x,1);

        [~, b2x, ~] = gradient(- u1y,1);
        [b2y, ~, ~] = gradient(u3z + u1x,1);
        [~, ~, b2z] = gradient(- u3y,1);

        [~, b3x, ~] = gradient( -u1z,1);
        [b3y, ~, ~] = gradient(Pn.*(-u2z),1);
        [~, ~, b3z] = gradient(u2y + u1x,1);

        del_tail_u1 = b1x+b1y+b1z;
        del_tail_u2 = b2x+b2y+b2z;
        del_tail_u3 = b3x+b3y+b3z;

        lap_a1 = Pnx + del_div_u1 + del_tail_u1;
        lap_a2 = Pny + del_div_u2 + del_tail_u2;
        lap_a3 = Pnz + del_div_u3 + del_tail_u3;

        a1=pois3fft(-lap_a1);
        a2=pois3fft(-lap_a2);
        a3=pois3fft(-lap_a3);
        
        %% form c
        [cv1_uy,~,cv1_uz] = gradient(cv1_u,1);
        [~,cv2_ux,cv2_uz] = gradient(cv2_u,1);
        [cv3_uy,cv3_ux,~] = gradient(cv3_u,1);
        lap_c1 = - cv3_uy + cv2_uz;
        lap_c2 =   cv3_ux - cv1_uz;
        lap_c3 = - cv2_ux + cv1_uy;

        c1=pois3fft(-lap_c1);
        c2=pois3fft(-lap_c2);
        c3=pois3fft(-lap_c3);

        del_f1 = a1 + c1;
        del_f2 = a2 + c2;
        del_f3 = a3 + c3;
    end
    %update f1n,f2n to calculate u1,u2
    f1n=f1-ts*del_f1(2:Npts+1,2:Npts+1,2:Npts+1);
    f2n=f2-ts*del_f2(2:Npts+1,2:Npts+1,2:Npts+1);
    f3n=f3-ts*del_f3(2:Npts+1,2:Npts+1,2:Npts+1);
    %update u1,u2
    u1=pois3fft(f1n);
    u2=pois3fft(f2n);
    u3=pois3fft(f3n);
    U1 = matrixpad3D(u1,0);
    U2 = matrixpad3D(u2,0);
    U3 = matrixpad3D(u3,0);
    %update JT,curlT
    [jd_u, cv1_u, cv2_u, cv3_u] = compute_JD_and_Curl3D(U1 , U2, U3, 1);
    [u1y,u1x,u1z] = gradient(U1,1);
    [u2y,u2x,u2z] = gradient(U2,1);
    [u3y,u3x,u3z] = gradient(U3,1);
    tail_u = u2y.*u3z - u3y.*u2z + u1x.*u3z - u3x.*u1z + u1x.*u2y - u2x.*u1y;
    div_u = u1x + u2y + u3z;
    Pn=(div_u + jd_u + tail_u);
    ssd_new=sum(sum(sum(Pn(2:Npts+1,2:Npts+1,2:Npts+1).^2+cv1_u(2:Npts+1,2:Npts+1,2:Npts+1).^2+cv2_u(2:Npts+1,2:Npts+1,2:Npts+1).^2+cv3_u(2:Npts+1,2:Npts+1,2:Npts+1).^2)));
    r=ssd_new/ssd_initial;
%     display([' ts: ',num2str(ts),' r: ',num2str(r), ' ei: ',num2str(ei), ' ti: ',num2str(ti)]);
    if (ssd_new<ssd_old)
        ts=ts*1.005;
        f1=f1n;
        f2=f2n;
        f3=f3n;
        ssd_old=ssd_new;
        better=1;
    else
        better=0;
        ts=ts*0.995;
    end 
end
display([' ts: ',num2str(ts),' r: ',num2str(r), ' ei: ',num2str(ei), ' ti: ',num2str(ti)]);
phi1=X+U1;
phi2=Y+U2;
phi3=Z+U3;
end
