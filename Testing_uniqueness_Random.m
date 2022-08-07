close all
clear
clc


N=51;k1=37;k2=25;
Npts=N-2; 
% number of points each direction
% h1=1/(N-1);
h=1;
[x1,x2,x3]=ndgrid(1:N,1:N,1:N);
phi1=x1;
phi2=x2;
phi3=x3;

Phi1=x1;
Phi2=x2;
Phi3=x3;

phi1_temp=x1;
phi2_temp=x2;
phi3_temp=x3;
for i=1:N
    for j=1:N
        for k=1:N
            x=phi1(i,j,k);
            y=phi2(i,j,k);
            z=phi3(i,j,k);
            [x_new, y_new, z_new]=cut_off_small_3D(x,y,z,N,-pi/9,0.65,0.35,0.55);
            phi1(i,j,k)=x_new;
            phi2(i,j,k)=y_new;
            phi3(i,j,k)=z_new;
        end
    end
end
for i=1:N
    for j=1:N
        for k=1:N
            x=phi1(i,j,k);
            y=phi2(i,j,k);
            z=phi3(i,j,k);
            [x_new, y_new, z_new]=cut_off_small_3D(x,y,z,N,-pi/9,0.35,0.65,0.25);
            phi1(i,j,k)=x_new;
            phi2(i,j,k)=y_new;
            phi3(i,j,k)=z_new;
        end
    end
end
for i=1:N
    for j=1:N
        for k=1:N
            x=phi1(i,j,k);
            y=phi2(i,j,k);
            z=phi3(i,j,k);
            [x_new, y_new, z_new]=cut_off_3D(x,y,z,N,pi/4);
            [x_new, y_new, z_new]=shift_var(x_new, y_new, z_new,N,-0.15,1);
            [x_new, y_new, z_new]=shift_var(x_new, y_new, z_new,N,0.15,2);
            [x_new, y_new, z_new]=shift_var(x_new, y_new, z_new,N,0.05,3);
            phi1_temp(i,j,k)=x_new;
            phi2_temp(i,j,k)=y_new;
            phi3_temp(i,j,k)=z_new;
        end
    end
end



Ntest=12; %number of tests with different gaussian noises

mag_Up_max=zeros(1,Ntest);
mag_Up_min=zeros(1,Ntest);
FormJDmin=zeros(1,Ntest);
max_JDT=zeros(1,Ntest);
max_CVT=zeros(1,Ntest);
max_mag_distT=zeros(1,Ntest);
diff_Pn_U_max=zeros(1,Ntest);
diff_cv_U_max=zeros(1,Ntest);
dist_U_max=zeros(1,Ntest);
sumTime=zeros(1,Ntest);

for pick=1:Ntest
phi1=phi1_temp;
phi2=phi2_temp;
phi3=phi3_temp;
    for i=1:N
        for j=1:N
            for k=1:N
                x=phi1(i,j,k);
                y=phi2(i,j,k);
                z=phi3(i,j,k);
                [x_new, y_new, z_new]=cut_off_randNoise3D(x, y, z, N, 0.35);
                phi1(i,j,k)=x_new;
                phi2(i,j,k)=y_new;
                phi3(i,j,k)=z_new;
            end
        end
    end
    
    pPhi1=phi1;
    pPhi2=phi2;
    pPhi3=phi3;
    Up1 = pPhi1-x1;
    Up2 = pPhi2-x2;
    Up3 = pPhi3-x3;
    mag_Up_max(1,pick)=max(max(max(sqrt(Up1.^2+Up2.^2+Up3.^2))));
    mag_Up_min(1,pick)=min(min(min(sqrt(Up1.^2+Up2.^2+Up3.^2))));
    
    % mag_Up_max =  10.220883703039281
    % mag_Up_min =  0
    
    [JD_Phi, CV_Phi1, CV_Phi2, CV_Phi3]=compute_JD_and_Curl3D(pPhi1,pPhi2,pPhi3,h);
    FormJDmin(1,pick)=min(min(min(JD_Phi)));
    h2=1;
    % figure(1)
    % gridplot3D_flexible(pPhi1,pPhi2,pPhi3,h2,h2,h2),xlabel('x'),ylabel('y'),zlabel('z');
%     figure(2)
%     gridplot3D_framed_axis(pPhi1,pPhi2,pPhi3,1,k1,1),
%     gridplot3D_framed_axis(pPhi1,pPhi2,pPhi3,1,k2,3), grid on
%     axis([0, N+1, 0, N+1, 0, N+1]);
%     xlabel('x'),ylabel('y'),zlabel('z');
%     figure(3)
%     quiver3(x1(:, :, k2), x2(:, :, k2), x3(:, :, k2), Up1(:, :, k2),Up2(:, :, k2),Up3(:, :, k2), 0,'g'), hold on
%     quiver3(x1(k1, :, :), x2(k1, :, :), x3(k1, :, :), Up1(k1, :, :),Up2(k1, :, :),Up3(k1, :, :), 0,'g'), grid on;
%     axis([0, N+1, 0, N+1, 0, N+1]);
%     xlabel('x'),ylabel('y'),zlabel('z');
    
    
    %% iterative method old
    tstep=1e-4;
    tic
    [phi1_id,phi2_id,phi3_id,U1,U2,U3]=Check_Id_from_given3D(N,pPhi1,pPhi2,pPhi3);
    toc
    sumTime(1,pick)=toc;
    
%     figure(4)
%     gridplot3D_framed_axisR(phi1_id,phi2_id,phi3_id,1,k1,1),
%     gridplot3D_framed_axisR(phi1_id,phi2_id,phi3_id,1,k2,3), grid on
%     axis([0, N+1, 0, N+1, 0, N+1]);
%     xlabel('x'),ylabel('y'),zlabel('z');
%     
%     figure(5)
%     quiver3(x1(:, :, k2), x2(:, :, k2), x3(:, :, k2), U1(:, :, k2),U2(:, :, k2),U3(:, :, k2), 0,'b'), hold on
%     quiver3(x1(k1, :, :), x2(k1, :, :), x3(k1, :, :), U1(k1, :, :),U2(k1, :, :),U3(k1, :, :), 0,'b'), grid on;
%     axis([0, N+1, 0, N+1, 0, N+1]);
%     xlabel('x'),ylabel('y'),zlabel('z');
    
    % view(0.0, 0.0);
    
    [JD_T, CV_T1, CV_T2, CV_T3]=compute_JD_and_Curl3D(phi1_id,phi2_id,phi3_id,1);
    
    diff_JDT=abs(JD_T-ones(N,N,N));
    diff_CV1T=abs(CV_T1-zeros(N,N,N));
    diff_CV2T=abs(CV_T2-zeros(N,N,N));
    diff_CV3T=abs(CV_T3-zeros(N,N,N));
    max_JDT(1,pick)=max(max(max(diff_JDT)));
    mag_CVT=(diff_CV1T.^2+diff_CV2T.^2+diff_CV3T.^2).^(.5);
    max_CVT(1,pick)=max(max(max(mag_CVT)));
    
    diff_x1T=phi1_id-x1;
    diff_x2T=phi2_id-x2;
    diff_x3T=phi3_id-x3;
    mag_distT=sqrt(diff_x1T.^2+diff_x2T.^2+diff_x3T.^2);
    max_mag_distT(1,pick)=max(max(max(mag_distT)));
    
    [JD_U, CV1_U, CV2_U, CV3_U]=compute_JD_and_Curl3D(U1,U2,U3,1);
    [u1y,u1x,u1z] = gradient(U1,1);
    [u2y,u2x,u2z] = gradient(U2,1);
    [u3y,u3x,u3z] = gradient(U3,1);
    tail_u = u2y.*u3z - u3y.*u2z + u1x.*u3z - u3x.*u1z + u1x.*u2y - u2x.*u1y;
    div_u = u1x+u2y+u3z;
    diff_Pn=abs(div_u+JD_U+tail_u);
    diff_CV1_U=abs(CV1_U-zeros(N,N,N));
    diff_CV2_U=abs(CV2_U-zeros(N,N,N));
    diff_CV3_U=abs(CV3_U-zeros(N,N,N));
    mag_CV_U=sqrt(diff_CV1_U.^2+diff_CV2_U.^2+diff_CV3_U.^2);
    diff_Pn_U_max(1,pick)=max(max(max(diff_Pn)));
    diff_cv_U_max(1,pick)=max(max(max(mag_CV_U)));
    dist_U_max(1,pick)=max(max(max(sqrt(U1.^2+U2.^2+U3.^2))));



end


Sum_Mean_mag_Up_max=mean(mag_Up_max)
Sum_Mean_mag_Up_min=mean(mag_Up_min)
Sum_Mean_FormJDmin=mean(FormJDmin)
Sum_Mean_max_JDT=mean(max_JDT)
Sum_Mean_max_CVT=mean(max_CVT)
Sum_Mean_max_mag_distT=mean(max_mag_distT)
Sum_Mean_diff_Pn_U_max=mean(diff_Pn_U_max)
Sum_Mean_diff_cv_U_max=mean(diff_cv_U_max)
Sum_Mean_dist_U_max=mean(dist_U_max)
Sum_Mean_sumTime=mean(sumTime)

Sum_Var_mag_Up_max=var(mag_Up_max)
Sum_Var_mag_Up_min=var(mag_Up_min)
Sum_Var_FormJDmin=var(FormJDmin)
Sum_Var_max_JDT=var(max_JDT)
Sum_Var_max_CVT=var(max_CVT)
Sum_Var_max_mag_distT=var(max_mag_distT)
Sum_Var_diff_Pn_U_max=var(diff_Pn_U_max)
Sum_Var_diff_cv_U_max=var(diff_cv_U_max)
Sum_Var_dist_U_max=var(dist_U_max)
Sum_Var_sumTime=var(sumTime)
