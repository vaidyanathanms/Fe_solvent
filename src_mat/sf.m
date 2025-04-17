clear all
close all


tot_atom=23370;
Mg_num=270;
TFSI_num=540;
water_num=5000;
box=6.46285;
R=2.5;
delta_q=0.1;

neutral_factor=dlmread('neutral.dat','',1,0);
ion_factor=dlmread('ions.dat','',1,0);
for i=1:7
    for j=1:7
        P{i,j}=dlmread(['./rdf_results/rdf_',num2str(i,'%d'),num2str(j,'%d'),'.xvg'],'',24,0);
    end
end

b1=ion_factor(:,[1,9]); %Li
b2=neutral_factor(:,[1,10]);%F
b3=neutral_factor(:,[1,17]);%S
b4=neutral_factor(:,[1,9]);%O
b5=neutral_factor(:,[1,8]);%N
b6=neutral_factor(:,[1,7]);%C
b7=neutral_factor(:,[1,2]);%H
x1=Mg_num/tot_atom;
x2=6*TFSI_num/tot_atom;
x3=2*TFSI_num/tot_atom;
x4=(4*TFSI_num+water_num)/tot_atom;
x5=TFSI_num/tot_atom;
x6=2*TFSI_num/tot_atom;
x7=2*water_num/tot_atom;
rho=tot_atom/box^3; %number/volume

b1_new=[[0:delta_q:20]',interp1(b1(1:56,1),b1(1:56,2),[0:delta_q/10:2])'];
b2_new=[[0:delta_q:20]',interp1(b2(1:56,1),b2(1:56,2),[0:delta_q/10:2])'];
b3_new=[[0:delta_q:20]',interp1(b3(1:56,1),b3(1:56,2),[0:delta_q/10:2])'];
b4_new=[[0:delta_q:20]',interp1(b4(1:56,1),b4(1:56,2),[0:delta_q/10:2])'];
b5_new=[[0:delta_q:20]',interp1(b5(1:56,1),b5(1:56,2),[0:delta_q/10:2])'];
b6_new=[[0:delta_q:20]',interp1(b6(1:56,1),b6(1:56,2),[0:delta_q/10:2])'];
b7_new=[[0:delta_q:20]',interp1(b7(1:56,1),b7(1:56,2),[0:delta_q/10:2])'];

denominator_temp1=x1.*b1_new(:,2)+x2.*b2_new(:,2)+x3.*b3_new(:,2)+x4.*b4_new(:,2)+x5.*b5_new(:,2)+x6.*b6_new(:,2)+x7.*b7_new(:,2);
denominator_temp2=denominator_temp1.*denominator_temp1;
denominator=[b1_new(:,1), denominator_temp2];

x_tot=[x1, x2, x3, x4, x5, x6, x7];
b_tot=[b1_new(:,2), b2_new(:,2), b3_new(:,2), b4_new(:,2), b5_new(:,2), b6_new(:,2), b7_new(:,2)];
S_tot=zeros(1,length(b1_new(:,1)));
for i=1:7
for j=1:7
    for q=0:delta_q:20
        index_q=round(q/delta_q)+1;
       	for r=0.001:0.002:2.501
           index_r=floor(r/0.002)+1;
           %temp(index_r)=4*pi*r^2*(P{i,j}(index_r,2)-1)*sin(q*r)/(q*r)*sin(pi*R)*R/(pi*r);
           temp(index_r)=4*pi*r^2*(P{i,j}(index_r,2)-1)*sin(q*r)/(q*r)*sin(pi*r/R)*R/(pi*r);
           %temp(index_r)=4*pi*r^2*(P{i,j}(index_r,2)-1)*sin(q*r)/(q*r);
       	end
        %S_test{i,j}(index_q)=rho*x_tot(i)*x_tot(j)*b_tot(index_q,i)*b_tot(index_q,j)*sum(temp.*[0:0.002:2.5])/denominator(index_q,2);
        S_test{i,j}(index_q)=rho*x_tot(i)*x_tot(j)*b_tot(index_q,i)*b_tot(index_q,j)*sum(temp.*0.002)/denominator(index_q,2);
    end
    S_tot=S_tot+S_test{i,j};
end
end

xaxis=b1_new(:,1)';
semilogx(xaxis/10,S_tot);
axis([0.1 2 -1.2 0.5]);

save x_01.dat xaxis -ascii
save y_01.dat S_tot -ascii

%axis([4 20 -7000 3000]);

figure(2) %Li
% for i=1:7
%     for j=1:7
%         plot(xaxis,S_test{i,j});hold on
%     end
% end
plot(xaxis,S_test{1,1},'-r');hold on
plot(xaxis,S_test{1,2},'-b');hold on
plot(xaxis,S_test{1,3},'-g');hold on
plot(xaxis,S_test{1,4},'-m');hold on
plot(xaxis,S_test{1,5},'-c');hold on
plot(xaxis,S_test{1,6},'-k');hold on
plot(xaxis,S_test{1,7},'-y');hold on

figure(3) %H
plot(xaxis,S_test{7,1},'-r');hold on
plot(xaxis,S_test{7,2},'-b');hold on
plot(xaxis,S_test{7,3},'-g');hold on
plot(xaxis,S_test{7,4},'-m');hold on
plot(xaxis,S_test{7,5},'-c');hold on
plot(xaxis,S_test{7,6},'-k');hold on
plot(xaxis,S_test{7,7},'-y');hold on

figure(4) %F
plot(xaxis,S_test{2,1},'-r');hold on
plot(xaxis,S_test{2,2},'-b');hold on
plot(xaxis,S_test{2,3},'-g');hold on
plot(xaxis,S_test{2,4},'-m');hold on
plot(xaxis,S_test{2,5},'-c');hold on
plot(xaxis,S_test{2,6},'-k');hold on
plot(xaxis,S_test{2,7},'-y');hold on

figure(5) %S
plot(xaxis,S_test{3,1},'-r');hold on
plot(xaxis,S_test{3,2},'-b');hold on
plot(xaxis,S_test{3,3},'-g');hold on
plot(xaxis,S_test{3,4},'-m');hold on
plot(xaxis,S_test{3,5},'-c');hold on
plot(xaxis,S_test{3,6},'-k');hold on
plot(xaxis,S_test{3,7},'-y');hold on

figure(6) %O
plot(xaxis,S_test{4,1},'-r');hold on
plot(xaxis,S_test{4,2},'-b');hold on
plot(xaxis,S_test{4,3},'-g');hold on
plot(xaxis,S_test{4,4},'-m');hold on
plot(xaxis,S_test{4,5},'-c');hold on
plot(xaxis,S_test{4,6},'-k');hold on
plot(xaxis,S_test{4,7},'-y');hold on

figure(7) %N
plot(xaxis,S_test{5,1},'-r');hold on
plot(xaxis,S_test{5,2},'-b');hold on
plot(xaxis,S_test{5,3},'-g');hold on
plot(xaxis,S_test{5,4},'-m');hold on
plot(xaxis,S_test{5,5},'-c');hold on
plot(xaxis,S_test{5,6},'-k');hold on
plot(xaxis,S_test{5,7},'-y');hold on

figure(8) %C
plot(xaxis,S_test{6,1},'-r');hold on
plot(xaxis,S_test{6,2},'-b');hold on
plot(xaxis,S_test{6,3},'-g');hold on
plot(xaxis,S_test{6,4},'-m');hold on
plot(xaxis,S_test{6,5},'-c');hold on
plot(xaxis,S_test{6,6},'-k');hold on
plot(xaxis,S_test{6,7},'-y');hold on

figure(9) %ALL
for i=1:7
    for j=1:7
        plot(xaxis,S_test{i,j});hold on
    end
end
