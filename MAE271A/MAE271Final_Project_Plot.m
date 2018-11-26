
%% compare estiamtion and true value
clc
close all
figure('Name','Estimation V.S True Value','rend','painters','pos',[10 10 1600 800]);
ax1_1 = subplot(3,1,1);
p1_1 = plot(ax1_1,ta,estimate_state(1,:),ta,p);
p1_1(1).Color = [0 0 128]./255;
p1_1(2).Color = [139 0 0]./255;
title(ax1_1,'Position');
ylabel(ax1_1,'position(m)');
xlabel(ax1_1, 'time(s)');
legend('Estimate','True Value');
set(ax1_1,'FontSize',14);

ax1_2 = subplot(3,1,2);
p1_2 = plot(ax1_2,ta,estimate_state(2,:),ta,v);
p1_2(1).Color = [0 0 128]./255;
p1_2(2).Color = [139 0 0]./255;
title(ax1_2,'Velocity');
ylabel(ax1_2,'velocity(m/s)');
xlabel(ax1_2, 'time(s)');
legend('Estimate','True Value');
set(ax1_2,'FontSize',14);

ax1_3 = subplot(3,1,3);
p1_3 = plot(ax1_3,ta,estimate_state(3,:),ta,b*ones(1,Na));
p1_3(1).Color = [0 0 128]./255;
p1_3(2).Color = [139 0 0]./255;
title(ax1_3,'Bias');
ylabel(ax1_3,'bias(m/s^2)');
xlabel(ax1_3, 'time(s)');
legend('Estimate','True Value');
set(ax1_3,'FontSize',14);

%% Plot Errors
%Error P
figure('Name','Error Estimate of P, V, Bias','rend','painters','pos',[10 10 1600 800]);
ax2_1 = subplot(3,1,1);
p2_1 = plot(ax2_1, ta, priori_E(1,:), tg, posteriori_E(1,:), tg, sqrt(abs(P_track(1,1:3*tsG/tsA:end))), tg, -sqrt(P_track(1,1:3*tsG/tsA:end)));
p2_1(1).Color = [0 0 128]./255;
p2_1(2).Color = [139 0 0]./255;
p2_1(3).Color = [34 139 34]./255;
p2_1(4).Color = [47 79 79]./255;
legend('Priori Error','Posteriori Error', 'one sigma bound', '-one sigma bound');
ylabel(ax2_1,'Error of Position');
xlabel(ax2_1, 'time(s)');
title(ax2_1,'Position Error');
set(ax2_1,'FontSize',14);
%Error V
ax2_2 = subplot(3,1,2);
p2_2 = plot(ax2_2, ta, priori_E(2,:), tg, posteriori_E(2,:), tg, sqrt(abs(P_track(2,2:3*tsG/tsA:end))), tg, -sqrt(P_track(2,2:3*tsG/tsA:end)));
p2_2(1).Color = [0 0 128]./255;
p2_2(2).Color = [139 0 0]./255;
p2_2(3).Color = [34 139 34]./255;
p2_2(4).Color = [47 79 79]./255;
legend('Priori Error','Posteriori Error', 'one sigma bound', '-one sigma bound');
ylabel(ax2_2,'Error of Velocity');
xlabel(ax2_2, 'time(s)');
title(ax2_2,'Velocity Error');
set(ax2_2,'FontSize',14);
%Error B
ax2_3 = subplot(3,1,3);
p2_3 = plot(ax2_3, ta, priori_E(3,:), tg, posteriori_E(3,:), tg, sqrt(abs(P_track(3,3:3*tsG/tsA:end))), tg, -sqrt(P_track(3,3:3*tsG/tsA:end)));
p2_3(1).Color = [0 0 128]./255;
p2_3(2).Color = [139 0 0]./255;
p2_3(3).Color = [34 139 34]./255;
p2_3(4).Color = [47 79 79]./255;
legend('Priori Error','Posteriori Error', 'one sigma bound', '-one sigma bound');
ylabel(ax2_3,'Error of Bias');
xlabel(ax2_3, 'time(s)');
title(ax2_3,'Bias Error');
set(ax2_3,'FontSize',14);

%% check ave_E = 0
figure('Name','Check average of el = 0','rend','painters','pos',[10 10 1600 800]);
ax3_1 = subplot(3,1,1);
p3_1 = plot(ax3_1, tg, e_ave(1,:),tg, zeros(1,Ng),'--');
p3_1(2).Color = [139 0 0]./255;
p3_1(1).Color = [0 0 128]./255;
legend('position e-avg','0')
ylabel(ax3_1,'position e-avg');
xlabel(ax3_1, 'time(s)');
title(ax3_1,'Position');
set(ax3_1,'FontSize',14);

ax3_2 = subplot(3,1,2);
p3_2 = plot(ax3_2, tg, e_ave(2,:),tg, zeros(1,Ng),'--');
p3_2(2).Color = [139 0 0]./255;
p3_2(1).Color = [0 0 128]./255;
legend('velocity e-avg','0')
ylabel(ax3_2,'velocity e-avg');
xlabel(ax3_2, 'time(s)');
title(ax3_2,'Velocity');
set(ax3_2,'FontSize',14);

ax3_3 = subplot(3,1,3);
p3_3 = plot(ax3_3, tg, e_ave(3,:),tg, zeros(1,Ng),'--');
p3_3(2).Color = [139 0 0]./255;
p3_3(1).Color = [0 0 128]./255;
legend('bias e-avg','0')
ylabel(ax3_3,'bias e-avg');
xlabel(ax3_3, 'time(s)');
title(ax3_3,'Bias');
set(ax3_3,'FontSize',14);

%% check P_avg - P = 0
figure('Name','Check P_ave - P = 0','rend','painters','pos',[10 10 1600 800]);
m = zeros(1,9);
for i = 1:9
    ax4 = subplot(3,3,i);
    set(ax4,'FontSize',14);
    if mod(i,3) == 0
        plot(tg, p_ave(ceil(i/3),3:3:end)-P_track(ceil(i/3),3:3*tsG/tsA:end),'color',[139 0 0]./255);
        m(i) = 3;
    else
        plot(tg, p_ave(ceil(i/3),mod(i,3):3:end)-P_track(ceil(i/3),mod(i,3):3*tsG/tsA:end),'color',[139 0 0]./255);
        m(i) = mod(i,3);
    end
    title(['Pave - P (',num2str(ceil(i/3)),',',num2str(m(i)),')']);
    xlabel('time(s)');
    ylabel(['Pave(',num2str(ceil(i/3)),',',num2str(m(i)),') - P(',num2str(ceil(i/3)),',',num2str(m(i)),')']); 
end

%% check orthogonality of the error
figure('Name','Check Orthogonality of the Error','rend','painters','pos',[10 10 1600 800]);
m = zeros(1,9);
for i = 1:9
    ax5 = subplot(3,3,i);
    set(ax5,'FontSize',14);
    if mod(i,3) == 0
        plot(tg, orth_ave(ceil(i/3),3:3:end),'color',[139 0 0]./255);
        m(i) = 3;
    else
        plot(tg, orth_ave(ceil(i/3),mod(i,3):3:end),'color',[139 0 0]./255);
        m(i) = mod(i,3);
    end
    title(['O(',num2str(ceil(i/3)),',',num2str(m(i)),')']);
    xlabel('time(s)');
    ylabel(['(',num2str(ceil(i/3)),',',num2str(m(i)),')']); 
end

%% check the residuals = 0
figure('Name','Check the Residuals = 0','rend','painters','pos',[10 10 1600 800]);
m = zeros(1,4);
for i = 1:4
    ax6 = subplot(2,2,i);
    set(ax6,'FontSize',14);
    if mod(i,2) == 0
        plot(tg, riXrm(ceil(i/2),2:2:end),'color',[139 0 0]./255);
        m(i) = 2;
    else
        plot(tg, riXrm(ceil(i/2),mod(i,2):2:end),'color',[139 0 0]./255);
        m(i) = mod(i,2);
    end
    title(['R(',num2str(ceil(i/2)),',',num2str(m(i)),')']);
    xlabel('time(s)');
    ylabel(['(',num2str(ceil(i/2)),',',num2str(m(i)),')']); 
end






