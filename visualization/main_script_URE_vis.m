clear all
close all
% Main-Script for visualization of the results
%==========================================================================
% profile on
%==========================================================================
%                          Open Reference Results
%========================================================================== 
%open reference solution obtained by the RadSymS
fileID = "euler_output_raref_2d.h5";
Q = h5read(fileID,"/vis_mat");

tplot = Q(1,:);
xplot = Q(2,:);
vfin = Q(3,:);
pfin = Q(4,:);
% nt = length(tplot);
mt = length(xplot);
vplot = Q(5:mt+4,:);
pplot = Q(mt+5:2*mt+4,:);

 
%set flag_GL to one for self similar solutions
flag_GL = 1;
if flag_GL == 1
    
    %open radially symmetric self similar solution of ODE [Lai2019]
    fileID_GL = "GL_output_raref_2d.h5";
    Q_GL = h5read(fileID_GL,"/vis_mat");

    r_GL = Q_GL(1,:);
    v_GL = Q_GL(2,:);
    p_GL = Q_GL(3,:);
            
end
%==========================================================================
%                          Open Trixi.jl Results
%==========================================================================
fileID = "../data/self_similar_expansion_P4estMesh2D_8_11_14.h5";
t_trixi = h5read(fileID,"/time");
r_trixi = h5read(fileID,"/radius");
p_trixi = h5read(fileID,"/pressure");
v_trixi = h5read(fileID,"/velocity");

%scale velocity
v_sc_trixi = v_trixi./sqrt(1 + v_trixi.^2);

% %correct output mistake in Trixi.jl code for Ex. 5 in 3d - maybe outdated
% v_trixi = v_trixi*sqrt(2)/sqrt(3);

[vr_idx,vt_idx] = size(v_trixi);
[pr_idx,pt_idx] = size(p_trixi);

%Limits for Colorbar - can be adjusted
p_up = ceil(max(max(max(p_trixi)),max(max(pplot))));
v_min = min(min(min(v_sc_trixi)),min(min(vplot))); %maybe use also floor()
v_up = ceil(max(max(max(v_sc_trixi)),max(max(vplot))));
%==========================================================================
%                              Output & Plots
set(0,'DefaultTextInterpreter','Latex')
%==========================================================================
figure
hp(1) = subplot(1,2,1);
surf(t_trixi,r_trixi,p_trixi,'EdgeColor','interp')
ylim([0 max(r_trixi)])
hp(1).FontSize = 14;
xlabel('$t$','FontSize',18)
ylabel('radius $x$','FontSize',18)
title('$p$ - Trixi.jl')
view(2)

hp(2) = subplot(1,2,2);
surf(tplot,xplot,pplot,'EdgeColor','interp')
ylim([0 max(r_trixi)])
hp(2).FontSize = 14;
xlabel('$t$','FontSize',18)
ylabel('radius $x$','FontSize',18)
title('$p$ - Radially Symmetric Scheme')
view(2)

%specify global colormap for both plots
set(hp(1), 'Colormap', grb_colormap(200), 'Clim', [0,p_up])
set(hp(2), 'Colormap', grb_colormap(200), 'Clim', [0,p_up])

cp = colorbar();

% set(hp(1),'Units','Points','Position',[50 50 200 300])
% set(hp(2),'Units','Points','Position',[290 50 200 300])
% set(cp,'Units','Points','Position',[505 50 35 300])


figure
hv(1) = subplot(1,2,1);
surf(t_trixi,r_trixi,v_sc_trixi,'EdgeColor','interp')
hv(1).FontSize = 40;
ylim([0 max(r_trixi)])
% colorbar
xlabel('$t$','FontSize',40)
ylabel('radius $x$','FontSize',40)
title('$v$ - Trixi.jl')
view(2)
   
hv(2) = subplot(1,2,2);
surf(tplot,xplot,vplot,'EdgeColor','interp')
hv(2).FontSize = 40;
ylim([0 max(r_trixi)])
xlabel('$t$','FontSize',40)
ylabel('radius $x$','FontSize',40)
title('$v$ - Radially Symmetric Scheme')
view(2)

%specify global colormap
set(hv(1), 'Colormap', grb_colormap(200), 'Clim', [v_min,v_up])
set(hv(2), 'Colormap', grb_colormap(200), 'Clim', [v_min,v_up])

cv = colorbar();
% set(hv(1),'Position',[0.04 0.07 0.425 0.89])
% set(hv(2),'Position',[0.51 0.07 0.425 0.89])
% set(cv,'Position',[0.945 0.07 0.025 0.89])


fig_v = figure;
hold on
plot(r_trixi,v_sc_trixi(:,vt_idx),'o',xplot,vfin,'.-')
fontsize(fig_v,30,"points")
xlim([0 max(xplot)])

if flag_GL == 1

    plot(r_GL,v_GL)
    legend('Trixi.jl','Radially Symmetric Scheme','ODE Solution','Location','NorthEast','FontSize', 30)

else
    
    legend('Trixi.jl','Radially Symmetric Scheme','Location','NorthEast','FontSize', 30)
    
end
hold off
xlabel('radius $x$','FontSize', 30)
ylabel('$v(t_{end},x)$','FontSize', 30)

fig_p = figure;
hold on
plot(r_trixi,p_trixi(:,vt_idx),'o',xplot,pfin,'.-')
fontsize(fig_p,30,"points")
xlim([0 max(xplot)])

if flag_GL == 1

    plot(r_GL,p_GL)
    legend('Trixi.jl','Radially Symmetric Scheme','ODE Solution','Location','NorthEast','FontSize', 30)

else
    
    legend('Trixi.jl','Radially Symmetric Scheme','Location','NorthEast','FontSize', 30)
    
end
hold off
xlabel('radius $x$','FontSize', 30)
ylabel('$p(t_{end},x)$','FontSize', 30)