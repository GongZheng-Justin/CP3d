clear all
close all
clc

%% check x direction

for ii=1:3
    VarStr={'uu','vv','ww'};
    lcs=load (['.\StatOutOC0550\LCS_',VarStr{ii},'_x.txt']);

    %%
    Retau=lcs(1,1);
    yplus=lcs(1,2:end)';
    y=yplus/Retau;
    kx=lcs(2:end,1)';
    lx=kx.^(-1)*2*pi;
    coh=lcs(2:end,2:end)';
    clear lcs
    %%
    Font_Size=10;
    Font_Name='Times';

    figure
    % inner scale
    ax1=subplot(1,2,1);
    contourf(lx(2:end)*Retau,y*Retau,coh(:,2:end),[0.02,0.1:0.1:1],'LineColor',[0.5,0.5,0.5],'LineStyle','--');
    colorbar(ax1,'Location','northoutside','Position',[0.1,0.87,0.83,0.05],'Ticks',0:0.2:1);

    cgray=hot;
    cgray=flipud(cgray);
    colormap(ax1,cgray);

    set([ax1],'XScale','log','YScale','log','FontName',Font_Name,'FontSize',Font_Size,...
        'XMinorTick','on','YMinorTick','on',...
        'Xlim',[20,8*pi*Retau],'XTick',[100,1000,10000,100000],...
        'Ylim',[yplus(14),Retau],'YTick',[10,100,1000],...
        'TickLength',[0.02,0.008],'XColor','black','YColor','black','LineWidth',0.5,...
        'Position',[0.1,0.22,0.39,0.61]);
    xlabel(ax1,'\lambda_x^+','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    ylabel(ax1,'y^+','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');

    % outer scale
    ax2=subplot(1,2,2);
    contourf(lx(2:end),y,coh(:,2:end),[0.02,0.1:0.1:1],'LineColor',[0.5,0.5,0.5],'LineStyle','--');
    cgray=hot;
    cgray=flipud(cgray);
    colormap(ax2,cgray);
    set([ax2],'XScale','log','YScale','linear','FontName',Font_Name,'FontSize',Font_Size,...
        'XMinorTick','on','YMinorTick','on',...
        'Xlim',[20/Retau,8*pi],'XTick',[0.1,1,10],...
        'Ylim',[y(14),1],'YTick',0:0.2:1,...
        'TickLength',[0.02,0.008],'XColor','black','YColor','black','LineWidth',0.5,...
        'Position',[0.58,0.22,0.39,0.61]);
    xlabel(ax2,'\lambda_x/h','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    ylabel(ax2,'y/h','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');

    set(gcf,'Color','w','PaperPositionMode', 'auto','Units','centimeters','Position',[4,4,15,6]);
    print(gcf,[cd,'\ccf_lcs_x_',VarStr{ii},'.tiff'],'-dtiff','-r400');

    %%
    figure
    semilogx(lx(2:end),coh(end,2:end),'ko');
    set(gca,'XScale','log','YScale','linear','FontName',Font_Name,'FontSize',Font_Size,...
        'XMinorTick','on','YMinorTick','on',...
        'Xlim',[20/Retau,8*pi],'XTick',[0.1,1,10],...
        'TickLength',[0.02,0.008],'XColor','black','YColor','black','LineWidth',0.5);
    xlabel('\lambda_x/h','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    ylabel(['\gamma_{',VarStr{ii},'}'],'FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    set(gcf,'Color','w','PaperPositionMode', 'auto','Units','centimeters','Position',[4,4,8,6]);
    print(gcf,[cd,'\ccf_lcs_x_',VarStr{ii},'_yh1.tiff'],'-dtiff','-r400');

    clear all
end

clear all
close all
clc

%% check z direction

for ii=1:3
    VarStr={'uu','vv','ww'};
    lcs=load (['.\StatOutOC0550\LCS_',VarStr{ii},'_z.txt']);

    %%
    Retau=lcs(1,1);
    yplus=lcs(1,2:end)';
    y=yplus/Retau;
    kx=lcs(2:end,1)';
    lx=kx.^(-1)*2*pi;
    coh=lcs(2:end,2:end)';
    clear lcs
    %%
    Font_Size=10;
    Font_Name='Times';

    figure
    % inner scale
    ax1=subplot(1,2,1);
    contourf(lx(2:end)*Retau,y*Retau,coh(:,2:end),[0.02,0.1:0.1:1],'LineColor',[0.5,0.5,0.5],'LineStyle','--');
    colorbar(ax1,'Location','northoutside','Position',[0.1,0.87,0.83,0.05],'Ticks',0:0.2:1);

    cgray=hot;
    cgray=flipud(cgray);
    colormap(ax1,cgray);

    set([ax1],'XScale','log','YScale','log','FontName',Font_Name,'FontSize',Font_Size,...
        'XMinorTick','on','YMinorTick','on',...
        'Xlim',[10,3*pi*Retau],'XTick',[100,1000,10000,100000],...
        'Ylim',[yplus(14),Retau],'YTick',[10,100,1000],...
        'TickLength',[0.02,0.008],'XColor','black','YColor','black','LineWidth',0.5,...
        'Position',[0.1,0.22,0.39,0.61]);
    xlabel(ax1,'\lambda_z^+','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    ylabel(ax1,'y^+','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');

    % outer scale
    ax2=subplot(1,2,2);
    contourf(lx(2:end),y,coh(:,2:end),[0.02,0.1:0.1:1],'LineColor',[0.5,0.5,0.5],'LineStyle','--');
    cgray=hot;
    cgray=flipud(cgray);
    colormap(ax2,cgray);
    set([ax2],'XScale','log','YScale','linear','FontName',Font_Name,'FontSize',Font_Size,...
        'XMinorTick','on','YMinorTick','on',...
        'Xlim',[10/Retau,3*pi],'XTick',[0.1,1,10],...
        'Ylim',[y(14),1],'YTick',0:0.2:1,...
        'TickLength',[0.02,0.008],'XColor','black','YColor','black','LineWidth',0.5,...
        'Position',[0.58,0.22,0.39,0.61]);
    xlabel(ax2,'\lambda_z/h','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    ylabel(ax2,'y/h','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');

    set(gcf,'Color','w','PaperPositionMode', 'auto','Units','centimeters','Position',[4,4,15,6]);
    print(gcf,[cd,'\ccf_lcs_z_',VarStr{ii},'.tiff'],'-dtiff','-r400');

    %%
    figure
    semilogx(lx(2:end),coh(end,2:end),'ko');
    set(gca,'XScale','log','YScale','linear','FontName',Font_Name,'FontSize',Font_Size,...
        'XMinorTick','on','YMinorTick','on',...
        'Xlim',[5/Retau,3*pi],'XTick',[0.1,1,10],...
        'TickLength',[0.02,0.008],'XColor','black','YColor','black','LineWidth',0.5);
    xlabel('\lambda_z/h','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    ylabel(['\gamma_{',VarStr{ii},'}'],'FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    set(gcf,'Color','w','PaperPositionMode', 'auto','Units','centimeters','Position',[4,4,8,6]);
    print(gcf,[cd,'\ccf_lcs_z_',VarStr{ii},'_yh1.tiff'],'-dtiff','-r400');

    clear all
end



