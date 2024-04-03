clear all
close all
clc

%%
for ii=1:4

    yrefind=14;
    VarStr={'uu','vv','ww','pp'};

    I=load(['LCS_',VarStr{ii},'_theta_x.txt']);

    %%
    Retau=I(1,1);
    yplus=I(1,2:end)';
    y=yplus/Retau;
    kx=I(3:end,1)';
    lx=kx.^(-1)*2*pi;
    theta=I(3:end,2:end)';

    %%

    Font_Size=10;
    Font_Name='Times';

    figure
    % inner scale
    ax1=subplot(1,2,1);
    contourf(lx*Retau,y*Retau,theta,0:2:30,'LineColor',[0.5,0.5,0.5],'LineStyle','--');
    colorbar(ax1,'Location','northoutside','Position',[0.1,0.87,0.83,0.05],'Ticks',0:2:30);

    cgray=hot;
    cgray=flipud(cgray);
    colormap(ax1,cgray);

    set([ax1],'XScale','log','YScale','log','FontName',Font_Name,'FontSize',Font_Size,...
        'XMinorTick','on','YMinorTick','on',...
        'Xlim',[20,8*pi*Retau],'XTick',[100,1000,10000,100000],...
        'Ylim',[yplus(yrefind),Retau],'YTick',[10,100,1000],...
        'TickLength',[0.02,0.008],'XColor','black','YColor','black','LineWidth',0.5,...
        'Position',[0.1,0.22,0.39,0.61]);
    xlabel(ax1,'\lambda_x^+','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    ylabel(ax1,'y^+','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');

    % outer scale
    ax2=subplot(1,2,2);
    contourf(lx,y,theta,0:2:30,'LineColor',[0.5,0.5,0.5],'LineStyle','--');
    cgray=hot;
    cgray=flipud(cgray);
    colormap(ax2,cgray);
    set([ax2],'XScale','log','YScale','linear','FontName',Font_Name,'FontSize',Font_Size,...
        'XMinorTick','on','YMinorTick','on',...
        'Xlim',[20/Retau,8*pi],'XTick',[0.1,1,10],...
        'Ylim',[y(yrefind),1],'YTick',0:0.2:1,...
        'TickLength',[0.02,0.008],'XColor','black','YColor','black','LineWidth',0.5,...
        'Position',[0.58,0.22,0.39,0.61]);
    xlabel(ax2,'\lambda_x/h','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    ylabel(ax2,'y/h','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');

    set(gcf,'Color','w','PaperPositionMode', 'auto','Units','centimeters','Position',[4,4,15,6]);
    print(gcf,[cd,'\ccf_theta_x_',VarStr{ii},'.tiff'],'-dtiff','-r400');

    %%
    yh05ind=find(abs(y-0.5)==min(abs(y-0.5)));
    figure
    semilogx(lx,theta(yh05ind,:),'ko-');
    hold on
    semilogx(lx,zeros(length(lx(:)),1))
    set(gca,'XScale','log','YScale','linear','FontName',Font_Name,'FontSize',Font_Size,...
        'XMinorTick','on','YMinorTick','on',...
        'Xlim',[20/Retau,8*pi],'XTick',[0.1,1,10],...
        'ylim',[-90,90],'ytick',-90:45:90,...
        'TickLength',[0.02,0.008],'XColor','black','YColor','black','LineWidth',0.5);
    xlabel('\lambda_x/h','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    ylabel(['\theta_{',VarStr{ii},'}'],'FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    set(gcf,'Color','w','PaperPositionMode', 'auto','Units','centimeters','Position',[4,4,8,6]);
    print(gcf,[cd,'\ccf_theta_x_',VarStr{ii},'_yh05.tiff'],'-dtiff','-r400');

    clear all
end

clear all
close all
clc

%%
for ii=1:4

    yrefind=14;
    VarStr={'uu','vv','ww','pp'};

    I=load(['LCS_',VarStr{ii},'_theta_z.txt']);

    %%
    Retau=I(1,1);
    yplus=I(1,2:end)';
    y=yplus/Retau;
    kx=I(3:end,1)';
    lx=kx.^(-1)*2*pi;
    theta=I(3:end,2:end)';

    %%

    Font_Size=10;
    Font_Name='Times';

    figure
    % inner scale
    ax1=subplot(1,2,1);
    contourf(lx*Retau,y*Retau,theta,0:2:30,'LineColor',[0.5,0.5,0.5],'LineStyle','--');
    colorbar(ax1,'Location','northoutside','Position',[0.1,0.87,0.83,0.05],'Ticks',0:2:30);

    cgray=hot;
    cgray=flipud(cgray);
    colormap(ax1,cgray);

    set([ax1],'XScale','log','YScale','log','FontName',Font_Name,'FontSize',Font_Size,...
        'XMinorTick','on','YMinorTick','on',...
        'Xlim',[10,3*pi*Retau],'XTick',[100,1000,10000,100000],...
        'Ylim',[yplus(yrefind),Retau],'YTick',[10,100,1000],...
        'TickLength',[0.02,0.008],'XColor','black','YColor','black','LineWidth',0.5,...
        'Position',[0.1,0.22,0.39,0.61]);
    xlabel(ax1,'\lambda_z^+','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    ylabel(ax1,'y^+','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');

    % outer scale
    ax2=subplot(1,2,2);
    contourf(lx,y,theta,0:2:30,'LineColor',[0.5,0.5,0.5],'LineStyle','--');
    cgray=hot;
    cgray=flipud(cgray);
    colormap(ax2,cgray);
    set([ax2],'XScale','log','YScale','linear','FontName',Font_Name,'FontSize',Font_Size,...
        'XMinorTick','on','YMinorTick','on',...
        'Xlim',[10/Retau,3*pi],'XTick',[0.1,1,10],...
        'Ylim',[y(yrefind),1],'YTick',0:0.2:1,...
        'TickLength',[0.02,0.008],'XColor','black','YColor','black','LineWidth',0.5,...
        'Position',[0.58,0.22,0.39,0.61]);
    xlabel(ax2,'\lambda_z/h','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    ylabel(ax2,'y/h','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');

    set(gcf,'Color','w','PaperPositionMode', 'auto','Units','centimeters','Position',[4,4,15,6]);
    print(gcf,[cd,'\ccf_theta_z_',VarStr{ii},'.tiff'],'-dtiff','-r400');

    %%
    yh05ind=find(abs(y-0.5)==min(abs(y-0.5)));
    figure
    semilogx(lx,theta(yh05ind,:),'ko-');
    hold on
    semilogx(lx,zeros(length(lx(:)),1))
    set(gca,'XScale','log','YScale','linear','FontName',Font_Name,'FontSize',Font_Size,...
        'XMinorTick','on','YMinorTick','on',...
        'Xlim',[10/Retau,3*pi],'XTick',[0.1,1,10],...
        'ylim',[-90,90],'ytick',-90:45:90,...
        'TickLength',[0.02,0.008],'XColor','black','YColor','black','LineWidth',0.5);
    xlabel('\lambda_z/h','FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    ylabel(['\theta_{',VarStr{ii},'}'],'FontName',Font_Name,'FontSize',Font_Size,'FontAngle','italic');
    set(gcf,'Color','w','PaperPositionMode', 'auto','Units','centimeters','Position',[4,4,8,6]);
    print(gcf,[cd,'\ccf_theta_z_',VarStr{ii},'_yh05.tiff'],'-dtiff','-r400');

    clear all
end
