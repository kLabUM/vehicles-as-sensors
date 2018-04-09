function [qx,qy,qz] = drawFigure(M,radarTSeries,dSamp)


F = scatteredInterpolant(M(:,1), M(:,2), M(:,3));
[qx, qy] = meshgrid(linspace(min(M(:,1)), max(M(:,1)), ceil(size(radarTSeries,1)/dSamp)), ...
                    linspace(min(M(:,2)), max(M(:,2)), ceil(size(radarTSeries,2)/dSamp)));
qz = F(qx, qy);
fig2 = figure('position',[100 100 1200 450],'Color',[1 1 1]);
surfc(qx, qy, qz,'EdgeColor','none','FaceAlpha',0.8);
colormap('jet');hold on;
set(gca,'XTick',[min(min(qx)) max(max(qx))])
set(gca,'YTick',[min(min(qy)) max(max(qy))])
% set(gca,'ZTick',[0 max(max(qz))+0.01])
set(gca,'ZTick',[0 1+0.01])
set(gca,'XTickLabel',{sprintf('%.2f',-84.10),sprintf('%.2f', -83.30)}) % longitude
set(gca,'YTickLabel',{sprintf('%.2f',42.10),sprintf('%.2f',42.40)});    % latitude
% set(gca,'ZTickLabel',{'0',sprintf('%f',max(max(qz)))});
set(gca,'ZTickLabel',{'0',sprintf('%f',1)});

xlabel('longitude');ylabel('latitude');
set(gca,'fontsize',16);
colorbar;
caxis([0 1]);
grid on;
view(0,90);
%view(134,80);
% axis([-0.02+min(min(qx)) 0.02+max(max(qx)) -0.02+min(min(qy)) 0.02+max(max(qy)) 0 max(max(qz))+0.01])
axis([-0.02+min(min(qx)) 0.02+max(max(qx)) -0.02+min(min(qy)) 0.02+max(max(qy)) 0 1+0.01])