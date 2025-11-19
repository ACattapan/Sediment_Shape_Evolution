function Convergence_Figure_check(i,j,k,m,QuantilesP,TD_E_Quant,TD_M_Quant,I_Ka_Local)

FigNum=get(gcf,'Number');
% h1=figure(FigNum)
% hold on
% yyaxis left
% scatter(mu,C, 'filled')
% yyaxis right
% for t=1:k
%     scatter(mu,TD_Local{t}, 'filled')
% end
% plot(mu, TD_Local{I_Ka_Local(i,j,m)},'r','LineWidth',2)


h2=figure(FigNum+1);
hold on;
plot(TD_E_Quant(:,m),QuantilesP,'g','LineWidth',2)
% cdfplot(TD_E_Quant(:,m))
for t=1:k
    plot(TD_M_Quant(:,t),QuantilesP,'-ob')
    % cdfplot(TD_M_Quant(:,t))
end
plot(TD_M_Quant(:,I_Ka_Local(i,j,m)),QuantilesP,'-or','LineWidth',2)
xlim([0.95*min([min(TD_E_Quant(:,m)) min(TD_M_Quant(:,I_Ka_Local(i,j,m)))]) 1.05*max(TD_M_Quant(:,I_Ka_Local(i,j,m)))])
xlabel('Travel Distance [m]')
ylabel('F(x) [m^{-1}]')
end
