function [fit_dissociation_2exp, gof_2edx] = doubleExpFit(cum_duration, duration, metadata, fit_dissociation)

frames = metadata.totalFrames;
dt = metadata.deltaT;
tot_time = frames*dt;

ft = fittype( 'a*exp(-b*x)+c*exp(-d*x)+e', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.StartPoint = [cum_duration(1)/2 fit_dissociation.b/3 cum_duration(1)/2 fit_dissociation.b*3 cum_duration(end)];
opts.Lower = [0 0 0 0 0];
[fit_dissociation_2exp, gof_dissociation] = fit(duration, cum_duration, ft, opts);
adjR2_dissociation = gof_dissociation.adjrsquare;
temp_dissociation = confint(fit_dissociation_2exp);
conf_int_Koff1 = temp_dissociation(:,2);
conf_int_frac1 = temp_dissociation(:,1);
conf_int_Koff2 = temp_dissociation(:,4);
conf_int_frac2 = temp_dissociation(:,3);
rev_value = fit_dissociation_2exp.a+fit_dissociation_2exp.c;
irr_fraction = fit_dissociation_2exp.e/(rev_value +fit_dissociation_2exp.e)*100;
% Calculation of error propagation from x=c/(a+c): dx =
% x*sqrt((dx/da)^2+(dx/dc)^2) -> dx = 1.41*ac/(a+c)^2
irr_error = irr_fraction*rev_value*sqrt(2)/(rev_value+fit_dissociation_2exp.c)^2;
figure
subplot(4,1,1:3);
plot(fit_dissociation_2exp,duration,cum_duration,'.');
title('Dissociation');
%legend(gca,'off');
ylabel('Bound particles');
xlim([0 tot_time/2]);
txt_dissociation = {['Fit function: ' num2str(fit_dissociation_2exp.a,3)...
    '*exp(-' num2str(fit_dissociation_2exp.b,'%.2e') '*x)+' num2str(fit_dissociation_2exp.c,3)... 
    '*exp(-' num2str(fit_dissociation_2exp.d,'%.2e') '*x)+' num2str(fit_dissociation_2exp.e,3)],...
    ['K_{off,1} = ' num2str(fit_dissociation_2exp.b,'%.2e') char(177) num2str(fit_dissociation_2exp.b-conf_int_Koff1(1),'%.2e') ' s^-1'],...
    ['K_{off,2} = ' num2str(fit_dissociation_2exp.d,'%.2e') char(177) num2str(fit_dissociation_2exp.d-conf_int_Koff2(1),'%.2e') ' s^-1'],...
    ['Irreversible fraction (IF) = ' num2str(irr_fraction,2) char(177) num2str(irr_error,2) '%'],...
    ['Adjusted R squared: ' num2str(adjR2_dissociation,3)]};
text(tot_time/2*.95,(cum_duration(1)-cum_duration(end))*.95+cum_duration(end),...
    txt_dissociation,'HorizontalAlignment','right','VerticalAlignment','top')

if nargin>3
    hold on
    plot(linspace(duration(1), duration(end),1000), fit_dissociation...
        (linspace(duration(1), duration(end),1000)), 'Color', [75 122 71]/255);
    legend(gca,'Data','Double Exp','Single Exp','Location','southeast','Box','off');
    hold off
end

subplot(4,1,4);
res_2exp = cum_duration-fit_dissociation_2exp(duration);
plot(duration, res_2exp,'.-r');
hold on
res_exp = cum_duration-fit_dissociation(duration);
plot(duration, res_exp,'.-','Color', [75 122 71]/255);
plot([duration(1) duration(end)], [0 0],'k');
xlim([0 tot_time/2]);
ylim([min([res_2exp; res_exp])-1 max([res_2exp; res_exp])+1])
xlabel('Time since attachment (s)');
ylabel('residuals'),
drawnow

gof_2edx.dissociation = gof_dissociation;
gof_2edx.dissociation.koff_95conf = conf_int_Koff1';
gof_2edx.dissociation.frac1_95conf = conf_int_frac1';
gof_2edx.dissociation.koff_95conf = conf_int_Koff2';
gof_2edx.dissociation.frac2_95conf = conf_int_frac2';
gof_2edx.dissociation.irr_frac_95conf = [irr_fraction-irr_error irr_fraction+irr_error];