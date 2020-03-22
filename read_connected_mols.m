%% input: .frames file from Trotsky run in connected molecules mode

if ~exist('rootpath','var')
    rootpath='F:\ptnPaint\data\20190325 SYNZIP 5-6\long movies\LeonResult.20190327114407';
end

[fname, rootpath]=uigetfile([rootpath '\*.frames'],'Choose .frames');

mol=load([rootpath '\' fname]);

%% get basic parameters
del_t=.1;  % unit: sec
fr_avr_factor=10;
frmax=max(mol(:));
tmax=frmax*del_t;
t_vct=(1:frmax)*del_t*fr_avr_factor;

%% binding and dissociation rate
N_bind=[];
N_diss=[];

fr_bind=mol(:,1);
fr_diss=mol(:,2);

t_tag=0;
for fri=1:fr_avr_factor:frmax-1
    t_tag=t_tag+1;
    N_bind(t_tag,1)=sum(fr_bind==fri);
    N_diss(t_tag,1)=sum(fr_diss==fri);
end

% k_bind=mean(N_bind(cutoff_fr+1:end))/del_t*fr_avr_factor;

t_vct=(1:fr_avr_factor:frmax-1)*del_t;
figure(103)
plot(t_vct,N_bind);hold on;
plot(t_vct,N_diss);hold off;
xlabel('Time (s)');
ylabel('Count');
legend('Binding','Dissociation');


%% For the NEWly binding molecules
cutoff_fr=3;
valid_id=mol(:,1)>cutoff_fr;
mol_selected=mol(valid_id,:);

% calc life time
fr_bind=mol_selected(:,1);
fr_diss=mol_selected(:,2);

lifetime=(fr_diss-fr_bind)*del_t;
lifetime=lifetime(lifetime~=0);

t_bin_size=del_t*fr_avr_factor;
% t_bin=t_bin_size/2:t_bin_size:max(lifetime)/10;
t_bin=t_bin_size/2:t_bin_size:50;

life_hist=[];
[life_hist(:,2),life_hist(:,1)]=hist(lifetime,t_bin);

k_diss=1/mean(lifetime);
fhd101=figure(101);clf
fhd101.Position(3)=300;
fhd101.Position(4)=300;

if 0
    % single exp fit
    x_data=life_hist(3:end,1);
    y_data=life_hist(3:end,2);

    lifetime_cutoff_min=1; % unit: frames
    x_data=x_data(lifetime_cutoff_min:end-1);
    y_data=y_data(lifetime_cutoff_min:end-1);

    fo1 = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0],...
                   'Upper',[Inf,Inf],...
                   'StartPoint',[500 0.8]);
    ft1 = fittype('a*exp(-x/tau)','options',fo1);
    [fcurve1,~] = fit(x_data,y_data,ft1);

    
    plot(x_data,y_data,'o');hold on;
    plot(fcurve1,'r');hold off;

    text(fcurve1.tau,fcurve1.a*2/3,['tau=' num2str(fcurve1.tau,'%.1f') 's'],'FontSize',12)
else
    % double exp fit
    x_data=life_hist(:,1);
    y_data=life_hist(:,2);

    lifetime_cutoff_min=1; % unit: frames
    x_data=x_data(lifetime_cutoff_min:end-1);
    y_data=y_data(lifetime_cutoff_min:end-1);

    fo1 = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,0,0],...
                   'Upper',[Inf,Inf,Inf,Inf],...
                   'StartPoint',[max(y_data) max(y_data)/2 1.5 20]);
    ft1 = fittype('a1*exp(-x/tau1) + a2*exp(-x/tau2)','options',fo1);

    [fcurve1,fff] = fit(x_data,y_data,ft1);

    plot(x_data,y_data,'o');hold on;
    plot(fcurve1,'r');hold off;
    
%     text(fcurve1.tau1,fcurve1.a1*2/3,['tau=' num2str(fcurve1.tau1,'%.1f') 's'],'FontSize',12)
%     text(fcurve1.tau2,fcurve1.a2*2/3,['tau=' num2str(fcurve1.tau2,'%.1f') 's'],'FontSize',12)
    tot_pop=fcurve1.a1+fcurve1.a2;
	text(fcurve1.tau1,fcurve1.a1*2/3,...
        ['tau=' num2str(fcurve1.tau1,'%.1f') 's' '(' num2str(fcurve1.a1/tot_pop*100,'%.1f') '%)' ],'FontSize',12)
    text(fcurve1.tau2,fcurve1.a2*2/3,...
        ['tau=' num2str(fcurve1.tau2,'%.1f') 's' '(' num2str(fcurve1.a2/tot_pop*100,'%.1f') '%)' ],'FontSize',12)
end
set(gca,'FontSize',15,'Position',[.3 .2 .65 .65]);
legend off;
xlabel('Lifetime (s)');
ylabel('Count');
xlim([0 30]);
title('newly comes');


%% collect initial frames 
cutoff_fr=1;
valid_id=mol(:,1)==cutoff_fr;
mol_selected=mol(valid_id,:);

% calc life time
fr_bind=mol_selected(:,1);
fr_diss=mol_selected(:,2);

lifetime=(fr_diss-fr_bind)*del_t;
lifetime=lifetime(lifetime~=0);

t_bin_size=del_t*fr_avr_factor;
% t_bin=t_bin_size/2:t_bin_size:max(lifetime)/10;
t_bin=t_bin_size/2:t_bin_size:200;

life_hist=[];
[life_hist(:,2),life_hist(:,1)]=hist(lifetime,t_bin);

k_diss=1/mean(lifetime);
fhd102=figure(102);clf;
fhd102.Position(3)=300;
fhd102.Position(4)=300;
if 1
    % single exp fit
    x_data=life_hist(3:end,1);
    y_data=life_hist(3:end,2);

    lifetime_cutoff_min=1; % unit: frames
    x_data=x_data(lifetime_cutoff_min:end-1);
    y_data=y_data(lifetime_cutoff_min:end-1);

    fo1 = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0],...
                   'Upper',[Inf,Inf],...
                   'StartPoint',[500 0.8]);
    ft1 = fittype('a*exp(-x/tau)','options',fo1);
    [fcurve1,~] = fit(x_data,y_data,ft1);

    plot(x_data,y_data,'o');hold on;
    plot(fcurve1,'r');hold off;  

    text(fcurve1.tau,fcurve1.a*2/3,['tau=' num2str(fcurve1.tau,'%.1f') 's'],'FontSize',12)
else
    % double exp fit
    x_data=life_hist(:,1);
    y_data=life_hist(:,2);

    lifetime_cutoff_min=1; % unit: frames
    x_data=x_data(lifetime_cutoff_min:end-1);
    y_data=y_data(lifetime_cutoff_min:end-1);

    fo1 = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,0,0],...
                   'Upper',[Inf,Inf,Inf,Inf],...
                   'StartPoint',[500 500 0.8 8]);
    ft1 = fittype('a1*exp(-x/tau1) + a2*exp(-x/tau2)','options',fo1);
    [fcurve1,~] = fit(x_data,y_data,ft1);

    plot(x_data,y_data,'o');hold on;
    plot(fcurve1,'r');hold off;

    text(fcurve1.tau1*2,fcurve1.a1,['tau=' num2str(fcurve1.tau1,'%.1f') 's'],'FontSize',12)
    text(fcurve1.tau2*2,fcurve1.a2,['tau=' num2str(fcurve1.tau2,'%.1f') 's'],'FontSize',12)
end

set(gca,'FontSize',15,'Position',[.3 .2 .65 .65]);
legend off;
xlabel('Lifetime (s)');
ylabel('Count');
xlim([0 30]);
title('initially bound');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Get N mol per frame
fname1=strrep(fname,'.frames','_leon_Nmol_all_frames.txt');

N_mol=load([rootpath '\' fname1]);

[t_max,~]=size(N_mol);
t_vct=(1:t_max)*del_t*fr_avr_factor;


%% single exp fit
fo1 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,-Inf],...
               'Upper',[Inf,Inf,Inf],...
               'StartPoint',[1000,1000, 150]);
ft1 = fittype('-a*exp(-x/tau)+offset','options',fo1);
[fcurve1,gof1] = fit(t_vct',N_mol(:,1),ft1);

figure(201)
plot(t_vct,N_mol(:,1),'o');hold on;
plot(fcurve1,'r');hold off;
xlabel('Time (s)');
ylabel('# mol. per image');
% ylim([0 500]);
% text(fcurve1.tau/5,fcurve1.a+fcurve1.offset,['tau=' num2str(fcurve1.tau,'%.0f') 's'],'FontSize',12)
text(fcurve1.tau/5,-fcurve1.a/2+fcurve1.offset,['tau=' num2str(fcurve1.tau,'%.0f') 's'],'FontSize',12)
legend off;
set(gca,'FontSize',15);
% %% two exp fit
% fo1 = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,0,0,0,-Inf],...
%                'Upper',[Inf,Inf,Inf,Inf,Inf],...
%                'StartPoint',[130 1,130 10,170]);
% ft1 = fittype('a1*exp(-x/tau1)+a2*exp(-x/tau2)+offset','options',fo1);
% [fcurve1,gof1] = fit(t_vct',N_mol(:,2),ft1);
% 
% figure(201)
% plot(t_vct,N_mol(:,2),'o');hold on;
% plot(fcurve1,'r');hold off;
% xlabel('Time (s)');
% ylabel('Count');
% text(fcurve1.tau1*5,fcurve1.a1/2+fcurve1.offset,num2str(fcurve1.tau1),'FontSize',12)
% text(fcurve1.tau2*5,fcurve1.a2/2+fcurve1.offset,num2str(fcurve1.tau2),'FontSize',12)