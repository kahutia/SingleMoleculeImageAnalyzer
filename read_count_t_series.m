if ~exist('rootpath','var')
    rootpath='C:\exp00 before washing\ch03 100pM DNA-bt 100pM TFEB 0hr Nuclear 532 0 min\LeonResult.20170910 17h38m03s';
end

[fname, rootpath]=uigetfile([rootpath '\*.*']);

mol=load([rootpath '\' fname]);

%% get basic parameters
del_t=0.5;  % unit: sec
frmax=max(mol(:,4));
tmax=frmax*del_t;
t_vct=(1:frmax)*del_t;


%% remove initial frames
cutoff_fr=0;
valid_id=mol(:,3)>cutoff_fr;
mol_selected=mol(valid_id,:);


%% remove short-lived molecules
% cutoff_t=1;
% valid_id=(mol_selected(:,4)-mol_selected(:,3))>cutoff_t;
% mol_selected=mol_selected(valid_id,:);

%% 
bind=mol_selected(:,3);
diss=mol_selected(:,4);

%% calc life time
lifetime=(diss-bind)*del_t;
x_bin_size=0.5;
x_bin=x_bin_size/2:x_bin_size:max(lifetime);

lifetime=lifetime(lifetime~=0);
life_hist=[];
[life_hist(:,2),life_hist(:,1)]=hist(lifetime,x_bin);

plot(life_hist(:,1),life_hist(:,2));
xlabel('Lifetime (s)');
ylabel('Count');

k_diss=1/mean(lifetime);

%% single exp fit
fo1 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[Inf,Inf],...
               'StartPoint',[500 0.8]);
ft1 = fittype('a*exp(-k*x)','options',fo1);
[fcurve1,gof1] = fit(life_hist(:,1),life_hist(:,2),ft1);

figure(101)
plot(life_hist(:,1),life_hist(:,2),'o');hold on;
plot(fcurve1,'r');hold off;
xlabel('Lifetime (s)');
ylabel('Count');


%% double exp fit
% fo2 = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,0,0,0],...
%                'Upper',[Inf,Inf,Inf,Inf],...
%                'StartPoint',[500, 100, 0.8, 10]);
% ft2 = fittype('a1*exp(-k1*x) + a2*exp(-k2*x) ','options',fo2);
% [fcurve2,gof2] = fit(life_hist(:,1),life_hist(:,2),ft2);
% 
% 
% figure(102)
% plot(life_hist(:,1),life_hist(:,2),'o');hold on;
% plot(fcurve2,'b');hold off;


%% binding and dissociation rate
N_bind=[];
N_diss=[];
for fri=1:frmax
    N_bind(fri,1)=sum(bind==fri);
    N_diss(fri,1)=sum(diss==fri);
end

k_bind=mean(N_bind(cutoff_fr+1:end))/del_t;


figure(103)
plot(t_vct,N_bind);hold on;
plot(t_vct,N_diss);hold off;
xlabel('Time (s)');
ylabel('Count');
legend('Binding','Dissociation');





