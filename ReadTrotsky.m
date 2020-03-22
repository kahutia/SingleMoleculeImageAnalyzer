% clear all;

%% Get root directory
if ~exist('rootpath','var')
    rootpath='D:\';
end
if ~ischar(rootpath)
    rootpath='D:\';
end
rootpath=uigetdir(rootpath);
save([userpath '\pgdata_ReadTrotsky.mat'],'rootpath');

    
%% define parameters
int_global_bg=0;    % data setting
y_int_min=105;  % plot setting


%% Get subdirectory info
tmp=dir(rootpath);
isdir_id=[0 0]; % the first two elements are for '.' and '..'
for sdi=3:length(tmp)
    isdir_id(sdi)=tmp(sdi).isdir;
end
subdirs_all=tmp(logical(isdir_id));
N_subdirs_all=length(subdirs_all);


%% select data set
data_set=1:N_subdirs_all;
% data_set=[1 2 3 4 5 6];
N_subdirs=length(data_set);
subdirs=subdirs_all(data_set);


%% Get Leon Results directories
for sdi=1:N_subdirs
    tmp=dir([rootpath '\' subdirs(sdi).name '\LeonResult*']);
    if ~isempty(tmp)
    	Leon_result_name{sdi}=tmp(end).name;
    else
        Leon_result_name{sdi}='';
    end
end

%% Get counts
N_mol=[];
N_mol_avr=[];
N_mol_avr_left=[];
N_mol_avr_right=[];

left_int=[];
right_int=[];
left_int_avr=[];
right_int_avr=[];



for sdi=1:N_subdirs
    co_local_mode(sdi)=0;
    if ~isempty(Leon_result_name{sdi})
        tmp_path=[rootpath '\' subdirs(sdi).name '\' Leon_result_name{sdi} '\'];
    else
        tmp_path=[rootpath '\' subdirs(sdi).name '\'];
    end
    
    mol_pos_files=dir([tmp_path '*Nmol.txt']);    
    N_movies=length(mol_pos_files);
    
    subdir_name{sdi}=subdirs(sdi).name;
    
    tmp_N_mol=[];
    tmp_left_int=[];
    tmp_right_int=[];
    tmp_N_mol_left=[];
    tmp_N_mol_right=[];

    for mvi=1:N_movies
        N_fid=fopen([tmp_path mol_pos_files(mvi).name],'r');
        
        kkk=fscanf(N_fid,'%s %s %s %s %t',4);
        tmp_left_int=[tmp_left_int fscanf(N_fid,'%g',1)];
        
        kkk=fscanf(N_fid,'%s %s %s %s %t',4);
        tmp_right_int=[tmp_right_int fscanf(N_fid,'%g',1)];
        
        kkk=fscanf(N_fid,'%s %s %t',2);
        
        if strcmp(kkk,'Co-localcount')
            co_local_mode(sdi)=1;
            kkk=fscanf(N_fid,'%s',1);
            tmp_N_mol=[tmp_N_mol fscanf(N_fid,'%g',1)];
            
            kkk=fscanf(N_fid,'%s %s %s %s',4);
            tmp_N_mol_left=[tmp_N_mol_left fscanf(N_fid,'%g',1)];
            
            kkk=fscanf(N_fid,'%s %s %s %s',4);
            tmp_N_mol_right=[tmp_N_mol_right fscanf(N_fid,'%g',1)];
        else
            tmp_N_mol=[tmp_N_mol fscanf(N_fid,'%g',1)];
        end
        
        fclose(N_fid);
    end
    
    tmp_left_int=tmp_left_int-int_global_bg;
    left_int{sdi}=tmp_left_int;
    left_int_avr(sdi,1)=mean(tmp_left_int);
    left_int_avr(sdi,2)=std(tmp_left_int);
    left_int_avr(sdi,3)=left_int_avr(sdi,2)/sqrt(N_movies);
    
    tmp_right_int=tmp_right_int-int_global_bg;
    right_int{sdi}=tmp_right_int;
    right_int_avr(sdi,1)=mean(tmp_right_int);
    right_int_avr(sdi,2)=std(tmp_right_int);
    right_int_avr(sdi,3)=right_int_avr(sdi,2)/sqrt(N_movies);
    
    N_mol{sdi}=tmp_N_mol;
    N_mol_avr(sdi,1)=mean(tmp_N_mol);
    N_mol_avr(sdi,2)=std(tmp_N_mol);
    N_mol_avr(sdi,3)=N_mol_avr(sdi,2)/sqrt(N_movies);
    
    
    N_mol_left{sdi}=tmp_N_mol_left;
    N_mol_avr_left(sdi,1)=mean(tmp_N_mol_left);
    N_mol_avr_left(sdi,2)=std(tmp_N_mol_left);
    N_mol_avr_left(sdi,3)=N_mol_avr_left(sdi,2)/sqrt(N_movies);
    
    N_mol_right{sdi}=tmp_N_mol_right;
    N_mol_avr_right(sdi,1)=mean(tmp_N_mol_right);
    N_mol_avr_right(sdi,2)=std(tmp_N_mol_right);
    N_mol_avr_right(sdi,3)=N_mol_avr_right(sdi,2)/sqrt(N_movies);
        

%             tmp_left_int=nan;
%             left_int{sdi}=nan;
%             left_int_avr(sdi,1)=nan;
%             left_int_avr(sdi,2)=nan;
%             left_int_avr(sdi,3)=nan;
% 
%             tmp_right_int=nan;
%             right_int{sdi}=nan;
%             right_int_avr(sdi,1)=nan;
%             right_int_avr(sdi,2)=nan;
%             right_int_avr(sdi,3)=nan;
% 
%             N_mol{sdi}=nan;
%             N_mol_avr(sdi,1)=nan;
%             N_mol_avr(sdi,2)=nan;
%             N_mol_avr(sdi,3)=nan;
%             
%                 N_mol_left{sdi}=nan;
%                 N_mol_avr_left(sdi,1)=nan;
%                 N_mol_avr_left(sdi,2)=nan;
%                 N_mol_avr_left(sdi,3)=nan;
%                 
%                 N_mol_right{sdi}=nan;
%                 N_mol_avr_right(sdi,1)=nan;
%                 N_mol_avr_right(sdi,2)=nan;
%                 N_mol_avr_right(sdi,3)=nan;

end



%% Get  Intensities of individual molecules

%%%%%%%% new scheme (Trotsky V6 and later)
Int_mol_all=[];
frames=10:15;

for sdi=1:N_subdirs   
    if ~isempty(Leon_result_name{sdi})
        tmp_path=[rootpath '\' subdirs(sdi).name '\' Leon_result_name{sdi} '\'];
    else
        tmp_path=[rootpath '\' subdirs(sdi).name '\'];
    end
    
%     tmp_path=[rootpath '\' subdirs(sdi).name '\' Leon_result_name{sdi} '\'];
%     if ~isempty(tmp_path)
        mol_staticI_files=dir([tmp_path '*.traces']);

        N_movies=length(mol_staticI_files);

        tmp_staticI=[];

        for mvi=1:N_movies
            %         mol_pos{mvi}=load([tmp_path mol_pos_files(mvi).name]);
            N_fid=fopen([tmp_path mol_staticI_files(mvi).name],'r');

            len=fread(N_fid,1,'int32');
            Ntraces=fread(N_fid,1,'int16');
            raw=fread(N_fid,Ntraces*len,'int16');
            fclose(N_fid);

            raw=reshape(raw,Ntraces,len);
            left_ch=raw(1:2:Ntraces,:);
    %         right_ch=raw(2:2:Ntraces,:);
            kkk=mean(left_ch(:,frames),2);

            tmp_staticI=[tmp_staticI kkk'];
        end

        staticI_all{sdi}=tmp_staticI';
%     end
end


%% %%%%%% old scheme (Trotsky V5 and earier)
% Int_mol_all=[];
% 
% 
% 
% for sdi=1:N_subdirs
%     tmp_path=[rootpath '\' subdirs(sdi).name '\' Leon_result_name{sdi} '\'];
%     mol_staticI_files=dir([tmp_path '*staticI.txt']);
%     
%     N_movies=length(mol_staticI_files);
%     
%     subdir_name{sdi}=subdirs(sdi).name;
%     tmp_staticI=[];
% 
%     for mvi=1:N_movies
%         %         mol_pos{mvi}=load([tmp_path mol_pos_files(mvi).name]);
%         N_fid=fopen([tmp_path mol_staticI_files(mvi).name],'r');
%         kkk=fscanf(N_fid,'%g');
%         tmp_staticI=[tmp_staticI kkk'];
%         fclose(N_fid);
%     end
%     
%     staticI_all{sdi}=tmp_staticI';
% end
% 


%% list the directories
for sdi=1:N_subdirs
    disp([num2str(sdi) ': ' subdir_name{sdi}])
end

%% Save results
fid=fopen([rootpath '\leon summary Nmol.txt'],'wt');
fprintf(fid,'%s\t%s\t%s\t%s\n','name','mean','std','se');
for sdi=1:N_subdirs
    fprintf(fid,'%s\t',subdirs(sdi).name);
    fprintf(fid,'%g\t%g\t%g\n',N_mol_avr(sdi,:));
end
fclose(fid);

fid=fopen([rootpath '\leon summary intensity.txt'],'wt');
fprintf(fid,'%s\t%s\t%s\t%s\n','name','mean','std','se');
for sdi=1:N_subdirs
    fprintf(fid,'%s\t',['left channel ' subdirs(sdi).name]);
    fprintf(fid,'%g\t%g\t%g\n',left_int_avr(sdi,:));
end
for sdi=1:N_subdirs
    fprintf(fid,'%s\t',['right channel ' subdirs(sdi).name]);
    fprintf(fid,'%g\t%g\t%g\n',right_int_avr(sdi,:));
end
fclose(fid);


%% Co-local counts Plot
if 1
    figure(11);clf;
    subplot(3,1,1);
    plothdl3=bar(1:N_subdirs,N_mol_avr(:,1),'FaceColor',[.8 .8 .8]);hold on;
    plothdl4=errorbar(1:N_subdirs,N_mol_avr(:,1),N_mol_avr(:,3),'.','Color','k');hold off;
    xlim([0.1 N_subdirs+0.9]);
    set(gca,'xTickLabel',num2cell(data_set));
    ylabel('Co-local Count');
    set(gca,'FontSize',18);

    subplot(3,1,2)
    plothdl1=bar(1:N_subdirs,N_mol_avr_left(:,1),'FaceColor',[0.3 .95 1]);hold on;
    plothdl11=errorbar(1:N_subdirs,N_mol_avr_left(:,1),N_mol_avr_left(:,3),'.','Color','k');hold off;
    xlim([0.1 N_subdirs+0.9]);
    set(gca,'xTickLabel',num2cell(data_set));
    ymax=ceil(1.2*max(N_mol_avr_left(:,1)));
    if isnan(ymax), ymax=1; end
    ylim([0 ymax]);
    ylabel('Green count');
    set(gca,'FontSize',18);

    subplot(3,1,3)
    plothdl1=bar(1:N_subdirs,N_mol_avr_right(:,1),'FaceColor',[0.95 .35 0]);hold on;
    plothdl11=errorbar(1:N_subdirs,N_mol_avr_right(:,1),N_mol_avr_right(:,3),'.','Color','k');hold off;
    xlim([0.1 N_subdirs+0.9]);
    set(gca,'xTickLabel',num2cell(data_set));
    ymax=ceil(1.2*max(N_mol_avr_left(:,1)));
    if isnan(ymax), ymax=1; end
    ylim([0 ymax]);
    ylabel('Red count');
    set(gca,'FontSize',18);
end

%% GFP Plot
figure(1);clf;
% subplot(1,2,1);
plothdl3=bar(1:N_subdirs,N_mol_avr(:,1),'FaceColor',[.9 .9 .9]);hold on;
plothdl3=plotSpread(N_mol,'distributionColors',[0.1 .6 .2]);hold on;
plothdl4=errorbar(1:N_subdirs,N_mol_avr(:,1),N_mol_avr(:,3),'o','MarkerSize',5,'Color','k');hold off;
set(gca,'xTickLabel',num2cell(data_set));
xlim([0.1 N_subdirs+0.9]);
ylabel('Count');
set(gca,'FontSize',18);

% subplot(1,2,2)
% plothdl1=bar(1:N_subdirs,left_int_avr(:,1),'FaceColor',[0.8 .95 .9]);hold on;
% plothdl3=plotSpread(left_int,'distributionColors',[0.1 .6 .2]);hold on;
% plothdl11=errorbar(1:N_subdirs,left_int_avr(:,1),left_int_avr(:,3),'.','Color','k');hold off;
% set(gca,'xTickLabel',num2cell(data_set));
% xlim([0.1 N_subdirs+0.9]);
% ylim([y_int_min ceil(1.2*max(left_int_avr(:,1)))]);
% ylabel('GFP intensity (a.u.)');
% set(gca,'FontSize',18);



%% RFP plot
figure(3);clf;
subplot(1,2,1);
plothdl3=bar(1:N_subdirs,N_mol_avr(:,1),'FaceColor',[.9 .9 .9]);hold on;
plothdl3=plotSpread(N_mol,'distributionColors',[.8 .1 .2]);hold on;
plothdl4=errorbar(1:N_subdirs,N_mol_avr(:,1),N_mol_avr(:,3),'o','MarkerSize',5,'Color','k');hold off;
set(gca,'xTickLabel',num2cell(data_set));
xlim([0.1 N_subdirs+0.9]);
ylabel('Count');
set(gca,'FontSize',18);

subplot(1,2,2)
plothdl1=bar(1:N_subdirs,right_int_avr(:,1),'FaceColor',[1 .6 .2]);hold on;
plothdl3=plotSpread(right_int,'distributionColors',[.8 .1 .2]);hold on;
plothdl11=errorbar(1:N_subdirs,right_int_avr(:,1),right_int_avr(:,3),'.','Color','k');hold off;
set(gca,'xTickLabel',num2cell(data_set));
xlim([0.1 N_subdirs+0.9]);
ylim([y_int_min ceil(1.2*max(right_int_avr(:,1)))]);
set(gca,'FontSize',18);
ylabel('RFP intensity (a.u.)');


%% special code for GFP time course 1
return;

data_set=[1:2];

figure(12137); clf; 

subplot(2,2,1); hold on;

for sdi=data_set
%     subplot(2,1,sdi); hold on;
	t_vct=((1:length(N_mol{sdi}))-1)*15;
	data=N_mol{sdi};
    
    fo1 = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,-Inf],...
                   'Upper',[Inf,Inf,Inf],...
                   'StartPoint',[1500 2 200]);
    ft1 = fittype('-a*exp(-x/tau)+offset','options',fo1);
    [fcurve1{sdi},~] = fit(t_vct',data',ft1);
    
    plot(t_vct,data,'o');
    plot(fcurve1{sdi},'r-');
%     text(20,fcurve1{sdi}.a+fcurve1{sdi}.offset,['tau=' num2str(fcurve1{sdi}.tau,'%.0f')]);
    text(20,-fcurve1{sdi}.a+fcurve1{sdi}.offset,['tau=' num2str(fcurve1{sdi}.tau,'%.0f')]);
    
    legend(subdirs(sdi).name);
%     ylim([0 500]);
    xlim([0 150]);
    ylabel('# molecules');
    xlabel('Time (s)');    
end

% legend off;
%%
subplot(2,2,2); hold on;
for sdi=data_set
	t_vct=((1:length(N_mol{sdi}))-1)*10;
    data=N_mol{sdi};

    tmp_out=data/max(data);
    plot(t_vct,tmp_out,'o');
end
ylabel('Population (norm.)');
xlabel('Time (s)');

%%
figure(12137); clf; 
% subplot(2,2,3); hold on;
data_set=[1 2];
for sdi=data_set
    subplot(2,1,sdi);hold on;
    t_vct=((1:length(left_int{sdi}))-1)*15;    
    
    fo1 = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,-Inf],...
                   'Upper',[Inf,Inf,Inf],...
                   'StartPoint',[1500 2 200]);
    ft1 = fittype('a*exp(-x/tau)+offset','options',fo1);
    [fcurve1{sdi},~] = fit(t_vct',left_int{sdi}',ft1);
    
    plot(t_vct,left_int{sdi},'o');
    plot(fcurve1{sdi},'r-');
    text(20,fcurve1{sdi}.a+fcurve1{sdi}.offset,['tau=' num2str(fcurve1{sdi}.tau,'%.0f')]);
    legend(subdirs(sdi).name);
end
ylabel('Intensity');
xlabel('Time (s)');
% legend off;

%%
subplot(2,2,4); hold on;
for sdi=data_set
    t_vct=((1:length(left_int{sdi}))-1)*10;
    tmp_out=(left_int{sdi}-min(left_int{sdi}))/(max(left_int{sdi}-min(left_int{sdi})));
    plot(t_vct,tmp_out,'o');
end
ylabel('Intensity (norm)');
xlabel('Time (s)');
legend(subdirs(data_set).name);




%% special code for GFP time course 2
t_vct=[(1:8)*15+30 8*15+30+(1:8*3)*60]/60;
t_vct=[(1:8*4)*30 8*30+(1:8*4)*60]/60;
cnt_per_movie=[];

for sdi=1:N_subdirs
    cnt_per_movie=[cnt_per_movie N_mol{sdi}];
end
tmp_N_movies=length(cnt_per_movie);

% uvid=[1, 10, 14, 19, 24, 26, 29];
uvid=[1, 5, 10, 15, 19, 20, 25, 30, 35, 40, 45, 49, 50, 55];
vid=true(1,tmp_N_movies);
for i=1:length(uvid)
    vid(uvid(i))=false;
end

t_vct=t_vct(vid);
cnt_per_movie=cnt_per_movie(vid);

% single exp fit
fo1 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,-Inf],...
               'Upper',[Inf,Inf,Inf],...
               'StartPoint',[1500 2 200]);
ft1 = fittype('a*exp(-x/tau)+offset','options',fo1);
[fcurve1,gof1] = fit(t_vct',cnt_per_movie',ft1);

figure(1112)
plot(t_vct,cnt_per_movie,'o');hold on;
plot(fcurve1,'r');hold off;
xlabel('Lifetime (s)');
ylabel('Count');


text(fcurve1.tau*5,fcurve1.a/2,num2str(fcurve1.tau),'FontSize',12)


%% GFP intensity sum
for i=1:N_subdirs
    staticI_sumall(i)=mean(staticI_all{i});
	staticI_stdall(i)=std(staticI_all{i});
end

% data_set=[2 4 6 8];

figure(12);clf;
plotSpread(staticI_all,'distributionColors',[.5 .7 .7]); hold on;
% plot(1:N_subdirs,staticI_sumall,'O','Color',[0 0 0]);
plothdl11=errorbar(1:N_subdirs,staticI_sumall,staticI_stdall,'o','Color',[0 .3 .0],'LineWidth',1);hold off;
set(gca,'xTickLabel',num2cell(data_set));
xlim([0.1 N_subdirs+0.9]);
ylim([0 35000]);
ylabel('Avr. spot intensity (a.u.)');
set(gca,'FontSize',18);
set(gca,'Box','on');
legend('individual spots','mean+std','Location','Northwest');legend boxoff;



%% intensity histograms
int_hist=[];
for i=1:N_subdirs
    xbin=(0:1000)*100;
    kkk=hist(staticI_all{i},xbin);  
    kkk=kkk(1:end-1);
%     int_hist(:,i)=kkk/max(kkk);
    int_hist(:,i)=kkk;
end

xbin=xbin(1:end-1);


figure(5);clf
subplot(5,1,1)
plot(xbin,int_hist(:,1),'Color',[1 0.3 0],'LineWidth',3);hold on;
plot(xbin,int_hist(:,2),'Color',[0 0.3 1],'LineWidth',3);hold off;
xlim([0 50000]);
legend({'G0','G12'},'FontSize',16);legend boxoff;

subplot(5,1,2)
plot(xbin,int_hist(:,3),'Color',[1 0.3 0],'LineWidth',3);hold on;
plot(xbin,int_hist(:,4),'Color',[0 0.3 1],'LineWidth',3);hold off;
xlim([0 3000]);
legend({'G0','G12'},'FontSize',16);legend boxoff;

subplot(5,1,3)
plot(xbin,int_hist(:,5),'Color',[1 0.3 0],'LineWidth',3);hold on;
plot(xbin,int_hist(:,6),'Color',[0 0.3 1],'LineWidth',3);hold off;
xlim([0 3000]);
legend({'G0','G12'},'FontSize',16);legend boxoff;

% subplot(5,1,4)
% plot(xbin,int_hist(:,7),'Color',[1 0.3 0],'LineWidth',3);hold on;
% plot(xbin,int_hist(:,8),'Color',[0 0.3 1],'LineWidth',3);hold off;
% xlim([0 3000]);
% legend({'G0','G12'},'FontSize',16);legend boxoff;

% subplot(5,1,5)
% plot(xbin,int_hist(:,9),'Color',[1 0.3 0],'LineWidth',3);hold on;
% plot(xbin,int_hist(:,10),'Color',[0 0.3 1],'LineWidth',3);hold off;
xlim([0 50000]);
legend({'G0','G12'},'FontSize',16);legend boxoff;

%% intensity histograms (special for two-layer ultra results)
int_hist=[];
xbin=0:50:10000;
xbin_plot=xbin(2:end-1)/100;

if N_subdirs<11
    color_scheme=[1 0 0; 1 .6 0; 1 1 0; .7 .8 0; 0 .8 .4; 0 0.6 1; .2 0 1; 0.2 0 .6; .8 0 .8; .8 0.5 .5];
else
    color_scheme=rand(N_subdirs,3);
end

figure(5);clf;hold on;
for i=1:N_subdirs
    kkk=hist(staticI_all{i},xbin);  
    kkk=kkk(2:end-1);
    int_hist(:,i)=kkk/max(kkk(:));
    plot(xbin_plot,int_hist(:,i),'Color',color_scheme(i,:),'LineWidth',3);
    legend_texts{i}=num2str(data_set(i));
end


x_lims=[0 10000]/100;
xlim(x_lims);
legend(legend_texts,'FontSize',12,'Location','NorthEast');legend boxoff;
xlabel('Spot intensity (a.u.)');
ylabel('Population');
set(gca,'FontSize',14);
box on;



%% plot multimers (special for two-layer ultra)
cutoffint=5000;
    
for i=1:N_subdirs
    mt_intsum(i)=sum(staticI_all{i}(staticI_all{i}>cutoffint)/100);
    mt_cnt(i)=sum(staticI_all{i}>cutoffint);
end

figure(512);clf;
subplot(2,1,1);
bar(mt_cnt);
xlim([0.1 N_subdirs+0.9]);
set(gca,'xTickLabel',num2cell(data_set));
ylabel(['Count (>' num2str(cutoffint/100) ')']);
subplot(2,1,2);
bar(mt_intsum);
xlim([0.1 N_subdirs+0.9]);
set(gca,'xTickLabel',num2cell(data_set));
ylabel(['Sopt int. sum (>' num2str(cutoffint/100) ')']);

%% intensity histograms (special for 5x5)


% subplot(5,1,1)
% plot(xbin,int_hist(:,1),'Color',[1 0.3 0],'LineWidth',3);hold on;
% plot(xbin,int_hist(:,2),'Color',[0 0.3 1],'LineWidth',3);hold off;
% xlim(x_lims);
% legend({'G0','G12'},'FontSize',16);legend boxoff;
% 
% subplot(5,1,2)
% plot(xbin,int_hist(:,3),'Color',[1 0.3 0],'LineWidth',3);hold on;
% plot(xbin,int_hist(:,4),'Color',[0 0.3 1],'LineWidth',3);hold off;
% xlim(x_lims);
% legend({'G0','G12'},'FontSize',16);legend boxoff;
% 
% subplot(5,1,3)
% plot(xbin,int_hist(:,5),'Color',[1 0.3 0],'LineWidth',3);hold on;
% plot(xbin,int_hist(:,6),'Color',[0 0.3 1],'LineWidth',3);hold off;
% xlim(x_lims);
% legend({'G0','G12'},'FontSize',16);legend boxoff;
% 
% subplot(5,1,4)
% plot(xbin,int_hist(:,7),'Color',[1 0.3 0],'LineWidth',3);hold on;
% plot(xbin,int_hist(:,8),'Color',[0 0.3 1],'LineWidth',3);hold off;
% xlim(x_lims);
% legend({'G0','G12'},'FontSize',16);legend boxoff;
% 
% subplot(5,1,5)
% plot(xbin,int_hist(:,9),'Color',[1 0.3 0],'LineWidth',3);hold on;
% plot(xbin,int_hist(:,10),'Color',[0 0.3 1],'LineWidth',3);hold off;
% xlim(x_lims);
% legend({'G0','G12'},'FontSize',16);legend boxoff;

