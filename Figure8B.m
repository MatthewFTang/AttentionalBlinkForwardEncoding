clearvars
SetUp
dataFolder= [dataFolder '/Exp2'];

participants =[1:4 6:15 17:25];


whiten=2
modelChannels=9
threshold =30
reconstruct=1
clear BigStore
partcount=1; 

title=(sprintf('BigIEM-whiten%2.0f-modelChannels%2.0f-threshold%2.0f-reconstruct%2.0f.mat',...
    whiten,modelChannels,threshold,reconstruct));
%%
for partN =1:23
    part=participants(partN);
    tic
    sprintf('Part %2.0f  whiten%2.0f modelChannels%2.0f threshold %2.0f reconstruct%2.0f',part,whiten,modelChannels,threshold,reconstruct)
    %% loads the cleaned EEG data
    cd([dataFolder '/EEG' ] )
    
    data=load('-mat',strcat(num2str(part),'-Cleaned.set'));
    EEG=data.EEG;
    
    cd([dataFolder '/Behaviour' ] )
        a=dir(strcat('Sub',num2str(part),'Run','*.csv'));
    X2=csvread(a.name);
    
    %%
    
    StoreBeh=zeros(length(X2),4);
    for trial =1:length(X2)
        T1pos= X2(trial,5);
        T1actual=mod(X2(trial,6+T1pos),180);
        T1report =X2(trial,27);
        
        
        T2pos= X2(trial,6);
        T2actual=mod(X2(trial,6+T2pos),180);
        T2report =X2(trial,28);
        
        StoreBeh(trial,:)=[ T1actual T1report T2actual T2report];
    end
    OrientationErrorT1= StoreBeh(:,1)-StoreBeh(:,2);
    OrientationErrorT1=mod((OrientationErrorT1+90),180)-90;
    OrientationErrorT2= StoreBeh(:,3)-StoreBeh(:,4);
    OrientationErrorT2=mod((OrientationErrorT2+90),180)-90;
    Lag = X2(:,4);
    
    
    CorrectT1=abs(OrientationErrorT1)<threshold;
    CorrectT2=abs(OrientationErrorT2)<threshold;
    
    
    
    
    Data =EEG.data;
    tio=EEG.times<0;
    Data = bsxfun(@minus,Data,mean(Data(:,tio,:),2));
    %
    OnsetTimes = 0:120:2300;
    StackedData = zeros(64,486,EEG.trials*20);
    StackedOri = zeros(EEG.trials*20,1);
    
    ItemCount= zeros(EEG.trials*20,1);
    TrialCount =zeros(EEG.trials*20,1);
    targetOrDistractorT1=zeros(EEG.trials*20,1);
    targetOrDistractorT2=zeros(EEG.trials*20,1);
    LagCount=zeros(EEG.trials*20,1);
    CorrectCountT1=zeros(EEG.trials*20,1);
    CorrectCountT2=zeros(EEG.trials*20,1);
    
    stackCount=1;
    
    
    for trial =1:EEG.trials
        T1pos= X2(trial,5);
        T2pos=X2(trial,6);
        for nitem =1:20
            Oris = X2(trial,6+nitem);
            ind=dsearchn(EEG.times',OnsetTimes(nitem));
            
            StackedData(:,:,stackCount)=Data(1:64,ind-128:ind+357,trial);
            StackedOri(stackCount)=Oris;
            ItemCount(stackCount)=nitem;
            TrialCount(stackCount)=trial;
            LagCount(stackCount)=Lag(trial);
            CorrectCountT1(stackCount)=CorrectT1(trial);
            CorrectCountT2(stackCount)=CorrectT2(trial);
            
            targetOrDistractorT1(stackCount)=T1pos==nitem;
            targetOrDistractorT2(stackCount)=T2pos==nitem;
            
            stackCount= stackCount+1;
        end
    end
    
    targetOrDistractor=(targetOrDistractorT1|targetOrDistractorT2);
    
    %%
    close all
    for type =1:2
        if type ==1
            Targets= find(targetOrDistractor==1); % train on targets
            trials =Targets;
        else
            Distractors= find(targetOrDistractor==0); % test on Distractors
            Distractors=Distractors(randperm(length(Distractors),length(Targets)));
            trials=Distractors;
        end
        
        
        Data= StackedData(:,:,trials);
        phi =StackedOri(trials);
        %%
        numN=size(Data,3);
        numT=size(Data,2);
        data=Data;
        angles =deg2rad(phi);
        angles=angles';
        n = length(angles);
        C1 = [cos(angles); sin(angles); ones(1,n)]; % chan x trials
        %%
        ORIG=zeros(size(data,1),size(data,2));
        parfor time = 1:size(data,2)
            dataTime  = squeeze(data(:,time,:));
            origTemp=zeros(1,size(data,1));
            for s = 1:size(data,1)
                thistrn = squeeze(dataTime(s,:)); % elec x trials
                W = C1'\thistrn';
                origTemp(s) = sqrt(W(1)^2 + W(2)^2);
            end
            ORIG(:,time)=origTemp;
        end
        
        nperm = 1000;
           %%
        RAND=zeros(size(data,1),size(data,2),nperm);
        parfor np = 1:nperm
            randTemp=zeros(size(data,1),size(data,2));
                fprintf('permuting sub :%2.0f | item : %2.0f |perm :%2.0f \n',part,nitem,np)

            for time = 1:size(data,2)
                dataTime  = squeeze(data(:,time,:));
                
                for s = 1:size(data,1)
                    thistrn = squeeze(dataTime(s,:)); % elec x trials
                    C2 = C1(:,randperm(n));
                    W = C2'\thistrn';
                    randTemp(s,time) = sqrt(W(1)^2 + W(2)^2);
                end
            end
            RAND(:,:,np)=randTemp;
        end
      %%
        Z=[];
        for s = 1:size(data,1)
            for t=1:size(data,2)
                perm = squeeze(RAND(s,t,:));
                actual = ORIG(s,t);
                p = sum(perm<actual)/nperm;
                Z(s,t) = norminv(p);
            end
        end
        Z = min(Z,10*ones(size(Z)));
        Z = max(Z,-10*ones(size(Z)));
        
        
        RandStore(type,:,:,:)=RAND;
        zStore(type,:,:)=Z;
        OrigStore(type,:,:)=ORIG;
        
        %%
    end
    cd([dataFolder ,'/Analysis'])
    save(sprintf('SelectivityTD %2.0f',part),'RandStore','zStore','OrigStore')
end
%%

for partN =1:23
    part=participants(partN);
sprintf('Loading :%2.0f\n',part)
    cd([dataFolder ,'/Analysis'])
    load(sprintf('SelectivityTD %2.0f',part))
zStore1(partN,:,:,:)=zStore;
origStore1(partN,:,:,:)=OrigStore;
end

    cd([dataFolder '/EEG' ] )

load('-mat','1-Cleaned.set')

EEG = pop_select(EEG,'channel',1:64);
%%
close all
ind= dsearchn(EEG.times',0);
NewTime = EEG.times(ind-128:ind+357);
c=1;
for t= 0:100:300
    
    tio= NewTime>t & NewTime <t+100;
    for type =1:2
        subplot(4,2,c)
        if type <3
        datM=squeeze(mean(zStore1(:,type,:,tio),[1 4]));
        
        else
                    datM1=collapse(OrigStore(:,1,:,tio),[1 4]);

                 datM2=collapse(OrigStore(:,2,:,tio),[1 4]);
   datM=datM1-datM2;
        end
        topoplot(datM,EEG.chanlocs,'maplimits',[-.5 2])
            c=c+1;

    end
    
end
%%
data=eeglab2fieldtrip(EEG,'preprocessing');
cfg= [];

e=ft_timelockanalysis(cfg,data);
e1{1}=e;
e1{2}=e;
cfg=[];
cfg.keepindividual ='yes';
tempGrand=ft_timelockgrandaverage(cfg,e1{:});
tempGrand.elec=data.elec;
ind= dsearchn(EEG.times',0);
newTime = EEG.times(ind-128:ind+357)/1000;
tempGrand.time = newTime;

%% correct trials
Cond1=tempGrand;
Cond1.individual=squeeze(zStore1(:,1,:,:));


Cond2=tempGrand;
Cond2.individual=squeeze(zStore1(:,2,:,:));
label='Correct trials';

cfg_neighb        = [];
cfg_neighb.method = 'distance';
neighbours        = ft_prepare_neighbours(cfg_neighb, Cond1);

cfg = [];
cfg.channel = 'EEG';
cfg.neighbours       = neighbours;



cfg.method  = 'montecarlo';
cfg.correctm ='cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsize';


cfg.statistic        = 'ft_statfun_depsamplesT';


cfg.minnbchan        = 1;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization =1500;
%cfg.spmversion = 'spm12';
subj = 23; % number of subjects
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

[stat] = ft_timelockstatistics(cfg,Cond1,Cond2);

%% gets the difference for plotting the effect
cfg= [];
cfg.parameter ='individual';
cfg.operation ='subtract';
Difference= ft_math(cfg,Cond1,Cond2);

differenceData=squeeze(mean(Difference.individual,1));
%% plots the results
close all
a= figure;
a.Name =label;
timesteps=.1;
Times2Plot=0:timesteps:.3; % the times that you want to plot in seconds

count=1;
for type =1:2
    if type ==1
        data = squeeze(mean(Cond1.individual,1));
    else
        data = squeeze(mean(Cond2.individual,1));
        
    end
    maplimits = [ -1 2];
    for time =Times2Plot % the time points that we want to plot
        subplot(3,length(Times2Plot)+1,count);
        hold on
        ind =dsearchn(stat.time', [time time+timesteps]'); % finds the
        %values in time that are of interested
        dat =(mean(data(:,ind(1):ind(2)),2));
        
        topoplot(dat,EEG.chanlocs,'maplimits',maplimits,'electrodes',...
            'off'); % plot significant
        hold on
        
        title(sprintf('%2.1f - %2.1f ',time,time+timesteps));
        count=count+1;
        
    end
            subplot(3,length(Times2Plot)+1,count);

    topoplot(dat,EEG.chanlocs,'maplimits',maplimits,'electrodes',...
        'off'); % plot significant
    colorbar
    count=count+1;
    
end

maplimits=[-1 1];
for time =Times2Plot % the time points that we want to plot
    subplot(3,length(Times2Plot)+1,count);
    hold on
    ind =dsearchn(stat.time', [time time+timesteps]'); % finds the
    %values in time that are of interested
    dat =(mean(differenceData(:,ind(1):ind(2)),2));
    marks = {'+','o'};
    sign = {'pos','neg'};
    for ii = 1:2 % cycle through pos / neg clusters
        if ~isempty(stat.([sign{ii} 'clusters']))
            sig_clust = find([stat.([sign{ii} 'clusters'])(:).prob] < stat.cfg.alpha);
            sig_clust_id = ismember(stat.([sign{ii} 'clusterslabelmat']), sig_clust);
            sig_clust_mark_any = any(sig_clust_id(:,ind(1):ind(2)), 2);
            
            topoplot(dat,EEG.chanlocs,'maplimits',maplimits,'electrodes',...
                'off','emarker2',{find(sig_clust_mark_any),marks{ii},'k',12,2}); % plot significant
            hold on
            
            title(sprintf('%2.1f - %2.1f ',time,time+timesteps));
            
        end
    end
    
    count=count+1;
end
subplot(3,length(Times2Plot)+1,count);
topoplot(dat,EEG.chanlocs,'maplimits',maplimits,'electrodes',...
    'off','emarker2',{find(sig_clust_mark_any),marks{ii},'k',6,1}); % plot significant
colorbar
colormap(linspecer(128))
cd(figFile)
set(gcf,'Color','w')
   set(gcf,    'Renderer','painters');
set(gcf,'Position',[   680   243   886   735])
print(['selectivityTD' '.eps'],'-depsc')
    