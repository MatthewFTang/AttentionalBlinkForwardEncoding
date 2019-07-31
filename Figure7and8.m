clearvars;
SetUp
dataFolder= [dataFolder '/Exp2'];

participants =[1:4 6:15 17:25];

AnalysisName='StackedtargetsNew';
whiten=2;
modelChannels=9;
threshold =30;
reconstruct=1;
clear BigStore
partcount=1;
for partN =1:2%length(participants)
    part=participants(partN);
    tic
    sprintf('Part %2.0f  whiten%2.0f modelChannels%2.0f threshold %2.0f reconstruct%2.0f',part,whiten,modelChannels,threshold,reconstruct)
    %% loads the cleaned EEG data
    cd([dataFolder '/EEG'] )
    data=load('-mat',strcat(num2str(part),'-Cleaned.set'));
    EEG=data.EEG;
    %% loads the behavioural data
    
    cd([dataFolder '/Behaviour' ] )
    a=dir(strcat('Sub',num2str(part),'Run','*.csv'));
    X2=csvread(a.name);
    
    
    %% baseline line correct EEG data
    Data =EEG.data;
    tio=EEG.times<0;
    Data = bsxfun(@minus,Data,mean(Data(:,tio,:),2));
    %% stacks the data into each presentation of the Gabors
    OnsetTimes = 0:120:2300;
    StackedData = zeros(64,486,EEG.trials*20);
    StackedOri = zeros(EEG.trials*20,1);
    
    ItemCount= zeros(EEG.trials*20,1);
    TrialCount =zeros(EEG.trials*20,1);
    
    stackCount=1;
    
    for nitem =1:20
        Oris = X2(:,6+nitem);
        ind=dsearchn(EEG.times',OnsetTimes(nitem));
        
        StackedData(:,:,stackCount:stackCount+EEG.trials-1)=Data(1:64,ind-128:ind+357,:);
        StackedOri(stackCount:stackCount+EEG.trials-1)=Oris;
        ItemCount(stackCount:stackCount+EEG.trials-1)=nitem;
        TrialCount(stackCount:stackCount+EEG.trials-1)=1:EEG.trials;
        
        stackCount= stackCount+EEG.trials;
    end
    
    
    
    %%
    Y=StackedData;
    if whiten ==1
        Y=zscore(Y,[],1);
    end
    Trials = 1:EEG.trials;
    
    phi =StackedOri;
    
    %% create the design matrix
    n_ori_chans= modelChannels;
    
    funType = @(xx,mu) (cosd(xx-mu)).^(n_ori_chans-mod(n_ori_chans,2));
    
    xx = linspace(1,180,180);
    basis_set = nan(180,n_ori_chans);
    chan_center = linspace(180/n_ori_chans,180,n_ori_chans);
    
    for cc = 1:n_ori_chans
        basis_set(:,cc) = funType(xx,chan_center(cc));
    end
    
    stim_mask = zeros(length(phi),length(xx));
    
    for tt = 1:size(stim_mask,1)  % loop over trials
        stim_mask(tt,phi(tt))=1;
        
    end
    
    % Generate design matrix
    design = (stim_mask*basis_set)';
    
    
    numN=size(Y,3);
    numT=size(Y,2);
    Xhat=nan(n_ori_chans,numT,numN,'single');
    
    Index =1:numN;
    Oris = StackedOri;
    
    phi=Oris ;
    
    %% makes the cross-validation folds have ~equal numbers of trials
    % for each orientation bin
    X =dsearchn(chan_center',Oris); %% 1.2. Cross-validation
    
    cfg = [];
    cfg.nFold = 20;
    folds = createFolds(cfg, X);
    
    %% does the IEM in one of two ways
    pim=1;
    if pim ==0
        parfor time = 1:numT
            XhatTemp = nan(n_ori_chans,numN);
            DataTime = squeeze(Y(:,time,:));
            for ii =1:length(folds)
                
                TestTrials =folds{ii};
                TrainTrials =find(~ismember(Index,TestTrials))';
                Y_train=DataTime(:,TrainTrials);
                Y_test= DataTime(:,TestTrials);
                DesignTrain =design(:,TrainTrials);
                
                W   = squeeze(Y_train*DesignTrain'*pinv(DesignTrain*DesignTrain')); %OLS
                XhatTemp(:,TestTrials)=(pinv(W'*W)*W'*Y_test);
            end
            Xhat(:,time,:)  =  XhatTemp; %multiply test data by weights from training data
            
        end
        
    else
        
        parfor time =1:numT
            
            DataTime = squeeze(Y(:,time,:));
            XhatTemp = nan(n_ori_chans,numN);
            for ii =1:length(folds)
                
                TestTrials =folds{ii};
                TrainTrials =find(~ismember(Index,TestTrials));
                Y_train=DataTime(:,TrainTrials);
                Y_test= DataTime(:,TestTrials);
                DesignTrain =design(:,TrainTrials);
                
                cfg = [];
                decoder=train_beamformerMT(cfg,DesignTrain,Y_train);
                XhatTemp(:,TestTrials) = decode_beamformer(cfg, decoder,Y_test);
            end
            Xhat(:,time,:)=XhatTemp;
        end
    end
    
    %% multiples the IEM results by the basis set to go from channels to orientations
    reconstruction=nan(size(Xhat,3),size(Xhat,2),size(basis_set,1),'single');
    for time =1:size(Xhat,2)
        reconstruction(:,time,:)=squeeze(Xhat(:,time,:))'*basis_set';
    end
    reconstruction=permute(reconstruction,[3 2 1]);
    
    %
    UniqueOri=unique(X);
    %% re-aligns all trials on the presented orientation
    
    reconstruction_algined = nan(size(reconstruction),'single');
    uniPhi=unique(phi);
    numC=length(uniPhi);
    center=90;
    
    for ic = 1:180
        shiftInd= mod(ic-center:ic+center-1,180)+1;
        reconstruction_algined( :,:,phi==ic)=...
            reconstruction(shiftInd,:, phi==ic);
    end
    
    %% temporally realigno targets
    X4= repmat(X2,20,1);
    T1pos= X4(:,5);
    Lag = X4(:,4);
    OnsetTimes = 0:120:2300;
    Start = 64;
    End = 256;
    %%
    Store=nan(size(reconstruction_algined,1),size(reconstruction_algined,2),12,length(X2));
    for trial =1:length(X2)
        T1=T1pos(trial);
        c=1;
        
        dat= reconstruction_algined(:,:,TrialCount ==trial);
        for ii=-3:9
            
            temp =dat(:,:,T1+ii);
            Store(:,:,c,trial)=temp;
            c=c+1;
        end
    end
    %% Gets the behavioural resuslts
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
    %%
    Lag = X2(:,4);
    
    
    CorrectT1=abs(OrientationErrorT1)<threshold;
    CorrectT2=abs(OrientationErrorT2)<threshold;
    
    %% breaks the results down by behaviour
    TempStore=zeros(2,2,size(Store,1),size(Store,2),size(Store,3));
    
    for lag =1:2
        CorrectCount=1;
        for correct =0:1
            trials  = lag==Lag&CorrectT1==1& CorrectT2==correct;
            TempStore(lag,CorrectCount,:,:,:)=mean(Store(:,:,:,trials),4);
            CorrectCount=CorrectCount+1;
        end
    end
    toc
    BigStore(partN,:,:,:,:,:)=TempStore;
    %%
    partcount=partcount+1;
end



%% fits Gaussians to the data for quantification
StoreBig=BigStore;
filtered= filtfast(StoreBig,5,[],'gaussian',3);
lb = [-10 -5 10 -3];
ub = [20   5  90 3];
f=@(c,xdata)c(1).*exp(-0.5.*((xdata-c(2))./c(3)).^2) +c(4);
x0=[.5 0 40 -.2];
options = optimoptions('lsqcurvefit','Display','off');
Orients =-89:90;
clear FitStore TempStore
for part =1:size(StoreBig,1)
    TempStore=zeros(2,2,size(StoreBig,5),size(StoreBig,6),4);
    for lag =1:2
        for correct=1:2
            
            for items =[4 7 11]
                
                parfor time =1:size(StoreBig,5)
                    
                    psy =  lsqcurvefit(f,x0,Orients',...
                        squeeze(filtered(part,lag,correct,:,time,items)),lb,ub,options);
                    
                    TempStore(lag,correct,time,items,:)=psy;
                end
            end
        end
        fprintf('Part : %2.0f | lag :%2.0f \n',part,lag)
    end
    FitStore(part,:,:,:,:,:)=TempStore;
end
%% Figure 6A
close all
for var=1:4
    a=figure;
    a.Color=[1 1 1];
    n=-3:9;
    ind= dsearchn(EEG.times',0);
    cols = linspecer(2);
    Newtime =  EEG.times(ind-128:ind+357);
    TIO = Newtime >-50 & Newtime <350; % the times to plot
    Newtime = Newtime(1:size(FitStore,4));
    c=1;
    LagReal=[3 7];
    a=figure;
    a.Color=[1 1 1];
    for lag =1:2
        if lag ==1
            b=0;
            items1 = [4 7 ];
        else
            b=13;
            items1 = [4 11];
            
        end
        for items=items1
            
            subplot(2,2,c)
            
            dat = squeeze(FitStore(:,lag,:,:,items,var));
            
            
            [datobs, datrnd] = cluster_test_helper(([squeeze(dat(:,2,:)-dat(:,1,:))])', 5000);
            [h_mem, p_mem, ~] = cluster_test(datobs,datrnd,0,0.05,0.05,'sum');
            
            
            
            for correct =1:2
                
                
                
                datM=squeeze(mean(dat(:,correct,:),1));
                datSE=squeeze(std(dat(:,correct,:),1))/sqrt(size(dat,1));
                
                hleng{correct}=shadedErrorBar(Newtime,datM,datSE,{'Color',cols(correct,:)},.5);
                hold on
                
            end
            
            title(sprintf('Lag :%2.0f \n item : %2.0f',LagReal(lag),n(items)))
            
            pclustu = unique(p_mem);
            npclust = nnz(pclustu < 0.05);
            for ipclust = 1:npclust % extract time range of each significant cluster and show in figure
                currind  = p_mem == pclustu(ipclust);
                SigTimes = Newtime(currind)
                sigInd= find(currind);
                fprintf('Item :%2.0f p = %2.4f\n',items,pclustu)
                plot([Newtime(sigInd(1))  Newtime(sigInd(end))], [-1 -1 ],'LineWidth',4,'Color','k')
                
            end
            
            box off
            if items ~=4
                yticks([])
            else
                if lag ==2
                    ylabel('Orientaiton selectivity (a.u.)')
                end
            end
            if lag ~=2
                
                xticks([])
            else
                xlabel('Time (ms')
                
            end
            
            c=c+1;
            if var ==1
                ylim([-1 4])
            elseif var ==2
                ylim([-20 20])
            elseif var ==3
                
                ylim([-10 80])
            else
                ylim([-2 2])
                
            end
            xlim([-50 350])
            plot([0 0],[Newtime(1) Newtime(end)],'k','LineWidth',1)
            plot([Newtime(1) Newtime(end)],[0 0],'k','LineWidth',1)
            
        end
    end
    set(gcf,'Renderer','Painters')
    set(gcf,'Position',[        1000         907         412         431])
    cd(figFile)
    print(['Figure6A-' num2str(var)],'-depsc')
    
end

%% Figure 6B
X=-89:90;
lb = [-10 -90 5 -3];
ub = [10 90 90 3];
f=@(c,xdata)c(1).*exp(-0.5.*((xdata-c(2))./c(3)).^2) +c(4);
x0=[.5 0 50 -.2];
Orients =-89:90;
options = optimoptions('lsqcurvefit','Display','off');

clear corr
c=1;
a=figure;
a.Color=[1 1 1];

SigTimes=Newtime >100 & Newtime <150;% times to average the data over
clear ax storeTemp
for lag =1:2;
    if lag ==1
        items =[7]
        
    else
        items =[11]
        
    end
    for correct =1:2
        ax(c)= subplot(2,1,c);
        dat = squeeze(mean(StoreBig(:,lag,correct,:,SigTimes,items),5));
        
        for p=1:size(dat,1)
            
            psy =  lsqcurvefit(f,x0,Orients',...
                squeeze(dat(p,:))',lb,ub,options);
            storeTemp(p,correct,:)=psy;
            
        end
        datM= mean(dat,1);
        datSE  =std(dat,1)/sqrt(size(dat,1));
        hleng{correct}=shadedErrorBar(X,datM,datSE,{'Color',cols(correct,:)},.5);
        
        hold on
        
        box off
        xlabel('Orientation (\circ)')
        ylabel('Orientation response (a.u.)')
        yticks(-1:2)
        xticks([-90:90:90])
    end
    c=c+1;
    for var =1:4
        
        t1= squeeze(storeTemp(:,1,var));
        t2= squeeze(storeTemp(:,2,var));
        mT1= mean(t1);
        mT2=mean(t2);
        sT1=std(t1);
        sT2=std(t2);
        r=corr(t1,t2);
        effect= (mT1-mT2)/(sqrt(sT1^2+sT2^2 -2*r*sT1*sT2));
        
        
        
        [H,P,CI,STATS] = ttest(t1,t2);
        fprintf('Items: %2.0f |var :%2.0f t %2.2f(%2.2f)p = %2.4f,d =%2.2f \n'...
            ,items,var,STATS.tstat,STATS.df, P,effect)
    end
    ylim([-1 2])
end
match_ylim(ax)
c=legend([hleng{1}.mainLine hleng{2}.mainLine ],'Incorrect','Correct');

set(gcf,'Position',[        1000         679         216         659])

cd(figFile)
print('Figure6B','-depsc')

%%
a=figure;
c=1;
X=-89:90;
tio = Newtime>100 & Newtime <150;
clear storeTemp
for lag =1:2
    for items =1:13
        
        for correct=1:2
            dat = squeeze(mean(StoreBig(:,lag,correct,:,tio,items),5));
            datM=squeeze(mean(dat,1));
            datSe = squeeze(std(dat,1))/sqrt(size(dat,1));
            for p=1:size(dat,1)
                
                psy =  lsqcurvefit(f,x0,Orients',...
                    squeeze(dat(p,:))',lb,ub,options);
                storeTemp(p,lag,items,correct,:)=psy;
            end
        end
        c=c+1;
    end
    
end
%% Figure 7A
a=figure;
a.Color=[1 1 1];
c=1;
clc
Cols2 =linspecer(4);
for var=1
    Cols = linspecer(2);
    for lag =1
        if lag ==1
            Dis=[1 2 3 5 6 8:13];
            
            
        else
            Dis=[1 2 3 5:10 12:13];
            
            
        end
        subplot(2,1,c)
        hold on
        plot([0 0],[-100 100],'--k')
        
        plot([3 3],[-100 100],':k')
        
        
        for correct =1:2
            dat = squeeze(storeTemp(:,lag,:,correct,var));
            datM=squeeze(mean(dat,1));
            datSe = squeeze(std(dat,1))/sqrt(size(dat,1));
            
            Distractor= squeeze(storeTemp(:,lag,Dis,:,var));
            Distractor=mean(Distractor,3);
            DistractorM= mean(mean(Distractor,2),1);
            DistractorSE= std(mean(Distractor,2),1)/sqrt(size(Distractor,1));
            %DistractorM=repmat(DistractorM,1,length(n));
            %DistractorSE=repmat(DistractorSE,1,length(n));
            hold on
            
            fh=fill([n(1),n(1),n(end),n(end)],[DistractorM+DistractorSE,...
                DistractorM-DistractorSE,DistractorM-DistractorSE,DistractorM+DistractorSE],[.5 .5 .5],'EdgeColor','none');
            %   h1=fill([n(1),n(1),n(end),n(end)],[DistractorM-DistractorSE,...
            %   DistractorM-DistractorSE,DistractorM+DistractorSE,DistractorM+DistractorSE],[0 0 0],'EdgeColor','none');
            set(fh,'facealpha',.2);
            
            
            % h{correct}=shadedErrorBar(n,DistractorM,DistractorSE,{'Color',[.5 .5 .5]},.5);
            errorbar(n,datM,datSe,'-o','Color',Cols(correct,:),'MarkerFaceColor'...
                ,Cols(correct,:),'CapSize',0,'LineWidth',2,'MarkerSize',8)
            hold on
            
            %% do stats comparison for targets vs distractors
            
            
            for target=1:2
                if target ==1
                    CompareDistractors=[2 3 5 6];
                    CompareTargets=4;
                else
                    CompareDistractors=[5 6 8 9];
                    CompareTargets=7;
                end
                
                distractorsCompare = mean(dat(:,CompareDistractors),2);
                
                targetsCompare = mean(dat(:,CompareTargets),2);
                
                t1= distractorsCompare;
                t2=targetsCompare;
                mT1= mean(t1);
                mT2=mean(t2);
                sT1=std(t1);
                sT2=std(t2);
                r=corr(t1,t2);
                effect= (mT1-mT2)/(sqrt(sT1^2+sT2^2 -2*r*sT1*sT2));
                
                
                
                [H,P,CI,STATS] =ttest(targetsCompare,distractorsCompare);
                fprintf('TargetsVsDistractors \n Var :%2.0f | correct :%2.0f | target :%2.0f | t= %2.2f, p=%2.4f,d=%2.2f\n',...
                    var,correct,target,STATS.tstat,  P,effect)
            end
            
            %% do stats comparison for distractors vs target+1 distractor
            
            
            for target=1:2
                if target ==1
                    CompareDistractors=[1:3 6:13];
                    CompareTargetPLus1=5;
                else
                    CompareDistractors=[1:7 9:13 ];
                    CompareTargetPLus1=8;
                end
                
                distractorsCompare = mean(dat(:,CompareDistractors),2);
                
                targetsCompare = mean(dat(:,CompareTargetPLus1),2);
                t1= distractorsCompare;
                
                t2=targetsCompare;
                mT1= mean(t1);
                mT2=mean(t2);
                sT1=std(t1);
                sT2=std(t2);
                r=corr(t1,t2);
                effect= (mT1-mT2)/(sqrt(sT1^2+sT2^2 -2*r*sT1*sT2));
                
                [H,P,CI,STATS] =ttest(targetsCompare,distractorsCompare);
                fprintf('Targets+T1 \n Var :%2.0f | correct :%2.0f | target :%2.0f | t= %2.2f, p=%2.4f,d=%2.2f\n',...
                    var,correct,target,STATS.tstat,  P,effect)
            end
            
            
        end
        
        
        box off
        yticks(-1:5)
        ylim([0.5 5.5])
        ylabel('Orientation response (a.u.)')
        
        
        
        xticks(-3:1:13)
    end
    c=c+1;
end
legend('T1','T2','','Incorrect','','Correct');

cd(figFile)
print('Figure7A','-depsc')
