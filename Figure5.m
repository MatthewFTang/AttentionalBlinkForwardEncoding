clearvars;
SetUp
dataFolder= [dataFolder '/Exp2'];

partcount=1;

%%
participants =[1:4 6:15 17:25];

for part =participants
    %% loads the cleaned EEG data
    cd([dataFolder '/EEG'] )
    load('-mat',strcat(num2str(part),'-Cleaned.set'))
    
    %% loads the behavioural data
    
    
    cd([dataFolder '/Behaviour' ] )
    a=dir(strcat('Sub',num2str(part),'Run','*.csv'));
    X2=csvread(a.name);
    
    
    Data =EEG.data(1:64,:,:);
    Y = Data;
    Trials = 1:EEG.trials;
    %%
    Activation=nan(180,EEG.pnts,20,'single');
    for nitem =1:20
        tic
        fprintf('Part : %2.0f | item : %2.0f \n',part,nitem)
        Oris = X2(:,6+nitem);
        
        phi=Oris ;
        
        chans=0:20:160;
        
        X =dsearchn(chans',Oris);
        
        %% 1.2. Cross-validation
        cfg = [];
        cfg.nFold = 10;
        folds = createFolds(cfg, X);
        
        %% create the design matrix
        n_ori_chans= 9;
        
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
        %% two
       
            
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
        
        
        %%
        reconstruction=nan(size(Xhat,3),size(Xhat,2),size(basis_set,1),'single');
        for time =1:size(Xhat,2)
            reconstruction(:,time,:)=squeeze(Xhat(:,time,:))'*basis_set';
        end
        reconstruction=permute(reconstruction,[3 2 1]);
        
        %
        UniqueOri=unique(X);
        
        %% re-align the reconstruction
        center=90;
        close all
        reconstruction_algined = nan(size(reconstruction),'single');
        uniPhi=unique(phi);
        for ii=1:length(uniPhi)
            shift = center-uniPhi(ii);
            ind = phi==uniPhi(ii);
            reconstruction_algined(:,:,ind)=circshift(reconstruction(:,:,ind),shift,1);
        end
        
        %%
        toc
        close all
        Activation(:,:,nitem)=mean(reconstruction_algined,3);
        
        
    end
    StoreParts(partcount,:,:,:,:)=Activation;
    partcount=partcount+1;
end

%%


filtered= filtfast(StoreParts,3,[],'gaussian',4);


clear FitStore

StoreParts=double(StoreParts);
filtered=double(filtered);
lb = [-4 -90 5 -1];
ub = [4 90 90 1];
f=@(c,xdata)c(1).*exp(-0.5.*((xdata-c(2))./c(3)).^2) +c(4);
x0=[.2 0 90 -0.2];
options = optimoptions('lsqcurvefit','Display','off');
Orients=-89:90;
for part =1:size(StoreParts,1)
    parfor items =1:20
            TempStore=zeros(size(StoreParts,3),4);

        for time =1:size(StoreParts,3)            
            
            psy =  lsqcurvefit(f,x0,Orients,...
                squeeze(filtered(part,:,time,items)),lb,ub,options);
            
            TempStore(time,:)=psy;
        end
        fprintf('Part : %2.0f | item :%2.0f \n',part,items)
            FitStore(part,items,:,:,:)=TempStore;

    end
end
%



%% Panel C
close all
AverageSub = squeeze(mean(filtered,1));
ISI=120;
Times = 0:ISI:2300;
tio =EEG.times>-100 & EEG.times<2700;
Orients=-89:90;
a=figure;
a.Color=[1 1 1];
var=1;
%
clim=[-.17 .17];
c=1;
Nt =EEG.times(tio);
for ii=1:20
    %% plots the color map
    ax(c)=subplot(20,1,c);
    
    imagesc(EEG.times(tio),Orients,AverageSub(:,tio,ii),clim)
    hold on
    plot([Times(ii) Times(ii)],[Orients(1) Orients(end)],'k')
    
    
    plot([Nt(1) Nt(end)],[0 0],':k','LineWidth',1)
    yticks([])
    c=c+1;
    
    if ii~=20
        xticks([])
    else
        xticks(0:1000:2000)
        xlabel('Time (ms)')
    end
    if ii==10
        ylabel('Orientation (\circ)')
        
    end
    
    box off
end
colormap(linspecer(125))

set(gcf,'Render','Painters')
set(gcf,'Position',[731   437   368   715])




%% Panel D
c=1;
a=figure;
a.Color=[1 1 1];
L=75;
N=20;
ColorOris=linspace(0,200,N);
colors=[L*cosd(ColorOris) ;L*sind(ColorOris)] ;
colors=[repmat(L,1,N); colors];
cols=lab2rgb(colors','OutputType','uint8');


for ii=1:20
    step=.4;
    %% plot the fits over time with p results
    dat = squeeze(FitStore(:,ii,:,var));
    datM = squeeze(mean(dat,1));
    datSE = squeeze(std(dat))/sqrt(size(dat,1));
    hold on
    
    plot([EEG.times(1) EEG.times(end)],[ -(step*ii) -(step*ii) ],'Color',[.7 .7 .7],'LineWidth',.5)
    xlim([EEG.times(1) EEG.times(end)])
    if ii==20
        xticks(0:1000:2000)
    else
        xticks([])
    end
    hold on
    
    yticks([])
    shadedErrorBar(EEG.times,datM-(step*ii),datSE,{'Color',cols(ii,:),'LineWidth',2},.5)
    [datobs, datrnd] = cluster_test_helper(dat', 2000);
    [h_mem, p_mem, ~] = cluster_test(datobs,datrnd,0,0.05,0.05);
    
    Time =EEG.times;
    pclustu = unique(p_mem);
    npclust = nnz(pclustu < 0.05);
    for ipclust = 1:npclust % extract time range of each significant cluster and show in figure
        currind  = p_mem == pclustu(ipclust);
        sigInd= find(currind);
        plot([EEG.times(sigInd(1))  EEG.times(sigInd(end))], [ -(step*ii) -(step*ii) ],'LineWidth',2,'Color',cols(ii,:))
        
    end
    c=c+1;
    box off
end
%
set(gcf,'Render','Painters')
set(gcf,'Position',[1000         886         342         452])
%% Panel A and B


clear Aligned
for ii=1:20
    
    ind = dsearchn(EEG.times',(Times(ii)));
    
    Aligned(:,:,:,ii)=filtered(:,:,ind-100:ind+256,ii);
    
end


ind = dsearchn(EEG.times',0);
clims=[-.17 .17];
NewTime  = EEG.times(ind)-100:EEG.times(ind)+256;
close all
datM=collapse(Aligned,[1 4]);
subplot(1,2,1)
imagesc(-89:90,NewTime,datM',clims)
colormap(linspecer(128))
c1=colorbar;
c1.Label.String='Orientation response (a.u.)';
c1.Ticks=[-.1:.1:.1]
ylabel('Time (ms)')
xticks(-45:45:45)
yticks([0:100:200])
subplot(1,2,2)
cols=linspecer(4);
c=1;
clear h
for t=[ 50 100 150]
    
    ind= NewTime >t & NewTime <t+50;
    dat=collapse(Aligned(:,:,ind,:),[3 4]);
    datM= squeeze(mean(dat,1));
    
    datSE = squeeze(std(dat,1))/sqrt(size(dat,1));
    h{c}=shadedErrorBar(-89:90,datM,datSE,{'Color',cols(c,:),'LineWidth',2},.5)
    hold on
    c=c+1;
end
c=legend([h{1}.mainLine h{2}.mainLine h{3}.mainLine],'50-100','100-150','150-200');


xticks(-45:45:45)
yticks([-.1:.1:.1])
box off
xlabel('Orientation (\circ)')




