clearvars;

SetUp
dataFolder= [dataFolder '/Exp2'];



cd([dataFolder '/Behaviour' ] )
partCount=1;
Cols= linspecer(2);
files =dir('*-AbGaborEEGv1.csv');
for part =[1:4 6:15 17:25]
    part
    
    
    cd([dataFolder '/Behaviour' ] )
    a=dir(strcat('Sub',num2str(part),'Run','*.csv'));
    X2=csvread(a.name);
    
    
    for lag =1:2
        trialC=1;
        clear Store
        
        
        for trial =1:length(X2)
            if X2(trial,4)==lag
                T1pos= X2(trial,5);
                T1actual=mod(X2(trial,6+T1pos),180);
                T1report =mod(X2(trial,27),180);
                
                T2pos= X2(trial,6);
                T2actual=mod(X2(trial,6+T2pos),180);
                T2report =mod(X2(trial,28),180);
                
                Store(trialC,:)=[ T1actual T1report T2actual T2report];
                trialC=trialC+1;
            end
        end
        
        
        %%
        for items =1:2
            %%
            
            if items ==1
                DifferenceOverall=(Store(:,1)-Store(:,2));
            else
                
                DifferenceOverall=(Store(:,3)-Store(:,4));
            end
            binSize=15;
            EDGES=-90:binSize:90;
            DifferenceOverall=mod(DifferenceOverall+90,180)-90;
            
            Difference= (histcounts(DifferenceOverall,EDGES));
            EDGES(end)=[];
            Difference=Difference/sum(Difference);
            
            StoreData(partCount,lag,items,:)=Difference;
            splitCorrectEEG(partCount,lag,items,:)=mean(abs(DifferenceOverall<30));
            
            X1=-90:90;
            lb = [0 -90 5 0];
            ub = [1 90 60 1];
            f=@(c,xdata)c(1).*exp(-0.5.*((xdata-c(2))./c(3)).^2) +c(4);
            x0=[.2 0 30 0];
            options = optimoptions('lsqcurvefit','Display','off');
            psy =  lsqcurvefit(f,x0,EDGES+(binSize/2),...
                Difference,lb,ub,options);
            Correct(partCount,lag,items) =mean(abs(DifferenceOverall)<(40));
            PsyStore(partCount,lag,items,:)=psy;
            
        end
    end
    partCount=partCount+1;
end


%% Panel A
close all
a=figure;
a.Color =[ 1 1 1];
name =[ 3 7];
for lag=1:2
    
    subplot(2,1,lag)
    for item =1:2
        dat =squeeze(StoreData(:,lag,item,:));
        datM=squeeze(mean(dat,1));
        datSE = squeeze(std(dat,1))/sqrt(size(StoreData,1));
        errorbar(EDGES,datM,datSE,'o','Color',Cols(item,:),'MarkerFaceColor'...
            ,Cols(item,:),'CapSize',0,'LineWidth',2);
        
        psy =  lsqcurvefit(f,x0,EDGES,...
            datM,lb,ub,options);
        hold on
        
        plot(X1,f(psy,X1),'Color',Cols(item,:),'MarkerFaceColor',Cols(item,:),'LineWidth',2)
        hold on
    end
    ylim([0.0 .3])
    yticks([0:.1:.3])
    xticks(-90:90:90)
    box off
    if lag ==1
        ylabel('Proportion of total responses','FontWeight','bold','FontAngle','italic')
    else
        %yticks([0 1])
    end
    if lag ==2
        xlabel({'Orientation error (\circ)'},'FontWeight','bold')
        
    end
    title(sprintf('Lag %2.0f',name(lag)));
    axes1=gca;
    set(axes1,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','off','YMinorTick','off'...
        ,'TickLength',[.015 .015],'Layer','bottom')
end
legend('T1','T1','T2','T2')

set(gcf,'Position',[   560   512   225   436])
cd(figFile)

print('TuningCurveExp2','-depsc')
%% Panel B
X= [ 3 7];
clc
a=figure;
a.Color=[1 1 1];
for var =1:4
    subplot(2,2,var)
    
    T1= PsyStore(:,:,1,var);
    T2= PsyStore(:,:,2,var);
    
    [H,P,CI,STATS] = ttest(T1,T2);
    for i=1:2
        fprintf('Var : %2.0f | items:%2.0f \n t(%2.0f) = %2.4f p = %2.4f\n',...
            var,i,STATS.df(i),STATS.tstat(i),P(i))
        
    end
    hold on
    
    
    T1mean=mean(T1,1);
    T1std=std(T1,1)/sqrt(length(T1));
    
    T2mean=mean(T2,1);
    T2std=std(T2,1)/sqrt(length(T1));
    
    
    errorbar(X,T1mean,T1std,'o-','Color',Cols(1,:),...
        'MarkerFaceColor',Cols(1,:),'LineWidth',2,'CapSize',0,...
        'MarkerSize',8)
    
    
    errorbar(X,T2mean,T2std,'o-','Color',Cols(2,:),...
        'MarkerFaceColor',Cols(2,:),'LineWidth',2','CapSize',0,...
        'MarkerSize',8)
    if var==1
        ylim([.07 .25])
        
        yticks(.1:.1:.3)
        ylabel('Amplitude (a.u.)','FontWeight','bold','FontAngle','italic')
    elseif var ==2
        ylim([-20 20])
        yticks(-60:20:30)
        
        ylabel('Center orientation (\circ)','FontWeight','bold','FontAngle','italic')
    elseif var ==3
        ylabel('Width (\circ)','FontWeight','bold','FontAngle','italic')
        ylim([10 40])
        yticks(-30:10:40)
        
    elseif var==4
        ylabel('Baseline (a.u.)','FontWeight','bold','FontAngle','italic')
        yticks(0:.02:.3)
        ylim([0.01 .06])
    end
    box off
    if var ==3 || var ==4
        xlabel('Lag','FontWeight','bold')
    end
    if var ==3
        legend('T1','T2')
    end
    if var ==1
        point = [ .08 .08];
        
    elseif var ==2
        point = [ -25 -25];
    elseif var ==3
        point = [ 12 12];
        
    else
        point = [0.02 0.02 ];
        
    end
    for i=1:length(P)
        
        if P(i)<.05/2
            plot([X(i) X(i)] ,[point],'k*')
        end
    end
    xticks(X)
    xlim([0.75 9.25])
    
    axes1=gca;
    set(axes1,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','off','YMinorTick','off'...
        ,'TickLength',[.015 .015],'Layer','bottom')
    
    
end
set(gcf,'Position',[   680   748   341   350])
cd(figFile)
print('TuningCurveExp2Fits','-depsc')
%%
var=3;
dat= (PsyStore(:,:,:,var));
dat=reshape(dat,size(dat,1),size(dat,2)*size(dat,3));
dat=round(dat*1000)/1000;
%% Panel C



Cols= linspecer(2);
X=[3 7];
for part =[1:4 6:15 17:25]
    
    
    cd([dataFolder '/Behaviour' ] )
    a=dir(strcat('Sub',num2str(part),'Run','*.csv'));
    X2=csvread(a.name);
    
    
    
    
    %%
    
    for lag =1:2
        trialC=1;
        clear Store
        
        if lag ==1
            name ='Lag 2';
        else
            name ='Lag 7';
            
        end
        for items =1:2
            clear Store
            for trial =1:length(X2)
                if X2(trial,4)==lag
                    if items ==1
                        ItemPosition= X2(trial,5);
                    else
                        ItemPosition= X2(trial,6);
                    end
                    ItemPosition=ItemPosition+6;
                    Orientations=X2(trial,ItemPosition-4:ItemPosition+7);
                    
                    
                    T1report =mod(X2(trial,27),180);
                    T2report =mod(X2(trial,28),180);
                    OrientationAverageT1= deg2rad(mod(X2(trial,X2(trial,5)+6),180));
                    OrientationAverageT2= deg2rad(mod(X2(trial,X2(trial,6)+6),180));
                    av_or=rad2deg(circ_mean(OrientationAverageT1,OrientationAverageT2));
                    Store(trialC,:)=[T1report T2report  Orientations ];
                    trialC=trialC+1;
                end
            end
            %%
            %             for report =1:2
            RegressionMatrix =exp(1i*(deg2rad(Store(:,3:end)))); % makes it a complex
            % number for easier regression
            if items ==1
                Response =exp(1i*(deg2rad(Store(:,1)))); %what orientation the participant
            else
                Response =exp(1i*(deg2rad(Store(:,2)))); %what orientation the participant
            end
            % reported
            Result=pinv(RegressionMatrix'*RegressionMatrix)*RegressionMatrix'*Response;
            % this does the regression
            order =randperm(size(RegressionMatrix,1));
            
            
            if items ==1
                Responserand =exp(1i*(deg2rad(Store(order,1)))); %what orientation the participant
            else
                Responserand =exp(1i*(deg2rad(Store(order,2)))); %what orientation the participant
            end
            Resultrand=pinv(RegressionMatrix'*RegressionMatrix)*RegressionMatrix'*Responserand;
            
            ResultStore(part,lag,items,:)=Result;
            ResultStorerand(part,lag,items,:)=Resultrand;
            
            %             end
        end
    end
end
%
count=1;
a=figure;
a.Color =[1 1 1];
X=-4:7;
lags =[3 7];
for lag  =1:2
    
    subplot(2,1,count)
    
    for items=1:2
        Temp =squeeze(ResultStore([1:1:end],lag,items,:));
        tempMean =abs(mean(Temp,1));
        
        tempSTD = abs(std(Temp,1)/sqrt(length(Temp)));
        errorbar(X,tempMean,tempSTD,'o-','Color',Cols(items,:),...
            'MarkerFaceColor',Cols(items,:),'CapSize',0,'LineWidth',2)
        % plot(tempMean)
        hold on
    end
    count=count+1;
    ylim([0 .45])
    yticks(0:.2:.4)
    xlabel('Item Position')
    ylabel('Regression weight')
    xlim([-5 8])
    if  lag ==5
        xlabel({'Item position relative to target' ;'(0 is target)'},'FontWeight','bold')
    end
    if lag ==3
        ylabel('Regression weight','FontWeight','bold')
    end
    title(sprintf('Lag %2.0f',lags(lag)))
    legend('T1','T2')
    box off
    
    axes1=gca;
    set(axes1,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','off','YMinorTick','off'...
        ,'TickLength',[.015 .015],'Layer','bottom')
end
set(gcf,'Position',[   560   122   444   826])
cd(figFile)
print('RegressionWeightsExp','-depsc')

