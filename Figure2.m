clearvars;

SetUp.m
dataFolder= [dataFolder '/Exp1'];

cd(dataFolder);
directory = dir('*.csv');
Cols= linspecer(2);
for part = 1:length(directory)
    X2=csvread(directory(part).name);
    
    for lag =1:5
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
            
            StoreData(part,lag,items,:)=Difference;
            
            
            X1=-90:90;
            lb = [0 -90 5 0];
            ub = [1 90 60 1];
            f=@(c,xdata)c(1).*exp(-0.5.*((xdata-c(2))./c(3)).^2) +c(4);
            x0=[.2 0 30 0];
            options = optimoptions('lsqcurvefit','Display','off');
            psy =  lsqcurvefit(f,x0,EDGES+(binSize/2),...
                Difference,lb,ub,options);
            PsyStore(part,lag,items,:)=psy;
            
        end
    end
end
%% Panel A
close all
a=figure;
a.Color =[ 1 1 1];
name =[ 1 2 3 5 7];
for lag=1:5
    
    subplot(5,1,lag)
    for item =1:2
        tempM =squeeze(mean(StoreData(:,lag,item,:),1));
        tempSe =squeeze(std(StoreData(:,lag,item,:),1))/sqrt(size(StoreData,1));
        
        errorbar(EDGES,tempM,tempSe,'o','Color',...
            Cols(item,:),'MarkerFaceColor',Cols(item,:),...
            'LineWidth',2,'CapSize',0);
        hold on
        psy =  lsqcurvefit(f,x0,EDGES,...
            tempM',lb,ub,options);
        hold on
        
        plot(X1,f(psy,X1),'Color',Cols(item,:),'MarkerFaceColor',Cols(item,:),'LineWidth',2)
        hold on
    end
   ylim([0 .25])
    yticks([0 .2])
    
    if lag ==1
        ylabel('Proportion of total responses','FontWeight','bold','FontAngle','italic')
    else
    end
    if lag ==5
        xlabel({'Orientation error (\circ)'},'FontWeight','bold')
        
    end
    xticks(-90:45:90)
    box off
    % title(sprintf('Lag %2.0f',name(lag)));
    axes1=gca
    set(axes1,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','off','YMinorTick','off'...
        ,'TickLength',[.015 .015],'Layer','bottom')
end
legend('T1','T1','T2','T2')
set(gcf,'Position',[   440   378   205   420])
print('Exp1TuningCurves.eps','-depsc')

%% Panel B
X= [ 1 2 3 5 7];
a=figure;
a.Color=[1 1 1];
for var =1:4
    subplot(2,2,var)
    
    T1= PsyStore(:,:,1,var);
    T2= PsyStore(:,:,2,var);
    
    [H,P,CI,STATS] = ttest(T1,T2);

   
    fprintf('Variable :%2.0f | L1:%2.4f |  L2:%2.4f |  L3:%2.4f | L5:%2.4f | L7:%2.4f\n',var,P)
    hold on
    
    
    T1mean=mean(T1,1);
    T1std=std(T1,1)/sqrt(size(PsyStore,1));
    
    T2mean=mean(T2,1);
    T2std=std(T2,1)/sqrt(size(PsyStore,1));
    
    
    errorbar(X,T1mean,T1std,'o-','Color',Cols(1,:),...
        'MarkerFaceColor',Cols(1,:),'LineWidth',2,'CapSize',0,...
        'MarkerSize',8)
    
    
    errorbar(X,T2mean,T2std,'o-','Color',Cols(2,:),...
        'MarkerFaceColor',Cols(2,:),'LineWidth',2','CapSize',0,...
        'MarkerSize',8)
    
    xlim([.5 7.5])
    if var ==3 || var ==4
        xlabel('Lag','FontWeight','bold')
    end
    if var ==3
        legend('T1','T2')
    end
    if var ==1
        point = [ .05 .05];
        
    elseif var ==2
        point = [ -25 25];
    elseif var==3
        point = [ 25 25];
        
        
    else
        point = [ .015 .015 ];
        
    end
    for i=1:length(P)
        
        if P(i)<.05/5
            plot([X(i) X(i)] ,[point],'k*')
        end
    end
    if var==1
     %   ylim([0 .24])
        yticks(0:.1:.3)
        ylim([.08 .2])
        ylabel('Amplitude (a.u.)','FontWeight','bold','FontAngle','italic')
    elseif var ==2
       ylim([-30 30])
        yticks(-20:20:20)
        
        ylabel('Center orientation (\circ)','FontWeight','bold','FontAngle','italic')
    elseif var ==3
        ylabel('Width (\circ)','FontWeight','bold','FontAngle','italic')
        yticks(10:10:50)
        ylim([10 40])
    elseif var==4
        ylabel('Baseline (a.u.)','FontWeight','bold','FontAngle','italic')
        yticks(0:.05:.1)
       ylim([0 .12])
    end
    xticks(X)
    axes1=gca;
    set(axes1,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','off','YMinorTick','off'...
        ,'TickLength',[.015 .015],'Layer','bottom')
    box off
end
print('Exp1TuningCurvesFits.eps','-depsc')


%% do the regression analysis

cd(dataFolder);

directory = dir('*.csv');
clear regressionResultStore %
for part = 1:length(directory)
    X2=csvread(directory(part).name);
    
      
    for lag =1:5
        clear Store
        
        if lag ==1
            name ='Lag 2';
        else
            name ='Lag 7';
            
        end
        for items =1:2
            trialC=1;
            
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
                    OrientationAverageT1= deg2rad(mod(X2(trial,X2(trial,5)+6),180))*2;
                    OrientationAverageT2= deg2rad(mod(X2(trial,X2(trial,6)+6),180))*2;
                    av_or=mod(rad2deg(circ_mean([OrientationAverageT1 OrientationAverageT2]')/2),180);
                    Store(trialC,:)=[T1report T2report  av_or Orientations ];
                    trialC=trialC+1;
                end
            end
            %%
            RegressionMatrix =exp(1i*(deg2rad(Store(:,3:end)))); % makes it a complex
            % number for easier regression
            if items ==1
                Response =exp(1i*(deg2rad(Store(:,1)))); %what orientation the participant
            elseif items ==2
                Response =exp(1i*(deg2rad(Store(:,2)))); %what orientation the participant
                
            else
                Response =exp(1i*(deg2rad(Store(:,3)))); %what orientation the participant
                
            end
            % reported
            Result=pinv(RegressionMatrix'*RegressionMatrix)*RegressionMatrix'*Response;
            % this does the regression
            
            regressionResultStore(part,lag,items,:)=Result;
  
        end
    end
end

%% Panel C
count=1;
a=figure;
a.Color =[1 1 1];
X=-5:7;
lags =[1 2 3 5 7];
Cols=linspecer(3);
for lag  =1:5
    
    subplot(5,1,count)
    
    for items=1:2
        Temp =squeeze(regressionResultStore(:,lag,items,:));
        tempMean =abs(mean(Temp,1));
        
        tempSTD = abs(std(Temp,1)/sqrt(size(Temp,1)));
        errorbar(X,tempMean,tempSTD,'o-','Color',Cols(items,:),...
            'MarkerFaceColor',Cols(items,:),'CapSize',0,'LineWidth',2)
        % plot(tempMean)
        hold on
    end
    count=count+1;
    %  ylim([0 1])
    yticks(0:.2:.4)
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

