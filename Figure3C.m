clearvars;
SetUp
dataFolder= [dataFolder '/Exp1'];
cd(dataFolder);
directory = dir('*.csv');

%%
clear Store
Store=nan(length(directory),120,5,2,7);
for part = 1:length(directory)
    disp(part)
    X2=csvread(directory(part).name);
    
    
    for lag =1:5
        trialC=1;
        for trial =1:length(X2)
            
            if X2(trial,4)==lag
                for items =1:2
                    if items ==1
                        ItemPosition= X2(trial,5);
                    else
                        ItemPosition= X2(trial,6);
                        
                    end
                    
                    T1pos= ItemPosition;
                    T1actual=mod(X2(trial,6+T1pos),180);
                    
                    if items ==1
                        T1report =mod(X2(trial,27),180);
                    else
                        T1report =mod(X2(trial,28),180);
                        
                    end
                    
                    ItemPosition=ItemPosition+6;
                    Orientations=X2(trial,ItemPosition-1:ItemPosition+4);
                    
                    Store(part,trialC,lag,items,:)=[ T1actual-T1report Orientations-T1actual ];
                end
                trialC=trialC+1;
            end
        end
    end
end
%%
close all

clear counts psyStoreTemp
X=-90:90;
lb = [-50 5  0 0];
ub = [+50 180 0 0];
x0=[ 0 60 0 0];
f=@(c,xdata)c(3)+c(1).*(1/c(2)).*(xdata-c(4)).*exp(-(((xdata-c(4)).^2)/(2*c(2).^2)));
options = optimoptions('lsqcurvefit','Display','off');

psyStoreTemp=zeros(size(Store,1),5,2,size(Store,5),4);
for target=1:2
   
    for lag =1:5
        subplot(1,5,lag)
        dat= squeeze(Store(:,:,lag,target,:));
   
                  error=mod(dat(:,:,1)+90,180)-90;


        parfor p=1:size(Store,1)
            for item =1:6
                difference=dat(p,:,1+item)
                difference = mod(difference+90,180)-90;
                psy =  lsqcurvefit(f,x0,difference,...
                    error(p,:),lb,ub,options);   
                psyStoreTemp(p,lag,target,item,:)=psy;
                
               %
            end

        end
    end
end

%%
close all
a=figure;
a.Color=[1 1 1];
n=-1:5;
    
clear ax
Cols=linspecer(2);
realLag=[1 2 3 5 7];
for lag =1:5
    ax(lag)=subplot(1,5,lag);
    hold on
    for t=1:2
        datReal= squeeze(psyStoreTemp(:,lag,t,:,1));
        

        [H,P,CI,STATS] = ttest(datReal);
        
        datM=collapse(datReal,1);
        datSE= std(datReal,1)/sqrt(size(datReal,1));
        errorbar(n,datM',datSE,'-o','Color',Cols(t,:),...
            'MarkerFaceColor',Cols(t,:),'LineWidth',2,'CapSize',0)
        
    
    for i=1:5
        if P(i) <.05/5
            
            plot(n(i),-40+t*5,'*','Color',Cols(t,:),'MarkerSize',8)
        end
        
    end
    end
    plot([realLag(lag) realLag(lag) ],[-30 30],'--','Color',Cols(1,:))
        plot([-realLag(lag) -realLag(lag) ],[-30 30],'--','Color',Cols(2,:))

    plot([0 0 ],[-30 30],'-.k')
    
    title(['Lag ',num2str(realLag(lag))])
    
    plot([- 4 9],[0 0 ],'k')
    xticks(-1:1:10)
    box off
    yticks([-40:20:40])
    if lag==1
        xlabel('Item posistion releative to target')
                ylabel('Bias magnitude')

    else
        
       set(gca,'Yticklabel',[])
        
    end
    
    
    ylim([ -40 10])
    xlim([-1.5 4])
    
    axes1=gca;
set(axes1,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
    'YAxisLocation','left','LineWidth',1.5,...
    'XMinorTick','off','YMinorTick','off'...
    ,'TickLength',[.015 .015],'Layer','bottom')
end
cd(figFile)


set(gcf,'Position',[   560   766   560   182])
print('BiasDistractors.eps','-depsc')
%%
 