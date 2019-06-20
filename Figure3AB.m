clearvars;
SetUp
dataFolder= [dataFolder '/Exp1'];
cd(dataFolder);
directory = dir('*.csv');
clear Store
Store=nan(length(directory),120,5,4);
for part = 1:length(directory)
    disp(part)
    X2=csvread(directory(part).name);
    
    for lag =1:5
        trialC=1;
        for trial =1:length(X2)
            
            if X2(trial,4)==lag
                T1pos= X2(trial,5);
                T1actual=mod(X2(trial,6+T1pos),180);
                T1report =mod(X2(trial,27),180);
                
                
                T2pos= X2(trial,6);
                T2actual=mod(X2(trial,6+T2pos),180);
                T2report =mod(X2(trial,28),180);
                
                Store(part,trialC,lag,:)=[ T1actual T1report T2actual T2report ];
                trialC=trialC+1;
            end
        end
    end
end
%%
%% Panel B
close all
lb = [-50 5  -5 -5];
ub = [+50 180 5 5];
x0=[ -1 60 0 0];
f=@(c,xdata)c(3)+c(1).*(1/c(2)).*(xdata-c(4)).*exp(-(((xdata-c(4)).^2)/(2*c(2).^2)));
options = optimoptions('lsqcurvefit','Display','off');

clear counts StoreE
for target=1:2
    if target ==1
        t=[1 2];
        
    else
        t=[3 4];
    end
    for lag =1:5
        subplot(1,5,lag)
        dat= squeeze(Store(:,:,lag,:));
        actual = dat(:,:,t(1));
        response=dat(:,:,t(2));
        actual=mod(actual,180);
        if target==1
            differenceT1T2=dat(:,:,1)-dat(:,:,3);
        else
            differenceT1T2=dat(:,:,3)-dat(:,:,1);
            
        end
        differenceT1T2=mod(differenceT1T2+90,180)-90;
        
        
        for p=1:length(directory)
            
            dResponse = response(p,:);
            dActual = actual(p,:);
            differenceT1T2a=differenceT1T2(p,:);
            error=dResponse-dActual;
            error=mod(error+90,180)-90;
            
            
            psy =  lsqcurvefit(f,x0,differenceT1T2a,...
                error,lb,ub,options);
            psyStoreTemp(p,lag,target,:)=psy;
            
            
            hold on
        end
    end
    
end

%% Panel B

a=figure;
a.Color=[1 1 1];
n=-1:5;

clear ax
Cols=linspecer(2);
realLag=[1 2 3 5 7];
hold on
for t=1:2
    datReal= squeeze(psyStoreTemp(:,:,t,1));
    
    
    [H,P,CI,STATS] = ttest(datReal);
    
    datM=collapse(datReal,1);
    datSE= std(datReal,1)/sqrt(size(datReal,1));
    errorbar(realLag,datM,datSE,'-o','Color',Cols(t,:),...
        'MarkerFaceColor',Cols(t,:),'LineWidth',2,'CapSize',0)
    
    for i=1:5
        if P(i) <.05/5
            plot(realLag(i),-0+t*5,'*','Color',Cols(t,:),'MarkerSize',8)
        end
        
    end
    plot([0 8],[0 0],'k')
    xticks(realLag)
    box off
    xlabel('Lag')
    ylabel('Bias magnitude')
    yticks([-40:20:0])
    xlim([0.5 7.5])
end
set(gcf,'Position',[1000        1136         212         202])

axes1=gca;
set(axes1,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
    'YAxisLocation','left','LineWidth',1.5,...
    'XMinorTick','off','YMinorTick','off'...
    ,'TickLength',[.015 .015],'Layer','bottom')

set(gcf,    'Renderer','painters');


%% Panel A
            bins= -90:20:90;

Cols = linspecer(2);
clear counts StoreE
for target=1:2
    if target ==1
        t=[1 2];
        
    else
        t=[3 4];
    end
    for lag =1:5
        subplot(1,5,lag)
        dat= squeeze(Store(:,:,lag,:));
        actual = dat(:,:,t(1));
        response=dat(:,:,t(2));
        actual=mod(actual,180);
        if target==1
            differenceT1T2=dat(:,:,1)-dat(:,:,3);
        else
            differenceT1T2=dat(:,:,3)-dat(:,:,1);
            
        end
        differenceT1T2=mod(differenceT1T2+90,180)-90;
        
        
        for p=1:length(directory)
            
            dResponse = response(p,:);
            dActual = actual(p,:);
            differenceT1T2a=differenceT1T2(p,:);
            error=dResponse-dActual;
            error=mod(error+90,180)-90;
            for i=1:length(bins)-1
                ind = differenceT1T2a>bins(i) & differenceT1T2a<bins(i+1);
                StoreE(p,lag,target,i,:)=median(error(ind));
                
            end
            
        end
    end
    
end


%% fit the non-permuted data
a=figure;
a.Color=[1 1 1];
X=-90:90;
lb = [-50 5  -10];
ub = [+50 180 10];
f=@(c,xdata)c(3)+c(1).*(1/c(2)).*(xdata-0).*exp(-(((xdata-0).^2)/(2*c(2).^2)));

plot(X,f([10 45 10],X))


x0=[-10 30 0];
options = optimoptions('lsqcurvefit','Display','off');

%
realLag=[1 2 3 5 7];

for l=1:5
    for t=1:2
        subplot(1,5,l)
        dat = squeeze(StoreE(:,l,t,:));
        
        datM= nanmean(dat);
        mean(datM)
        datSE =std(dat)/sqrt(size(StoreE,1));
        
        psy =  lsqcurvefit(f,x0,bins(1:end-1),...
            datM,lb,ub,options);
        
        psy
        errorbar(bins(1:end-1),datM,datSE,'o','Color',Cols(t,:),...
            'MarkerFaceColor',Cols(t,:),'LineWidth',2,'CapSize',0)
        
        fitLine = f(psy,X);
        
        hold on
        
        plot(X,fitLine,'Color',Cols(t,:),'LineWidth',2)
        ylim([-30 30])
        hold on
        xticks([-90:90:180])
        
    end
    yticks([-30:10:30])
    title(['Lag ',num2str(realLag(l))])
    box off
    if l==1
        
        ylabel('Orientation error (\circ)')
    elseif l==3
        xlabel('T1 orientation relative to T2 (\circ)')
        set(gca,'Yticklabel',[])
        
        
    else
        set(gca,'Yticklabel',[])
        
    end
    if l==5
               legend('T1','','T2','')
 
    end
end
set(gcf,    'Renderer','painters');
set(gcf,'Position',[   560   766   560   182])
