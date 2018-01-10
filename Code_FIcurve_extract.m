%{
To extract F-I curve data from .mat file (from spike2 file)
and save the results into ['DATA_',filename,'.mat']

Note
1) data for batch anlaysis should be collected from same protocol (same step amplitude,sampling rate, etc)

%}
%%
clear all;clc

Fidx{1}=dir('Quan_*.mat');

%%
for kk=1:length(Fidx)
    for k=1:length(Fidx{kk})
        %% ------------- step#1: To load data from .mat file ---------------------------
        filename=Fidx{kk}(k).name;
        load(filename)
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',num2str(kk),'/',num2str(length(Fidx)),'--',filename])
        %% ------------- step#2: To read data from structure ---------------------------
        clear(('*Ch31*'));
        
        varia=who('*Ch*');
        
        for i=1:length(varia)
            if isfield(eval(varia{i}),'units')
                str=eval([varia{i},'.units']);
                
                if ismember('pA',str)
                    Im=eval([varia{i} '.values']);
                    Im=Im-mean(Im(1:10));
                    dtI=eval([varia{i} '.interval']);
                end
                
                if ismember('mV',str)
                    Vm=eval([varia{i} '.values']);
                    dtV=eval([varia{i} '.interval']);
                end
                
            end
        end
        
        clear (varia{:})
        
        %% ------------- step#3:  To find stimulus time point ---------------------------
        TimeStim = eventedgefinder(Im,2,1/dtI,0.2,0.5,1,0).*dtI;
        StimDura=range(TimeStim,2);
        TimeStim=TimeStim(:,1);
        
        %% ------------- step#4:  To extract waveform  ---------------------------
        pretime=0.1;
        postime=0.8;
        
        waveformV=[];
        waveformI=[];
        Vrest=[];
        for i=1:length(TimeStim)
            dataidxV=round(TimeStim(i)/dtV-pretime/dtV:TimeStim(i)/dtV+postime/dtV);
            baseidxV=round(TimeStim(i)/dtV-pretime/dtV:TimeStim(i)/dtV);
            dataidxI=round(TimeStim(i)/dtI-pretime/dtI:TimeStim(i)/dtI+postime/dtI);
            baseidxI=round(TimeStim(i)/dtI-pretime/dtI:TimeStim(i)/dtI);
            
            Vrest=[Vrest;mean(Vm(baseidxV))];
            waveformV=[waveformV,Vm(dataidxV)];
            waveformI=[waveformI,Im(dataidxI)-mean(Im(baseidxI))];
        end
        tspanV=linspace(-pretime,postime,size(waveformV,1));
        tspanI=linspace(-pretime,postime,size(waveformI,1));
        
        %% ------------- step#4:  To make some necessary Tags   ---------------------------
        % Tags #1: firing rate
        spkFreq=[];
        
        for i=1:size(waveformV,2)
            peaktime=peakfinder(waveformV(:,i),20,Vrest(i)+10, 1, 0).*dtV;
            
            % to exclude wrong detection
            worngidx=[];
            for j=1:length(peaktime)
                idx=round(peaktime(j)/dtV-0.01/dtV:peaktime(j)/dtV);
                temp=min(waveformV(idx,i));
                if temp>waveformV(idx(end),i)-10;
                    worngidx=[worngidx;j];
                end
            end
            peaktime(worngidx)=[];
            
            % to exclude spontaneous  firing
            if ~isempty(peaktime);peaktime(peaktime>StimDura(i)+0.01+pretime|peaktime<pretime-0.01)=[];end
            
            spkFreq=[spkFreq;length(peaktime)./StimDura(i)];
            
            % to  make sure whether we have find the right SPIKE
            figure(1),clf
            plot(tspanV,waveformV(:,i)),hold on
            plot(tspanV(round(peaktime/dtV)),waveformV(round(peaktime/dtV),i),'ro')
            drawnow
            
        end
        
        % Tags #2: stimulus amplitude
        % (modified it according to your specific protocol)
        StimAmp=mean(waveformI(tspanI>0.1&tspanI<0.4,:))-mean(waveformI(tspanI<0,:));
        %         StimAmp(StimAmp<200)=round(StimAmp(StimAmp<200)/50)*50;
        %         StimAmp(StimAmp>=200)=round(StimAmp(StimAmp>=200)/100)*100;
        StimAmp(StimAmp<200)=round(StimAmp(StimAmp<200)/10)*10;
        StimAmp(StimAmp>=200)=round(StimAmp(StimAmp>=200)/100)*100;
        
        save(['DATA_',filename],'tspanV','tspanI','waveformV','waveformI','TimeStim','Vrest','spkFreq','StimAmp')
    end
    
end


% (c) 2018 He quansheng
% State Key Laboratory of Cognitive Neuroscience and Learning & IDG/McGovern 
% Institute for Brain Research, Beijing Normal University
% quanshe@mail.bnu.edu.cn