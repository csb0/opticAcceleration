function [] = calcCurlKurtosis(eyeStruct)
%CALCCURLKURTOSIS calculates kurtosis of retinal curl for each recorded
%frame

varNames = fieldnames(eyeStruct);
for i=1:length(varNames)
    eval([varNames{i} '=eyeStruct.' varNames{i} ';']);
end

varNames = fieldnames(simDeets);
for i=1:length(varNames)
    eval([varNames{i} '=simDeets.' varNames{i} ';']);
end
%%

for fr = eyeStruct.startFrame+1: eyeStruct.endFrame
    
    if eyeStruct.headFlowBool
        thisEyeD = 'head';
    else
        if twoEyeBool
            thisEyeD = 'lEye';
        else
            thisEyeD = 'eye';
        end
    end
    
    curlClim = .01;
    
    thisCurl = eyeStruct.(thisEyeD).curlFr(fr,:,:);
    
    pd = fitdist(thisCurl(:),'Normal');
    
    mu(fr) = pd.mu;
    sig(fr) = pd.sigma;
    
    
    x = linspace(-curlClim, curlClim, 100);
    
    kurt(fr) = kurtosis(thisCurl(:),0);
    
    figure(324873);clf
    subplot(141)
    histogram(thisCurl,x,'Normalization','pdf')
    hold on
    plot(x, pdf(pd,x))
%     ylim([0 .1])
    xlim([-curlClim, curlClim])
    
    subplot(142)
    plot(kurt,'.-')
    hold on
    xlim([eyeStruct.startFrame eyeStruct.endFrame])
    
    subplot(143)
    plot(mu,'o-')
    xlim([eyeStruct.startFrame eyeStruct.endFrame])
    
    subplot(144)
    plot(sig,'p-')
    xlim([eyeStruct.startFrame eyeStruct.endFrame])
    drawnow
end


