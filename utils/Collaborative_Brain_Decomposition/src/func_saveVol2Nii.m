function func_saveVol2Nii(resFile,resFileName,maskName,outDir,saveName)


res = load(resFileName);
initV = res.(resFile)(:);

maskNii = load_untouch_nii(maskName);

if ~exist(outDir,'dir')
    mkdir(outDir);
end

% smInd = initV ./ max(repmat(max(initV),size(initV,1),1),eps) < 1e-2;
% initV(smInd) = 0;

K = size(initV,2);
for ki=1:K
    if ki<10
        kStr = ['00',num2str(ki)];
    elseif ki<100
        kStr = ['0',num2str(ki)];
    else
        kStr = num2str(ki);
    end
    disp(['save nii -- icn ',saveName]);
    
    kNii = maskNii;
    kNii.img(maskNii.img~=0) = initV(:,ki);
    
%     outName = [outDir,filesep,'icn_',kStr,'.nii.gz'];
    outName = [outDir,filesep,saveName,'.nii'];
    save_untouch_nii(kNii,outName);
end


end
