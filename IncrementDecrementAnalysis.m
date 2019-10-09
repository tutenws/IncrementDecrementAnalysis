% This code was written to analyze increment-decrement experimental data
% collected on the UPENN AOSLO; you can select multiple data files and they
% will be combined prior to analysis, but take care that they share certain
% key parameters (subject, background, stimulus size, and stimulus
% spacings). There is some error handling to check for this. The code will
% generate plots of proportion seen for the various stimulus combinations.
% It also fits psychometric functions to the increments and decrements that
% are presented in isolation, and it uses these fits to generate
% predictions for independent-detector and summation-only models. Subject
% performance is compared to these models.
%
% NOTE #1: PF fitting here uses the Palamedes toolbox, so that will need to be
% installed and accessible on your Matlab path.
%
% NOTE #2: Over the evolution of this experiment, the controller repsonse that
% coresponding to "NO"/"NOT SEEN" may have changed. There is currently some
% checking around LINE 77 that goes on that reflects the extant spectrum of
% responses that have been used by the tested subjects, but if the plots
% don't make sense, check that the no response flag is coded correctly.
%
% 10-8-19       wst cleaned up existing code

%% Housekeeping
close all
clc

%% Select grouped data files
[fNames, pName] = uigetfile('*.mat', 'Select data files', 'MultiSelect', 'on');

% Check that you selected something
if isempty(fNames)
    error('No data files found; check folder');
end

% Conver single file selection to cell format
if ischar(fNames)
    fNames = {fNames};
end

% Pre-allocate prior to combining mat files
expData.testSeq = [];
expData.YesNoResponseMatrix = [];
expData.CFG.backgroundIntensity = [];
expData.CFG.subjectID = [];
expData.CFG.stimSize = [];
expData.CFG.stimCenterToCenterSpacings = [];

%% Load in the and compile data file(s)
for fileNum = 1:size(fNames,2)
    tempData = load(fullfile(pName, fNames{fileNum})); % Load data into a "temporary" structure
    expData.testSeq = [expData.testSeq; tempData.testSeq];
    expData.YesNoResponseMatrix = [expData.YesNoResponseMatrix; tempData.YesNoResponseMatrix];
    if fileNum == 1
        % These are key parameters which should be the same across all
        % analyzed files
        expData.CFG.backgroundIntensity = [expData.CFG.backgroundIntensity; tempData.CFG.backgroundIntensity];
        expData.CFG.subjectID = [expData.CFG.subjectID; tempData.CFG.subjectID];
        expData.CFG.stimSize = [expData.CFG.stimSize; tempData.CFG.stimSize];
        expData.CFG.stimCenterToCenterSpacings = [expData.CFG.stimCenterToCenterSpacings; tempData.CFG.stimCenterToCenterSpacings];
    else
        % Check that subsequent data files have the same subject ID, background intensity, stimulus size, and stimulus spacing 
        if expData.CFG.backgroundIntensity ~= tempData.CFG.backgroundIntensity | expData.CFG.subjectID ~= tempData.CFG.subjectID ...
                | expData.CFG.stimSize ~= tempData.CFG.stimSize | expData.CFG.stimCenterToCenterSpacings ~= tempData.CFG.stimCenterToCenterSpacings
            error('Data files contain different conditions. Check your selection.');
        end        
    end
end

%% Analyze the data

% First, find intensities for spot 1 and spot 2
spot1intensities = unique(expData.testSeq(:,1));
spot2intensities = unique(expData.testSeq(:,2));

% Pre-allocate output matrices
numSeenMatrix = nan(numel(spot2intensities), numel(spot1intensities),1);
numPresentedMatrix = numSeenMatrix;

% Response mapping (update manually)
if length(unique(expData.YesNoResponseMatrix)) == 2
    noResp = 1; % David data
else
    noResp = 4; % Early Will data
end

bump = 1;
for n = 1:numel(spot1intensities)
    for j = bump:numel(spot2intensities)
        stimRows = find(expData.testSeq(:,1)==spot1intensities(n) & expData.testSeq(:,2)==spot2intensities(j));
        if ~isempty(stimRows)
            
            % Add data to output matrices
            numPresentedMatrix(n,j) = length(stimRows);
            numSeenMatrix(n,j) = length(find(expData.YesNoResponseMatrix(stimRows)~=noResp));
            
            % Print to screen
            fprintf('Stim pair: [%.2f + %.2f]\n', spot1intensities(n), spot2intensities(j));
        end        
    end
    bump = bump+1;
end

% Plot the data
f1 = figure; hold on
set(f1, 'units', 'inches', 'position', [.25 .25 16 8]);
f1.Color = [1 1 1];

% Set the alpha shading so the matrix-style plot looks nicer (nan values
% won't map onto the color scheme; instead they will assume the color of
% the axis background)
imAlpha = ones(size(numSeenMatrix));
imAlpha(isnan(numSeenMatrix))=0;

%% Panel 1 -- Plot the raw yes/no detection data using imagesc

proportionSeenMatrixRaw = numSeenMatrix./numPresentedMatrix;
ax1 = subplot('Position', [0.025 0.525 0.875/4 0.875/4]); hold on
imagesc(proportionSeenMatrixRaw, 'AlphaData', imAlpha);

%-----------------------------PLOTTING-------------------------------------
% Axis stuff for plotting
axis square
ax1.Color = [1 1 1]; % Set background color to be white
ax1.XTick = 1:size(numSeenMatrix,2);
ax1.YTick = 1:size(numSeenMatrix,1);
ax1.TickLength = [0 0]; % Turn off tick marks
ax1.YAxisLocation = 'right';
ax1.XAxisLocation = 'top';
ax1.YDir = 'reverse';
xlabel('Stim 1 modulation', 'FontSize', 10)
ylabel('Stim 2 modulation', 'FontSize', 10, 'Rotation', 270, 'VerticalAlignment', 'bottom');
ax1.XRuler.Axle.LineStyle = 'none'; % Turn off axis lines without losing tick labels
ax1.YRuler.Axle.LineStyle = 'none'; % Turn off axis lines without losing tick labels
ax1.XLim = [0.5 size(numSeenMatrix,2)+0.5];
ax1.YLim = [0.5 size(numSeenMatrix,1)+0.5];
ax1.XTickLabel = cellstr(num2str(spot1intensities-expData.CFG.backgroundIntensity));
ax1.YTickLabel = cellstr(num2str(spot2intensities-expData.CFG.backgroundIntensity));
caxis([0 1]);

% Plot a box around the matrix elements
xValsBox = sort(repmat(0.5:1:size(numSeenMatrix,2)+.5, [1 2]));
xValsBox(end+1) = 0.5;
yValsBox = sort(repmat(0.5:1:size(numSeenMatrix,1)+.5, [1 2]));
yValsBox(1) = [];
yValsBox(end+1:end+2) = 0.5;
plot(xValsBox, yValsBox, 'k-', 'LineWidth', 1)

% Text label
title('(A) Raw data', 'FontSize', 12);

% Color bar stuff
c1 = colorbar('Location', 'westoutside');
c1.Label.String = 'Proportion seen';
c1.Label.FontSize = 10;
%-----------------------------PLOTTING-------------------------------------

%% Panel 2 -- Plot the detection data after correcting for false alarms

blankCol = find(spot1intensities==expData.CFG.backgroundIntensity);
blankRow = find(spot2intensities==expData.CFG.backgroundIntensity);
falseAlarmRate = proportionSeenMatrixRaw(blankRow, blankCol);

proportionSeenMatrixCorrectedForGuessing = (proportionSeenMatrixRaw-falseAlarmRate)./(1-falseAlarmRate);

ax2 = subplot('Position',[0.05+(0.875/4) 0.525 0.875/4 0.875/4]);
hold on
imagesc(proportionSeenMatrixCorrectedForGuessing, 'AlphaData', imAlpha);

%-----------------------------PLOTTING-------------------------------------
% Axis stuff for plotting
axis square
ax2.Color = [1 1 1]; % Set background color to be white
ax2.XTick = 1:size(numSeenMatrix,2);
ax2.YTick = 1:size(numSeenMatrix,1);
ax2.TickLength = [0 0]; % Turn off tick marks
ax2.YAxisLocation = 'right';
ax2.XAxisLocation = 'top';
ax2.YDir = 'reverse';
xlabel('Stim 1 modulation', 'FontSize', 10)
ylabel('Stim 2 modulation', 'FontSize', 10, 'Rotation', 270, 'VerticalAlignment', 'bottom');
ax2.XRuler.Axle.LineStyle = 'none'; % Turn off axis lines without losing tick labels
ax2.YRuler.Axle.LineStyle = 'none'; % Turn off axis lines without losing tick labels
ax2.XLim = [0.5 size(numSeenMatrix,2)+0.5];
ax2.YLim = [0.5 size(numSeenMatrix,1)+0.5];
ax2.XTickLabel = cellstr(num2str(spot1intensities-expData.CFG.backgroundIntensity));
ax2.YTickLabel = cellstr(num2str(spot2intensities-expData.CFG.backgroundIntensity));
caxis([0 1]);

% Plot a box around the matrix elements
plot(xValsBox, yValsBox, 'k-', 'LineWidth', 1)

% Put markers on the boxes that will be fit with the PF
plot(repmat(blankCol, [1 blankRow]), 1:1:blankRow, 'r:o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(blankCol:size(numSeenMatrix,2), repmat(blankRow, [1 1+size(numSeenMatrix,2)-blankCol]), 'k:+', 'LineWidth', 2, 'MarkerSize', 8);

% Text label
title('(B) Corrected for guessing', 'FontSize', 12);
%-----------------------------PLOTTING-------------------------------------

%% Fit the psychometric function
decMods = spot2intensities(spot2intensities<=expData.CFG.backgroundIntensity,:)-expData.CFG.backgroundIntensity;
decProp = proportionSeenMatrixCorrectedForGuessing(1:length(decMods), blankCol);
incMods = spot1intensities(spot1intensities>=expData.CFG.backgroundIntensity,:)-expData.CFG.backgroundIntensity;
incProp = proportionSeenMatrixCorrectedForGuessing(blankRow, blankCol:end)';

% Handle negative values of proportion seen that arise after correcting for
% guessing
decProp(decProp<0) = 0;
decProp(decProp>1) = 1;
incProp(incProp<0) =  0;
incProp(incProp>1) = 1;

% PF fit is on log scale, so set floor to be low but non-zero
decMods(decMods==0) = 0.001;
incMods(incMods==0) = 0.001;

% PF stuff
PF = @PAL_Logistic; % Choose your own adventure here
paramsFree = [1 1 0 0]; % Constraint flags for [threshold slope guess-rate lapse-rate]; 0 = fixed; 1 = free;

% Starting guess for the PF params
searchGridInc = [-1.5 4 0 0];
searchGridDec = searchGridInc;

% Do the fits
[paramsFitted_inc, LL_inc, exitflag_inc, ~] = PAL_PFML_Fit(log10(abs(incMods)), round(incProp.*1000), 1000*ones(size(incProp)), searchGridInc, paramsFree, PF);
[paramsFitted_dec, LL_dec, exitflag_dec, ~] = PAL_PFML_Fit(log10(abs(decMods)), round(decProp.*1000), 1000*ones(size(decProp)), searchGridDec, paramsFree, PF);

% Evaluate the fits
xEval = linspace(0,1,1000);
yEvalInc = PF(paramsFitted_inc, log10(xEval));
yEvalDec = PF(paramsFitted_dec, log10(xEval));

%% Panel 3 -- Two detector/probability summation plot

% Pre-allocate
twoDetectorCorrectedForGuessingMatrix = nan(size(numSeenMatrix));

bump = 1;
for n = 1:numel(spot1intensities)
    for j = bump:numel(spot2intensities)
        
        spot1Mod = spot1intensities(n)-expData.CFG.backgroundIntensity;
        spot2Mod = spot2intensities(j)-expData.CFG.backgroundIntensity;
        % Determine spot 1 detection rate
        if spot1Mod < 0 % Decrement
             spot1ProbCorrected = PF(paramsFitted_dec, log10(abs(spot1Mod)));
        elseif spot1Mod > 0 % Increment
             spot1ProbCorrected = PF(paramsFitted_inc, log10(abs(spot1Mod)));
        else
            spot1ProbCorrected = proportionSeenMatrixCorrectedForGuessing(blankRow, blankCol);
        end
        
        % Determine spot 2 detection rate
        if spot2Mod < 0 % Decrement
            spot2ProbCorrected = PF(paramsFitted_dec, log10(abs(spot2Mod)));
        elseif spot2Mod > 0 % Increment
            spot2ProbCorrected = PF(paramsFitted_inc, log10(abs(spot2Mod)));
        else
            spot2ProbCorrected = proportionSeenMatrixCorrectedForGuessing(blankRow, blankCol);
        end
        
        twoDetectorProbNotSeen = (1-spot1ProbCorrected).*(1-spot2ProbCorrected);
        twoDetectorCorrectedForGuessingMatrix(n,j) = 1-twoDetectorProbNotSeen; % probSeen is 1-probNotSeen
    end
    bump = bump+1;
end

ax3 = subplot('Position',[0.075+2*(0.875/4) 0.525 0.875/4 0.875/4]);
hold on
imagesc(twoDetectorCorrectedForGuessingMatrix, 'AlphaData', imAlpha);

%-----------------------------PLOTTING-------------------------------------
% Axis stuff for plotting
axis square
ax3.Color = [1 1 1]; % Set background color to be white
ax3.XTick = 1:size(numSeenMatrix,2);
ax3.YTick = 1:size(numSeenMatrix,1);
ax3.TickLength = [0 0]; % Turn off tick marks
ax3.YAxisLocation = 'right';
ax3.XAxisLocation = 'top';
ax3.YDir = 'reverse';
xlabel('Stim 1 modulation', 'FontSize', 10)
ylabel('Stim 2 modulation', 'FontSize', 10, 'Rotation', 270, 'VerticalAlignment', 'bottom');
ax3.XRuler.Axle.LineStyle = 'none'; % Turn off axis lines without losing tick labels
ax3.YRuler.Axle.LineStyle = 'none'; % Turn off axis lines without losing tick labels
ax3.XLim = [0.5 size(numSeenMatrix,2)+0.5];
ax3.YLim = [0.5 size(numSeenMatrix,1)+0.5];
ax3.XTickLabel = cellstr(num2str(spot1intensities-expData.CFG.backgroundIntensity));
ax3.YTickLabel = cellstr(num2str(spot2intensities-expData.CFG.backgroundIntensity));
caxis([0 1]);

% Plot a box around the matrix elements
plot(xValsBox, yValsBox, 'k-', 'LineWidth', 1)

% Text label
title('(C) Two-detector prediction', 'FontSize', 12);
%-----------------------------PLOTTING-------------------------------------

%% Panel 4, difference between data and predicition for two-detector

ax4 = subplot('Position',[0.1+3*(0.875/4) 0.525 0.875/4 0.875/4]);
hold on

% Subtract prediction from observed data
twoDetectorVsData = proportionSeenMatrixCorrectedForGuessing-twoDetectorCorrectedForGuessingMatrix;
imagesc(twoDetectorVsData, 'AlphaData', imAlpha);

%-----------------------------PLOTTING-------------------------------------
% Axis stuff for plotting
axis square
ax4.Color = [1 1 1]; % Set background color to be white
ax4.XTick = 1:size(numSeenMatrix,2);
ax4.YTick = 1:size(numSeenMatrix,1);
ax4.TickLength = [0 0]; % Turn off tick marks
ax4.YAxisLocation = 'right';
ax4.XAxisLocation = 'top';
ax4.YDir = 'reverse';
xlabel('Stim 1 modulation', 'FontSize', 10)
ylabel('Stim 2 modulation', 'FontSize', 10, 'Rotation', 270, 'VerticalAlignment', 'bottom');
ax4.XRuler.Axle.LineStyle = 'none'; % Turn off axis lines without losing tick labels
ax4.YRuler.Axle.LineStyle = 'none'; % Turn off axis lines without losing tick labels
ax4.XLim = [0.5 size(numSeenMatrix,2)+0.5];
ax4.YLim = [0.5 size(numSeenMatrix,1)+0.5];
ax4.XTickLabel = cellstr(num2str(spot1intensities-expData.CFG.backgroundIntensity));
ax4.YTickLabel = cellstr(num2str(spot2intensities-expData.CFG.backgroundIntensity));
caxis([-1 1]);

% Plot a box around the matrix elements
plot(xValsBox, yValsBox, 'k-', 'LineWidth', 1)

% Text label
title('(D) Difference (B - C)', 'FontSize', 12);

% Color bar stuff
c2 = colorbar('Location', 'eastoutside');
c2.Label.String = 'Difference';
c2.Label.FontSize = 10;
%-----------------------------PLOTTING-------------------------------------

%% Panel 6, frequency of seeing plots + PF

ax6 = subplot('Position',[0.05+(0.875/4) 0.1 0.875/4 0.875/4]);
hold on
% Plot the PF
plot(abs(decMods), decProp, 'ro', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(xEval, yEvalDec, 'r:', 'LineWidth', 2);
plot(abs(incMods), incProp, 'k+', 'LineWidth', 2, 'MarkerSize', 10)
plot(xEval, yEvalInc, 'k:', 'LineWidth', 2);
hLeg = legend('Decrements', '   fit',  'Increments', '   fit', 'Location', 'SouthEastOutside');
hLeg.Box = 'off';

%-----------------------------PLOTTING-------------------------------------
% Set the axis to be square
axis square
box off
ax6.Color = [1 1 1]; % Set background color to be white
ax6.XColor = [0 0 0];
ax6.YColor = [0 0 0];
ax6.LineWidth = 1.25;
ax6.XTick = 0:0.2:1;
ax6.YTick = 0:0.2:1;
ax6.TickLength = [0.025 0.025];
ax6.TickDir = 'in';
xlabel('Absolute stimulus modulation', 'FontSize', 12)
ylabel('Proportion seen', 'FontSize', 12);
ax6.XLim = [0 1];
ax6.YLim = [0 1];

% Text label
title('(E) Frequency of seeing', 'FontSize', 12);
%-----------------------------PLOTTING-------------------------------------

%% Panel 7 -- Summation-only plot, derived from PF fits

summationOnlyCorrectedForGuessingMatrix = nan(size(numSeenMatrix));
bump = 1;
for n = 1:numel(spot1intensities)
    for j = bump:numel(spot2intensities)
        
        % Determine modulations
        spot1Mod = spot1intensities(n)-expData.CFG.backgroundIntensity;
        spot2Mod = spot2intensities(j)-expData.CFG.backgroundIntensity;
        % Add them together
        summedMod = spot1Mod + spot2Mod;
        
        if summedMod < 0 % Use decrement PF
            % Feed summed mod into the PF fit; note that the sum may have
            % not actually been tested as an isolated stimulus -- that's
            % why we need to do the PF fit
            propSeen =  PF(paramsFitted_dec, log10(abs(summedMod)));
        elseif summedMod > 0 % Use increment PF
            propSeen =  PF(paramsFitted_inc, log10(abs(summedMod)));
        else
            propSeen = 0;
        end        
        summationOnlyCorrectedForGuessingMatrix(n,j) = propSeen;
    end
    bump = bump+1;
end

%-----------------------------PLOTTING-------------------------------------
ax7 = subplot('Position',[0.075+2*(0.875/4) 0.10 0.875/4 0.875/4]);
hold on
imagesc(summationOnlyCorrectedForGuessingMatrix, 'AlphaData', imAlpha);

% Axis stuff for plotting
axis square
ax7.Color = [1 1 1]; % Set background color to be white
ax7.XTick = 1:size(numSeenMatrix,2);
ax7.YTick = 1:size(numSeenMatrix,1);
ax7.TickLength = [0 0]; % Turn off tick marks
ax7.YAxisLocation = 'right';
ax7.XAxisLocation = 'top';
ax7.YDir = 'reverse';
xlabel('Stim 1 modulation', 'FontSize', 10)
ylabel('Stim 2 modulation', 'FontSize', 10, 'Rotation', 270, 'VerticalAlignment', 'bottom');
ax7.XRuler.Axle.LineStyle = 'none'; % Turn off axis lines without losing tick labels
ax7.YRuler.Axle.LineStyle = 'none'; % Turn off axis lines without losing tick labels
ax7.XLim = [0.5 size(numSeenMatrix,2)+0.5];
ax7.YLim = [0.5 size(numSeenMatrix,1)+0.5];
ax7.XTickLabel = cellstr(num2str(spot1intensities-expData.CFG.backgroundIntensity));
ax7.YTickLabel = cellstr(num2str(spot2intensities-expData.CFG.backgroundIntensity));

% Color axis limits
caxis([0 1]);

% Plot a box around the matrix elements
plot(xValsBox, yValsBox, 'k-', 'LineWidth', 1)

% Text label
title('(F) Summation-only prediction', 'FontSize', 12);
%-----------------------------PLOTTING-------------------------------------

%% Panel 8, difference between data and predicition for summation-only
ax8 = subplot('Position',[0.1+3*(0.875/4) 0.1 0.875/4 0.875/4]);
hold on
summationVsData = proportionSeenMatrixCorrectedForGuessing-summationOnlyCorrectedForGuessingMatrix;

%-----------------------------PLOTTING-------------------------------------
imagesc(summationVsData, 'AlphaData', imAlpha);

% Axis stuff for plotting
axis square
ax8.Color = [1 1 1]; % Set background color to be white
ax8.XTick = 1:size(numSeenMatrix,2);
ax8.YTick = 1:size(numSeenMatrix,1);
ax8.TickLength = [0 0]; % Turn off tick marks
ax8.YAxisLocation = 'right';
ax8.XAxisLocation = 'top';
ax8.YDir = 'reverse';
xlabel('Stim 1 modulation', 'FontSize', 10)
ylabel('Stim 2 modulation', 'FontSize', 10, 'Rotation', 270, 'VerticalAlignment', 'bottom');
ax8.XRuler.Axle.LineStyle = 'none'; % Turn off axis lines without losing tick labels
ax8.YRuler.Axle.LineStyle = 'none'; % Turn off axis lines without losing tick labels
ax8.XLim = [0.5 size(numSeenMatrix,2)+0.5];
ax8.YLim = [0.5 size(numSeenMatrix,1)+0.5];
ax8.XTickLabel = cellstr(num2str(spot1intensities-expData.CFG.backgroundIntensity));
ax8.YTickLabel = cellstr(num2str(spot2intensities-expData.CFG.backgroundIntensity));
caxis([-1 1]);

% Plot a box around the matrix elements
plot(xValsBox, yValsBox, 'k-', 'LineWidth', 1)

% Text label
title('(G) Difference (B - F)', 'FontSize', 12);
% Color bar stuff
c2 = colorbar('Location', 'eastoutside');
c2.Label.String = 'Difference';
c2.Label.FontSize = 10;
%-----------------------------PLOTTING-------------------------------------

%% Plot stars on best predictions for each square in the difference maps
predictionMatrix = nan(size(numSeenMatrix,1), size(numSeenMatrix,2), 2);
predictionMatrix(:,:,1) = twoDetectorVsData;
predictionMatrix(:,:,2) = summationVsData;

[c,i] = min(abs(predictionMatrix), [],3);

i(isnan(c)) = nan;
i(:,blankCol)= nan;
i(blankRow,:) = nan;

[twoRow, twoCol] = find(i==1);
[sumRow, sumCol] = find(i==2);

axes(ax4); % Switch to panel 4
hold on
plot(twoCol, twoRow, 'w*');

axes(ax8); % Switch to panel 8
hold on
plot(sumCol, sumRow, 'w*');


%% Write some text into the spare panel
ax5 = subplot('Position',[0.025 0.1 0.875/5 0.875/4]);

ax5.Color = [1 1 1]; % Set background color to be white
ax5.XColor = [1 1 1]; % Set background colors to be white
ax5.YColor = [1 1 1]; % Set background colors to be white

ax5.XLim = [0 1];
ax5.YLim = [0 1];

text(0.9, 0.9, ['Subject: ' expData.CFG.subjectID], 'HorizontalAlignment', 'right', 'FontSize', 14)
text(0.9, 0.7, 'Combined folders', 'HorizontalAlignment', 'right', 'FontSize', 14, 'Interpreter', 'none')
text(0.9, 0.5, ['Stim size (pixels): ' num2str(expData.CFG.stimSize)], 'HorizontalAlignment', 'right', 'FontSize', 14)
text(0.9, 0.3, ['Stim spacing factor: ' num2str(expData.CFG.stimCenterToCenterSpacings)], 'HorizontalAlignment', 'right', 'FontSize', 14)
text(0.9, 0.1, ['Background level (0-1): ' num2str(expData.CFG.backgroundIntensity)], 'HorizontalAlignment', 'right', 'FontSize', 14)

%% New figure, showing data and predictions along the negative diagonal, per DHB suggestion on 10/8/19

f2 = figure;
set(f2, 'Units', 'inches', 'Position', [2 2 7.5 5], 'Color', [1 1 1]);
hold on
% Plot data
plot(diag(proportionSeenMatrixCorrectedForGuessing), 'ko-', 'MarkerFaceColor', 'w', 'LineWidth', 2, 'MarkerSize', 6);
% Plot independent detector prediction
plot(diag(twoDetectorCorrectedForGuessingMatrix), 'rs-', 'MarkerFaceColor', 'w', 'LineWidth', 2, 'MarkerSize', 6);
% Plot summation-only prediction
plot(diag(summationOnlyCorrectedForGuessingMatrix), '+:', 'LineWidth', 2, 'MarkerSize', 10);

legend('Data', 'Independent detector', 'Summation-only', 'Location', 'North')

set(gca, 'LineWidth', 1.25, 'TickDir', 'out')
ylabel('Proportion seen', 'FontSize', 14);

% X axis label
modulationsString = (num2str(spot2intensities-expData.CFG.backgroundIntensity));
xLabels = cellstr([modulationsString repmat(', ', [length(modulationsString),1]) modulationsString]);
set(gca, 'XTickLabel', xLabels, 'XTickLabelRotation', 45, 'TickLength', [0.02 0.02], 'XColor', [0 0 0], 'YColor', [0 0 0]);
xlabel('Stimulus pair', 'FontSize', 14);



%% Save the figures

saveDir = [pwd '\IncrementDecrementPlots\'];
if ~isdir(saveDir)
    mkdir(saveDir);
end

if length(fNames)>1
    folderString = '_combined_';
else
    fStr = fNames{1};
    [strInd] = strfind(fStr, '_');
    if length(strInd) > 2
        folderString = fStr(strInd(1):strInd(end));
    else
        
    end
end
saveNameMultiPanel = fullfile(saveDir, [expData.CFG.subjectID '_spacing_' num2str(expData.CFG.stimCenterToCenterSpacings) folderString 'detectionAnalysis.png']);
print(f1, saveNameMultiPanel, '-dpng');

saveNameDiagonal = fullfile(saveDir, [expData.CFG.subjectID '_spacing_' num2str(expData.CFG.stimCenterToCenterSpacings) folderString 'modelEvaluation_likeType.png']);
print(f2, saveNameDiagonal, '-dpng');



