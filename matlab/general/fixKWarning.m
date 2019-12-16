%-------------------------------------------------------------------------
function fixedK = fixKWarning()
% The function is to ask further confirmation of running fixed K clustering
button = questdlg('You have selected a clustering analysis with the Fixed-K Mode. Are you sure to continue?',...
    'Fixed-K Mode');
switch button
    case 'Yes'
        fixedK = 1;
    case 'No'
        fixedK = 0;
        goToFixedK
    case 'Cancel'
        fixedK = 0.5;
end