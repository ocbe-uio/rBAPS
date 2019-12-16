%-------------------------------------------------------------------------

function goToFixedK
% The function to do fix K algorithm

h0 = findobj('Tag','fixk_menu');
h1 = findobj('Tag','mix_text');
old = get(h0, 'UserData');
if old == 0
    set(h0,'UserData',1, 'label', 'Disable Fixed-K Clustering');
    set(h1,'String', 'Population mixture analysis (Fixed-K Mode)');
    % disable all the non-relevant buttons
    set(findobj('Tag','partitioncompare_menu'),'Enable','off');
%     set(findobj('Tag','file_menu'),'Enable','off');
%     set(findobj('Tag','admix_text'),'Enable','off');
%     set(findobj('Tag','admix1_button'),'Enable','off');
%     set(findobj('Tag','admix2_button'),'Enable','off');
    disp('Fixed-K Mode is enabled.');
else
    set(h0,'UserData',0, 'label', 'Enable Fixed-K Clustering');
    set(h1,'String', 'Population mixture analysis');
    set(findobj('Tag','partitioncompare_menu'),'Enable','on');
%     set(findobj('Tag','file_menu'),'Enable','on');
%     set(findobj('Tag','admix_text'),'Enable','on');
%     set(findobj('Tag','admix1_button'),'Enable','on');
%     set(findobj('Tag','admix2_button'),'Enable','on');
    disp('Fixed-K Mode is disabled.');
end
