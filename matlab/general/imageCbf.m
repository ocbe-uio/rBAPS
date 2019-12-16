function imageCbf(action)
%Tahan funktioon ohjataan image_figuren callback:it.

switch action

case 'save_image'
    saveImage;
    
case 'export_bmp'
    export('bmp');
    
case 'export_jpg'
    export('jpg');
    
end

function saveImage
%Saves information needed to reconstruct the image later.

[filename,pathname] = uiputfile('*.mat','Save Figure');
if (filename == 0) & (pathname == 0)
    %Cancel was clicked.
    return;
end;
image_file_name = [pathname filename];
tiedot = get(gcbf,'UserData');
% tiedot on tietue, joka sis‰lt‰‰ info:n ja popnames:in

% save(image_file_name,'tiedot');
save(image_file_name,'tiedot','-v7.3'); % added by Lu Cheng, 08.06.2012


%-------------------------------------------------------------------


function export(format)
%Saves a figure in a format which has been given
%as a parameter. Exported images cannot be opened using BAPS.

[filename,pathname] = uiputfile(['*.' format], ['Export to ' format]);
if (filename == 0) & (pathname == 0)
    %Cancel was pressed:
    return;
end;
filename = checkTheFormat(format,filename);
resultfilename = [pathname filename];
print(resultfilename);

%---------------------------------------------------------

function newfilename = checkTheFormat(format,oldfilename)
%Checks if the 'oldfilename' has ending *.'format'. If not, ending
%will be added to newfilename.

if length(oldfilename) < 4
    newfilename = [oldfilename '.' format];
elseif isequal(oldfilename(end-3: end), ['.' format])
    newfilename = oldfilename;
elseif any(oldfilename == '.')
    n = 1;
    while ~isequal(oldfilename(n),'.')
        n = n+1;
    end;
    newfilename = [oldfilename(1:n) format];
else
    newfilename = [oldfilename '.' format];
end;