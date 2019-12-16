function out = CSEFlagDialog(items, varargin)
%CSEFLAGDIALOG makes a GUI dialog to select from options of diverse types 
% (checkboxes, radiobuttons, text-inputs and popup lists). Some options can 
% be linked, i.e. be mutually exclusive or be only enabled according to the 
% value of another option.
% It is an extension of CSEOptionDialog.
%
% out = CSEFlagDialog(items)
% out = CSEFlagDialog(items, title)
% out = CSEFlagDialog(items, title, msg)
% out = CSEFlagDialog(items, title, msg, options)
%
% items is an array structure containing the options, in the following fields:
% .name      is a string (the name of the option)
% .values    is a cell array of possible answer.
%             if empty, it is considered to contain "Yes" and "No", i.e.
%               the option is a checkbox or a radio-button.
%             if it contains only one numeric element, it is considered a
%               header, and not an active control.
%             if it contains only one text element, it is considered to be a
%               text input field, with the .values as default value.
%             Otherwise, it is a popup list with the choices given in values.
% .linked    [optional field] is an array of index in items indicating which
%            options are linked with the option. Linked options
%            will be grayed out when the option is set to "No".
%            If .linked contains negative elements, those will be grayed
%            out when the option is set to "Yes".
% .exclusive [optional field] is an array of index in items indicating which 
%            options are mutually exclusive with each option. When the option 
%            is set to "Yes", mutually exclusive options are set to "No". If 
%            the field does not exist, or is empty, the control is a checkbox
%            otherwise it is a radio-button.
%            Both .linked and .exclusive behaviour are only implemented for
%            "Yes"/"No" fields - but any field can be in "linked" and will
%            be grayed out.
% .default   [optional field] is an integer indicating the default value for 
%            each option (0: No, 1: Yes, -1: grayed out). Note that in the 
%            case of text input field, this field is irrelevant and the default 
%            is given in the "values" field; In the case of popup lists, the
%            default is an index in .values, or -1 (grayed out).
%            "Linked" and "Exclusive" are NOT evaluated in the initial layout
%            of the dialog, hence the default must be consequent with the rules,
%            e.g. linked fields of a "No" field must be grayed out.
% .indent    [optional field] allows an indentation (proportional to the value
%            of this field) from the position of the control. May be used
%            to graphically make "groups" as no "frames" are used.
% .help      [optional field] contains tooltips help texts. Can contain
%            "\n" to make multi-line help fields.
% title     is the window title (string)
% msg       is a string appearing on the top of the dialog
% options   is an optional structure containing options:
%           .center = 0|1 (center text msg)
%           .bold = 'light'|'normal'|'demi'|'bold'
%                   indicates how headers (see .values) must be printed out.
%           .fixed = 0|1: FixedWidth font for lists
%           [more to come in future versions]
%
% The controls will be display in the order they appear in items.
%
% out contains an array of structure of answers:
%   out(i).answer = value of the control.
%   values are:
%       1 or 0 for "Yes"/"No" controls,
%       the index of the chosen item for list controls.
%       the text given for text inputs.
%   out is empty if cancel was chosen, negative integer if error.
%
% Written by L.Cavin, 07.12.2003, (c) CSE
% This code is free to use and modify for non-commercial purposes.
% Web address: http://ltcmail.ethz.ch/cavin/CSEDBLib.html#FLAGDIALOG
%
% =====================================================================
% An example of usage is: (in items, elements not mentionned are empty;
% e.g. items(1).linked is empty).
% 	items(1).name = 'Contact:';
% 	items(1).default = 0;
% 	items(1).values = {'email@address'};
% 	items(1).help = 'Enter your email address.';
% 	items(2).name = 'I will be coming!';
% 	items(2).default = 1;
% 	items(2).linked = [3 4 5 6];
% 	items(3).name = 'With my family';  
% 	items(3).default = 1;
% 	items(3).exclusive = 4;
% 	items(3).indent = 1;
% 	items(4).name = 'Alone';
% 	items(4).default = 0;
% 	items(4).exclusive = 3;
% 	items(4).indent = 1;
% 	items(5).name = 'Transportation:';
% 	items(5).indent = 1;
% 	items(5).values = {1};
% 	items(6).name = 'Coming by';
% 	items(6).default = 1;
% 	items(6).indent = 2;
% 	items(6).values = {'Train'; 'Bus'; 'Foot'};
% 	items(6).help = 'Cars are polluting.\nUse public transportation whenever possible!';
% 	items(7).name = 'I''ll sure give a phone call!';
% 	items(7).default = 0;
% 
%   title = 'Birthday party incsription';
% 
%   msg = sprintf(['Dear friends,\nAs you all know, I am turning 30 next april.\nThis ' ...
%   'seems a worthy occasion to party a bit!\n\nWill you be able to ' ...
%   'attend?']);
% 
%   out = CSEFlagDialog(items, title, msg)


persistent handles;

if ischar(items)
    % this ain't no cell, so must be a callback
    callback_type = items;
    out = [];
elseif ~isa(items, 'struct')
    % this is an error
    warn('CSE:Dialog', 'The first parameter must be a stucture of items.');
    out = -1;
    return;
else
    % we must generate the dialog
    callback_type = 'create_dialog';
end

% first we do the callbacks, then we handle the initial call.
switch callback_type
    case 'item_click' % this is an action on a control
        ctrl_idx = get(gcbo, 'UserData');
        % should we "zero" other controls?
        if ~isempty(handles.items(ctrl_idx).exclusive) & get(gcbo, 'Value')==1
            for i = 1:length(handles.items(ctrl_idx).exclusive)
                set(handles.items_obj(handles.items(ctrl_idx).exclusive(i)), 'Value', 0);
            end
        end
        % should we disable or enable other elements
        if ~isempty(handles.items(ctrl_idx).linked)
            if get(gcbo, 'Value')==1
                % enable
                en_st = 'on';
                ds_st = 'off';
            else
                % disable
                en_st = 'off';
                ds_st = 'on';
            end
            for i = 1:length(handles.items(ctrl_idx).linked)
                tp = handles.items(ctrl_idx).linked(i);
                if tp > 0
                    set(handles.items_obj(tp), 'Enable', en_st);
                    if handles.head_obj(tp) > 0
                        set(handles.head_obj(tp), 'Enable', en_st);
                    end
                else
                    tp = -1 * tp;
                    set(handles.items_obj(tp), 'Enable', ds_st);
                    if handles.head_obj(tp) > 0
                        set(handles.head_obj(tp), 'Enable', ds_st);
                    end
                end
            end
        end
    case 'ok_click' % this is a click on the OK button
        for i = 1:length(handles.items_obj)
            if strcmp(get(handles.items_obj(i), 'Style'), 'edit')
                out(i).answer = cell2mat(get(handles.items_obj(i), 'String'));
            else
                out(i).answer = get(handles.items_obj(i), 'Value');
            end
        end
        set(gcbo,'userdata',out);
    case 'cancel_click' % click on cancel button
        out = [];
        if ~isempty(handles)
            set(handles.ok_button,'userdata',out);
        else handles = get(0,'CurrentFigure');
            delete(handles);
        end
    case 'create_dialog' % now this is actually the initial call
        options = [];
        msg = '';
        title = 'CSE Flag Dialog';
        if nargin > 1
            title = varargin{1};
            if nargin > 2
                msg = varargin{2};
                if nargin > 3
                    options = varargin{3};
                end
            end
        end
        handles.options = options;
        % precompute width and height of dialog:
        % width is 10 + the longest chain of characters + 10;
        % height is 1 + the number of lines in msg + 2 + number of items*2
        %           + 2 + 2 + 2 (sum: ... + 9)
        longest_chain = 0;
        item_length = 0;
        for i = 1:length(items)
            if ~isfield(items(i), 'indent') | isempty(items(i).indent)
                items(i).indent = 0;
            end
            if ~isfield(items(i), 'help')
                items(i).help = [];
            end
            if ~isfield(items(i), 'linked')
                items(i).linked = [];
            end
            if ~isfield(items(i), 'default')
                items(i).default = [];
            end
            if ~isfield(items(i), 'exclusive')
                items(i).exclusive = [];
            end
            val_length = 0;
            if length(items(i).values) > 0
                for j = 1:length(items(i).values)
                    if length(items(i).values{j}) > val_length
                        val_length = length(items(i).values{j});
                    end
                end
                if length(items(i).values) == 1
                    val_length = val_length + 10;
                end
            end
            items(i).mxlgt = val_length;
            if length(items(i).name)+val_length+items(i).indent*4 > item_length
                item_length = length(items(i).name)+val_length+items(i).indent*4;
            end
            if length(items(i).name)+val_length+items(i).indent*4 > longest_chain-10
                longest_chain = length(items(i).name)+10+val_length+items(i).indent*4;
            end
        end
        handles.items = items;
        % a = [regexp(msg, '\n') length(msg)]; % for > R13
        a = [strfind(msg,sprintf('\n')) length(msg)];  % for R12
        for i = 2:length(a)
            if a(i)-a(i-1) > longest_chain
                longest_chain = a(i)-a(i-1);
            end
        end
        dial_width = max(longest_chain + 20, 50);
        if isfield(handles.options, 'fixed')
            % much wider on average...
            dial_width = dial_width * 1.25;
        end
        dial_height = 9 + length(a)-1 + length(items)*2;
        if length(msg) == 0
            dial_height = dial_height -2;
        end
        item_length = (dial_width-item_length-10)/2;
        
        % A) create window: (invisible, for now)
        handles.dialog = dialog(    'Visible',          'off', ...
                                    'Units',            'characters', ...
                                    'Position',         [10 10 dial_width dial_height], ...
                                    'CloseRequestFcn',  'CSEFlagDialog(''cancel_click'');', ...
                                    'WindowStyle',      'modal', ...
                                    'Name',             title);
        % COMPATIBILITY TOWARDS R12: repeating instructions for safety:
        set(handles.dialog, 'Units', 'characters');
        set(handles.dialog, 'Position', [10 10 dial_width dial_height]);
                
        % B) create buttons:
        handles.ok_button = uicontrol(  'Units',    'characters', ...
                                        'Parent',   handles.dialog, ...
                                        'Position', [dial_width/2-20 2 10 2], ...
                                        'Callback', 'CSEFlagDialog(''ok_click'');', ...
                                        'String',   'OK' );
        handles.cancel_button = uicontrol(  'Units',    'characters', ...
                                            'Parent',   handles.dialog, ...
                                            'Position', [dial_width/2+10 2 10 2], ...
                                            'Callback', 'CSEFlagDialog(''cancel_click'');', ...
                                            'String',   'Cancel' );
        
        % C) create items:
        for i = 1:length(items)
            pos = length(items)-i+1;
            if isempty(items(i).values)
                ftype = 'radiobutton';
                if isempty(items(i).exclusive)
                    ftype = 'checkbox';		    
                end
                handles.head_obj(i) = -1;
                handles.items_obj(i) = uicontrol(   'Style',    ftype, ...
                                                    'Parent',   handles.dialog, ...
                                                    'Units',    'characters', ...
                                                    'Position', [items(i).indent*4+item_length 4+2*pos length(items(i).name)+10 1], ...
                                                    'Callback', 'CSEFlagDialog(''item_click'');', ...
                                                    'UserData', i, ...
                                                    'String',   items(i).name );
                if ~isempty(items(i).default) & items(i).default > 0
                    if strcmp(get(handles.items_obj(i), 'Style'), 'popupmenu')
                        if items(i).default > length(get(handles.items_obj(i), 'String'))
                            warn('CSE:Dialog', 'Initial value impossible for listbox %s.', items(i).name);
                            items(i).default = length(get(handles.items_obj(i), 'String'));
                        end
                    end
                    set(handles.items_obj(i), 'Value', items(i).default);
                end
            else
                if length(items(i).values) == 1
                    if isnumeric(items(i).values{1})
                        ftype = 'none';
                    else
                        ftype = 'edit';
                    end
                else
                    ftype = 'popupmenu';
                end
                % needs a "header":
                handles.head_obj(i) = uicontrol(    'Style',    'text', ...
                                                    'Parent',   handles.dialog, ...
                                                    'Units',    'characters', ...
                                                    'Position', [items(i).indent*4+item_length 4+2*pos (length(items(i).name)+2)*1.5 1], ...
                                                    'String',   items(i).name, ...
                                                    'HorizontalAlignment', 'left');
                if ~strcmp(ftype, 'none')
                    handles.items_obj(i) = uicontrol(   'Style',    ftype, ...
                                                        'Parent',   handles.dialog, ...
                                                        'BackgroundColor', 'white', ...
                                                        'Units',    'characters', ...
                                                        'Position', [items(i).indent*4+item_length+length(items(i).name)+3 4+2*pos-0.2 items(i).mxlgt+10 1.4], ...
                                                        'Callback', 'CSEFlagDialog(''item_click'');', ...
                                                        'UserData', i, ...
                                                        'String',   items(i).values, ...
                                                        'HorizontalAlignment', 'left' );
                    if ~isempty(items(i).default) & items(i).default > 0
                        if strcmp(get(handles.items_obj(i), 'Style'), 'popupmenu')
                            if items(i).default > length(get(handles.items_obj(i), 'String'))
                                warn('CSE:Dialog', 'Initial value impossible for listbox %s.', items(i).name);
                                items(i).default = length(get(handles.items_obj(i), 'String'));
                            end
                        end
                        set(handles.items_obj(i), 'Value', items(i).default);
                    end
                    if isfield(handles.options, 'fixed')
                        set(handles.items_obj(i), 'FontName', 'FixedWidth');
                        set(handles.items_obj(i), 'Position', [items(i).indent*4+item_length+length(items(i).name)+3 4+2*pos-0.2 (items(i).mxlgt+10)*1.3 1.4]);
                    end
                else
                    % was actually the "header" :-)
                    if isfield(handles.options, 'bold')
                        bld = handles.options.bold;
                    else
                        bld = 'normal';
                    end
                    set(handles.head_obj(i), 'FontWeight', bld);
                    handles.items_obj(i) = handles.head_obj(i);
                end
            end
            if ~isempty(items(i).default) & items(i).default < 0
                set(handles.head_obj(i), 'Enable', 'off'); % RECENTLY ADDED
                set(handles.items_obj(i), 'Enable', 'off');
            end
            if ~isempty(items(i).help)
                set(handles.items_obj(i), 'Tooltip', sprintf(items(i).help));
            end
        end
        
        % D) create message:
        if length(msg) > 0
            if isfield(handles.options, 'center') & handles.options.center == 1
                algn = 'center';
            else
                algn = 'left';
            end
            handles.text = uicontrol(  'Units',                 'characters', ...
                                       'HorizontalAlignment',   algn, ...
                                       'Style',                 'text', ...
                                       'Parent',                handles.dialog, ...
                                       'Position',              [5 dial_height-length(a)-1 longest_chain+10 length(a)], ...
                                       'String',                msg );
        end
        
        % E) "center" and show dialog:
        movegui(handles.dialog, 'center');
        set(handles.dialog, 'Visible', 'on');
        drawnow % NB! for R14
        % F) run dialog:
        waitfor(handles.ok_button,'userdata');
        
        % G) finish call:
        out = get(handles.ok_button,'userdata');
        delete(handles.dialog);
        clear handles;
        
    otherwise
        % this is a mistake
        warn('CSE:Dialog', 'Unknown Callback. Even if there is only one option, it must be passed as a cell.');
        out = -1;
end

        
%=====HELPER FUNCTION
function warn(tag, msg, varargin)
% patch for the warning function to make R13 calls compatible with R12

if str2num(version('-release')) < 13
    warning(sprintf(msg, varargin{:}));
else
    warning(tag, msg, varargin{:});
end
        