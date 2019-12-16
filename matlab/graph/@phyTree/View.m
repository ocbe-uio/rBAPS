function View(tr,sel,propsForFigure)
%VIEW views a phylogenetic tree in phytreetool.
%
%   VIEW(TREE) shows a phylogenetic tree object. The significant distances
%   between branches and nodes are in horizontal direction, vertical
%   coordinates are accommodated only for display purposes. Tree
%   Edit/Analysis tools are accessible through the mouse left/right buttons
%   and also using the 'Tree' menu.  
%
%   VIEW(TREE,SEL) starts the viewer with an initial selection of nodes
%   specified by SEL. SEL can be a logical array of any of the following
%   sizes: [NUMLEAVES+NUMBRANCHES x 1], [NUMLEAVES x 1], or [NUMBRANCHES x
%   1]. SEL may also be a list of indices. 
% 
%   Examples:
%      
%       tr = phytreeread('pf00002.tree')
%       view(tr)
%
%   See also PHYTREE, PHYTREE/PLOT, PHYTREEREAD, PHYTREETOOL, SEQLINKAGE,
%   SEQNEIGHJOIN.

% Copyright 2003-2006 The MathWorks, Inc.
% $Revision: 1.1.6.19.2.1 $ $Author: batserve $ $Date: 2006/07/24 13:54:28 $

if numel(tr)~=1
     error('Bioinfo:phytree:view:NoMultielementArrays',...
           'Phylogenetic tree must be an 1-by-1 object.');
end

tr = doBasicCalculations(tr);

nodeIndex   = 1:tr.numLabels;
leafIndex   = 1:tr.numLeaves;
branchIndex = tr.numLeaves+1:tr.numLabels;

% check empty names
for ind = nodeIndex
    if isempty(tr.names{ind}) 
        if ind > tr.numLeaves
            tr.names{ind} = ['Branch ' num2str(ind-tr.numLeaves)]; 
        else
            tr.names{ind} = ['Leaf ' num2str(ind)]; 
        end
    end
end

% initial drawing
if nargin<3
    propsForFigure.Name = ['Phylogenetic Tree Tool ' getphytreetoolnumber];
end
propsForFigure.PruneWarning = getacceptedwarningfromothertools;
fig = figure('Renderer','ZBuffer','Name',propsForFigure.Name,...
           'NumberTitle','off','IntegerHandle','off','tag','PhyTreeTool');
setappdata(fig,'propsForFigure',propsForFigure)  
setappdata(fig,'backupTree',tr) 
tr.ha = axes; hold on;
set(tr.ha,'Position',[.05 .05 .7 .9],'YTick',leafIndex,'FontSize',9,'Ydir','reverse',...
          'YAxisLocation','Right','YTickLabel',char(tr.names{leafIndex}))
tr.hlines = plot( ...
  tr.x([nodeIndex;repmat([tr.par(1:tr.numLabels-1) tr.numLabels],2,1)]),...
  tr.y([repmat(nodeIndex,2,1);[tr.par(1:tr.numLabels-1) tr.numLabels]]),...
                '-k');
tr.hpathline = plot(1,1,'--r','LineWidth',2,'Visible','off');  
tr.hdragbox = plot(1,1,':k','LineWidth',1,'Visible','off'); 
tr.hdots(1,1) = plot(tr.x(branchIndex),tr.y(branchIndex),'o',...
               'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','b');
tr.hdots(1,2) = plot(tr.x(leafIndex),tr.y(leafIndex),'square',...
               'MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor','w');
tr.hseldots(1,1) = plot(tr.x(branchIndex),tr.y(branchIndex),'o',...
               'MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r');
tr.hseldots(1,2) = plot(tr.x(leafIndex),tr.y(leafIndex),'square',...
               'MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r');
tr.hldots(1,1) = plot(tr.x(branchIndex),tr.y(branchIndex),'o',...
               'MarkerSize',5,'MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor',[.6 .6 1]);
tr.hldots(1,2) = plot(tr.x(leafIndex),tr.y(leafIndex),'square',...
               'MarkerSize',4,'MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor','w');     
set(tr.hldots(1),'Xdata',[],'Ydata',[])
set(tr.hldots(2),'Xdata',[],'Ydata',[])
tr.axhold = plot([-eps -eps],[0 0],'.','MarkerSize',eps,'Color','w');
tr.datatip = text(0,1,1,'k','Tag','TreeTag','BackgroundColor',[1 1 .93],...
               'Color', [0 0 0],'EdgeColor', [0.8 0.8 0.8],...
               'VerticalAlignment','Top','Clipping','off',...
               'Visible','off','Fontsize',8,'Interpreter','none'); 

if nargin == 1 || isempty(sel)
   tr.selected = false(tr.numLabels,1);        % selected nodes
else
    % validate sel
    if islogical(sel)
        if numel(sel)==tr.numLabels 
            sel = sel(:)==true;
        elseif numel(sel)==tr.numLeaves
            sel = [sel(:);false(tr.numBranches,1)];
        elseif numel(sel)==tr.numBranches
            sel = [false(tr.numLeaves,1);sel(:)];
        else
        close(fig)    
        error('Bioinfo:phytree:view:IncorrectLogical',...
              'Logical vector must have the same number of elements as nodes in the Phylogenetic Tree');
        end
    elseif isnumeric(sel) && isreal(sel) && all(sel>=1) && all(sel<=tr.numLabels)
        tem(tr.numLabels)=false;
        tem(floor(sel))=true;
        sel=tem(:);
    else
        close(fig)    
        error('Bioinfo:phytree:view:IncorrectTypeofArguments',...
              'Invalid value for NODES');
    end
    tr.selected =sel;
end              
               
% save more figure data needed for the gui functionality
tr.activeNodes    = true(tr.numLabels,1);   % active nodes
tr.activeBranches = true(tr.numBranches,1); % active Branches
tr.sel2root = false(tr.numLabels,1);        % path sel-node to root
tr.editMode = 'Select';                            % initial edit mode
tr.indicativeMode = false;                         % data-tip flag
tr.lastThresholdValue = [];                        % remembers last cut

% create uicontrols (will appear as needed, initially invisible)
tr.editBox =  uicontrol(fig,'Background',[1 1 1],'style','edit',...
                        'visible','off','callback',@doneRenaming);                     
tr.slider  =  uicontrol(fig,'style','slider','SliderStep',[.1 .1],...
                        'visible','off','callback',@sliderCallback);
tr.slidertx = uicontrol(fig,'style','text','visible','off');
tr.sliderok = uicontrol(fig,'style','pushbutton','visible','off',...
                        'string','OK','callback',@doThresholdCut);

% setup callback for click over nodes                    
set([tr.hseldots,tr.hdots,tr.hldots],'ButtonDownFcn',@toggleNode)
% setup figure callback functions
set(fig,'WindowButtonDownFcn',@mouseClickOnFigure);
set(fig,'WindowButtonUpFcn',@mouseRelease);
set(fig,'WindowButtonMotionFcn',@localWindowButtonMotion);

% setup UIMenus, context menus and toolbar
tr.hToggleUIMenu  = makePhyTreeViewerUIMenus(fig);                    
tr.hToggleToolbar = makePhyTreeViewerToolbar(fig);
[tr.hToggleContextMenu,tr.hAxisContextMenu,tr.hDotsContextMenu] = ...
                                        makePhyTreeViewerContextMenus(fig);
% activate Context Menus
set(tr.ha,'UIContextMenu',tr.hAxisContextMenu);
set([tr.hdots tr.hldots tr.hseldots],'UIContextMenu',tr.hDotsContextMenu);

set(fig,'UserData',tr)           % save figure data

correctFigureSize(fig, 15 * tr.numLeaves);        % resize figure if needed
setupYLabelsListeners;           % listeners for YLabels
updateTree(fig,[],[])            % updates figure after all initializations 
set(gca,'xLim',[0  max(tr.x)] + max(tr.x) * [-.1 .05]); 
tr.yLim = get(tr.ha,'Ylim');tr.xLim = get(tr.ha,'Xlim');
set(fig,'UserData',tr)           % save figure data
toolsmenufcn(fig,'PanY')         % set zoom mode to vertical constraining
toolsmenufcn(fig,'ZoomY')        % set pan  mode to vertical constraining
set(fig,'HandleVisibility','callback') % after all init, make it invisible

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ylabelsListener(hSrc,event,hf,ha) %#ok
% Auto sizes the ylabels
ratio = max(get(hf,'Position').*[0 0 0 1])/diff(get(ha,'YLim'));
set(ha,'Fontsize',min(9,ceil(ratio/1.7)));    % the gold formula
% Also verify if we need to re-position the slidebar of threshold cut
tr=get(hf,'Userdata');
if any(strcmp(tr.editMode,{'Distance to Leaves','Distance to Root'}))
    wS = get(hf,'Position');  % window dimensions
    aP = get(tr.ha,'Position'); % axes position
    set(tr.slider,  'Position',[aP(1)*wS(3) wS(4)-20 aP(3)*wS(3) 20])
    set(tr.slidertx,'Position',[sum(aP([1 3]))*wS(3) wS(4)-20 60 20])
    set(tr.sliderok,'Position',[sum(aP([1 3]))*wS(3)+60 wS(4)-20 30 20])
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseClickOnFigure(h,varargin)
% This callback function is activated when a mouse button is pressed in any
% location of the figure and under any of my edit modes
tr = get(gcbf,'Userdata');
switch tr.editMode
    case 'Renaming';           doneRenaming(h,varargin);
    case 'Distance to Leaves'; cancelThresholdCut(h,varargin);
    case 'Distance to Root';   cancelThresholdCut(h,varargin);
    case 'Select';       
        switch get(gcbf,'SelectionType')
            case {'normal','extend'}
                tr = get(gcbf,'userdata');
                cp = get(tr.ha,'CurrentPoint');
                xPos = cp(1,1); yPos = cp(1,2); 
                set(tr.hdragbox,'Visible','on',...
                    'Xdata',repmat(xPos,5,1),'Ydata',repmat(yPos,5,1))
            case 'open'
                autoFit(h)
        end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggleNode(h,varargin)
% This callback function is activated when a mouse button is pressed over
% any of the displayed nodes under any of my edit modes
hideActiveIndicators(h,varargin)
tr = get(gcbf,'Userdata');
switch get(gcbf,'SelectionType')
    case 'normal'
        switch tr.editMode
            case 'Select';             selectNode(h,varargin);
            case 'Inspect';            inspectNode(h,varargin);
            case 'Collapse/Expand';    collapseExpand(h,varargin);
            case 'Rotate Branch';      rotateBranch(h,varargin);
            case 'Rename';             renameNode(h,varargin);
            case 'Renaming';           doneRenaming(h,varargin);
            case 'Prune';              pruneTree(h,varargin);
            case 'Distance to Leaves'; cancelThresholdCut(h,varargin);    
            case 'Distance to Root';   cancelThresholdCut(h,varargin);        
        end
    case 'extend'
        switch tr.editMode
            case 'Select';             selectNode(h,varargin);
        end
    case 'alt'
    case 'open'    
 end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeEditMode(h,varargin) %#ok
% Callback function to change the edit mode, this function is
% called from the toolbar, the context menu or the uimenu.
tr = get(gcbf,'Userdata');
myModes = {'Inspect','Collapse/Expand','Rotate Branch','Rename','Prune'};

% first, disable any present edit mode
switch tr.editMode
    case {myModes{:},'Select'};  
                      disableMyContextMenus(h)
                      disableMyWindowButtonActions(h)
                      ind = strmatch(tr.editMode,myModes); %#ok
                      set(tr.hToggleToolbar(ind),    'State','off')
                      set(tr.hToggleUIMenu(ind),     'Checked','off')
                      set(tr.hToggleContextMenu(ind),'Checked','off')
    %case '&Zoom In';  toolsmenufcn(gcbf,'ZoomIn');
    case '&Zoom In';  zoom(gcbf,'off')
    %case 'Zoom &Out'; toolsmenufcn(gcbf,'ZoomOut');
    case 'Zoom &Out'; zoom(gcbf,'off')
    %case '&Pan';      toolsmenufcn(gcbf,'Pan');
    case '&Pan';      pan(gcbf,'off');
    case {'Distance to Leaves','Distance to Root'}
                      enableAllUI(h)
                      disableMyWindowButtonActions(h)
end
 
% depending on the caller instance, determine the new edit mode
switch get(h,'Type')
    case 'uimenu'; newEditMode = get(h,'Label');
    case 'uitoggletool'
        newEditMode = get(h,'Tag');
        switch newEditMode
            case 'Exploration.ZoomIn';  newEditMode = '&Zoom In';
            case 'Exploration.ZoomOut'; newEditMode = 'Zoom &Out';
            case 'Exploration.Pan';     newEditMode = '&Pan';
        end
    otherwise; newEditMode = 'Select';
end
%disp( [tr.editMode ' -->  ' newEditMode]  )
% if new mode is the same then we are toggling off
if strcmp(newEditMode,tr.editMode) 
    newEditMode = 'Select'; 
end

% if changing to Prune, verify the warnign has been accepted
if strcmp(newEditMode,'Prune') 
    propsForFigure = getappdata(gcbf,'propsForFigure');
    if isequal(propsForFigure.PruneWarning,'NotDone')
        warndlg(['Pruning nodes cannot be undone. Before continuing,',...
             ' you may want to export the current tree to a new tool.'],...
             'Warning','modal')
        setacceptedwarningtoothertools
    end
end

switch newEditMode
    case '&Zoom In';  toolsmenufcn(gcbf,'ZoomIn'); 
    case 'Zoom &Out'; toolsmenufcn(gcbf,'ZoomOut');
    %case '&Pan';      toolsmenufcn(gcbf,'Pan'); 
    case '&Pan';      pan(gcbf,'on');
    case myModes;     enableMyContextMenus(h)
                      enableMyWindowButtonActions(h)
                      ind = strmatch(newEditMode,myModes); %#ok
                      set(tr.hToggleToolbar(ind),    'State','on')
                      set(tr.hToggleUIMenu(ind),     'Checked','on')
                      set(tr.hToggleContextMenu(ind),'Checked','on')
    case 'Select';    enableMyContextMenus(h)
                      enableMyWindowButtonActions(h)
    case {'Distance to Leaves','Distance to Root'}
                      disableAllUI(h)
                      enableMyWindowButtonActions(h)
                      set(gcbf,'WindowButtonMotionFcn',[])
                      set(gcbf,'WindowButtonUpFcn',[])
                      set([tr.hseldots,tr.hdots,tr.hldots],'ButtonDownFcn',[])
end

switch newEditMode
    case 'Inspect';  if sum(tr.selected(:)) ~= 1
                         tr.selected(:) = false;
                         tr.selected(end) = true;
                     end
    otherwise
end    

tr.sel2root = path2root(tr, tr.selected);
tr.editMode = newEditMode;
set(gcbf,'userdata',tr)
updateTree(gcbf,[],[])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hideActiveIndicators(h,varargin) %#ok
tr = get(gcbf,'userdata');
set([tr.hpathline,tr.datatip],'visible','off')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disableAllUI(h,varargin) %#ok
hw = findall(gcbf,'Type','uimenu','Parent',gcbf);
gw = findall(gcbf,'Type','UIToggleTool');
set([hw;gw],'Enable','off')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enableAllUI(h,varargin) %#ok
hw = findall(gcbf,'Type','uimenu','Parent',gcbf);
gw = findall(gcbf,'Type','UIToggleTool');
set([hw;gw],'Enable','on')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disableMyContextMenus(h,varargin) %#ok
tr = get(gcbf,'userdata');
set(tr.ha,'UIContextMenu',[]);
set([tr.hdots tr.hseldots],'UIContextMenu',[]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enableMyContextMenus(h,varargin) %#ok
tr = get(gcbf,'userdata');
set(tr.ha,'UIContextMenu',tr.hAxisContextMenu);
set([tr.hdots tr.hldots tr.hseldots],'UIContextMenu',tr.hDotsContextMenu);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disableMyWindowButtonActions(h,varargin) %#ok
set(gcbf,'WindowButtonDownFcn',[]);
set(gcbf,'WindowButtonUpFcn',[])
set(gcbf,'WindowButtonMotionFcn',[]);
tr = get(gcbf,'userdata');
set([tr.hseldots,tr.hdots,tr.hldots],'ButtonDownFcn',[])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enableMyWindowButtonActions(h,varargin) %#ok
set(gcbf,'WindowButtonDownFcn',@mouseClickOnFigure);
set(gcbf,'WindowButtonUpFcn',@mouseRelease);
set(gcbf,'WindowButtonMotionFcn',@localWindowButtonMotion);
tr = get(gcbf,'userdata');
set([tr.hseldots,tr.hdots,tr.hldots],'ButtonDownFcn',@toggleNode)
set(gcbf,'KeyPressFcn',[]);
set(gcbf,'Pointer','arrow');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localWindowButtonMotion(h,varargin) %#ok
% Callback function activated when moving over the axes, checks location of
% the mouse and puts datatip if over an active node.

tr = get(h,'userdata');
% set a virtual grid to get the point
xThres=diff(get(tr.ha,'Xlim'))/100;
yThres=diff(get(tr.ha,'Ylim'))/100;
cp = get(tr.ha,'CurrentPoint');
xPos = cp(1,1); yPos = cp(1,2);
hp = tr.x<(xPos+xThres) & tr.x>(xPos-xThres) & ...
     tr.y<(yPos+yThres) & tr.y>(yPos-yThres);
hp = find (hp & tr.activeNodes);

% shortcut out when dragging a box in select mode
if strcmp(get(tr.hdragbox,'Visible'),'on')
    xdata = get(tr.hdragbox,'XData');xdata([3,4]) = xPos;
    ydata = get(tr.hdragbox,'YData');ydata([2,3]) = yPos;
    set(tr.hdragbox,'XData',xdata,'YData',ydata)
% shortcut out when turning off 'indicative' mode
elseif tr.indicativeMode && isempty(hp) %&& isempty(tr.highligth)
    set([tr.datatip tr.hpathline],'visible','off')
    set(tr.hlines,'color','black')
    set(tr.hldots(1),'Xdata',[],'Ydata',[])
    set(tr.hldots(2),'Xdata',[],'Ydata',[])
% turn on or update 'indicative' mode
elseif numel(hp) % && isempty(tr.highligth)
   % find leaves (children) below this branch
   children = false(1,tr.numLabels);
   children(hp(1)) = true;
   for ind = hp(1)-tr.numLeaves:-1:1
       if children(ind+tr.numLeaves)
           children(tr.tree(ind,:))=true;
       end
   end 
   
   % find and draw path to selected
   if strcmp(tr.editMode,'Inspect')
       [pathA,pathB] = path2sel(tr,hp(1)); 
       dis2sel = tr.x(find(pathA,1))+tr.x(find(pathB,1))...
                 -2*tr.x(find(pathA,1,'last'));
       if any(pathB)
          xx = [tr.x(pathA);NaN;tr.x(pathB)];
          yy = [tr.y(pathA);NaN;tr.y(pathB)];
          hh=zeros(2*numel(xx),1); hh(1:2:end)=1; hh=cumsum(hh);
          set(tr.hpathline,'XData',xx(hh(2:end)),...
             'YData',yy(hh(1:end-1)),'Visible','on');
       end
   end
   
   % place text 
   name = [tr.names{hp(1)} ' '];
   name(name=='_')=' ';
   children(hp(1)) = false;
   numChil = sum(children(1:tr.numLeaves));
   childrenNames = char(tr.names(children(1:tr.numLeaves)));
   childrenNames(childrenNames=='_')=' ';
   childrenNames=[repmat('   ',size(childrenNames,1),1) childrenNames];
   switch tr.editMode
       case 'Inspect'
           if numChil
               set(tr.datatip,'string',char( [ {name; ...
                   ['Dist to parent: '   num2str(tr.dist(hp(1)))];...
                   ['Dist to root:     ' num2str(tr.x(hp(1))-tr.x(end))];...
                   ['Path length:     '  num2str(dis2sel)];...
                   ['Samples: ' num2str(numChil)]};...
                    mat2cell(childrenNames,ones(size(childrenNames,1),1),...
                    size(childrenNames,2))]))
                   extraLines = 5;
           else
                   set(tr.datatip,'string',char(  {name; ...
                   ['Dist to parent: '   num2str(tr.dist(hp(1)))];...
                   ['Dist to root:     ' num2str(tr.x(hp(1))-tr.x(end))];...
                   ['Path length:     '  num2str(dis2sel)]}))
                   extraLines = 4;
           end
       case {'Collapse/Expand','Rotate Branch','Rename','Prune','Select'}
           if numChil
               set(tr.datatip,'string',char([ 
                   {[name '  (' num2str(numChil) ' samples)']};...
                   mat2cell(childrenNames,ones(size(childrenNames,1),1),...
                   size(childrenNames,2))]))
           else
               set(tr.datatip,'string', name)
           end
           extraLines = 1;
       otherwise % all other modes
   end

   %compute some values before adjusting data tip
   fp = get(gcbf,'Position'); % fig position in points
   fh = fp(4);%fw = fp(3);     % fig size (height & width) in points
   ap  = get(tr.ha,'Position');          % axis position normalized
   yl  = ylim(tr.ha); yl  = yl - ...
         [ap(2) ap(2)+ap(4)-1]*diff(yl)/ap(4); % fig height limits in axis units
   xl  = xlim(tr.ha); xl  = xl - ...
         [ap(1) ap(1)+ap(3)-1]*diff(xl)/ap(3); % fig width limits in axis units
   yPosPt = (-yPos -4*yThres + yl(2))*fh/diff(yl); % datatip position in pts
   reqPt  = (numChil+extraLines)*14+2;  % required datatip height in pts
                                        % adjust if other fontsize is used 
                                         
   %adjust string of datatip if it will not fit (i.e. remove names)
   if reqPt > fh
       str = get(tr.datatip,'String');
       set(tr.datatip,'String',str(1:extraLines,:));
       reqPt = extraLines*14+2;
   end
   
   %adjust vertical position of datatip just below cp
   topEdge = yl(2)-min(fh,max(yPosPt,reqPt))*diff(yl)/fh;
   switch tr.editMode
       case {'Collapse/Expand','Rotate Branch','Prune'} 
       % datatip usually to the left of cp to see shadowing of branches    
           datatipExtent = get(tr.datatip,'Extent');
           datatipWidth = datatipExtent(3);
           rightEdge = max(xPos-3*xThres,xl(1)+datatipWidth);
           % is the datatip over cp ?
           if rightEdge>xPos && topEdge<yPos
               % then try to put it above cp
                topEdge = yPos - 3 * yThres - reqPt*diff(yl)/fh; 
               % does datatip fit above cp ?
               if topEdge<yl(1)
                   % then adjust string by removing names of species
                   str = get(tr.datatip,'String');
                   set(tr.datatip,'String',str(1:extraLines,:));
                   reqPt = extraLines*14+2;
                   topEdge = yl(2)-min(fh,max(yPosPt,reqPt))*diff(yl)/fh;
               end
           end
           set(tr.datatip,'Position',[rightEdge-datatipWidth,topEdge,1])
           set(tr.datatip,'horizontalalignment','left')
       case {'Inspect','Rename','Select'}
       % datatip usually to the right of cp to minimize problems on the
       % left edge
           set(tr.datatip,'Position',[xPos+3*xThres,topEdge,1])
           set(tr.datatip,'horizontalalignment','left')
       otherwise
   end

   switch tr.editMode
       case {'Collapse/Expand','Rotate Branch','Prune'}
           % de-color branches to rotate or collapse
           uncoloredNodes = false(1,tr.numLabels);
           uncoloredNodes(hp(1)) = true;
           for ind =  hp(1)-tr.numLeaves:-1:1
               if uncoloredNodes(ind+tr.numLeaves)
                   uncoloredNodes(tr.tree(ind,:))=true;
               end
           end
           if ~strcmp(tr.editMode,'Prune')
               uncoloredNodes(hp(1)) = false;
           end
           uncoloredNodes = uncoloredNodes & tr.activeNodes';
           set(tr.hlines(uncoloredNodes),'color',[.87 .87 .87])
           ind=find(uncoloredNodes(tr.numLeaves+1:tr.numLabels))+tr.numLeaves;
           set(tr.hldots(1),'Xdata',tr.x(ind),'Ydata',tr.y(ind))
           ind=find(uncoloredNodes(1:tr.numLeaves));
           set(tr.hldots(2),'Xdata',tr.x(ind),'Ydata',tr.y(ind))
       otherwise
   end
  
   set(tr.datatip,'Visible','on')
   tr.indicativeMode = true;
   set(gcbf,'userdata',tr)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseRelease(h,varargin) %#ok
tr = get(gcbf,'userdata');
xdata = get(tr.hdragbox,'Xdata');
ydata = get(tr.hdragbox,'Ydata');
hp = (tr.x<max(xdata) & tr.x>min(xdata) & ...
          tr.y<max(ydata) & tr.y>min(ydata) & tr.activeNodes) ;
if (strcmp(get(gcbf,'SelectionType'),'normal') && ...
    strcmp(get(tr.hdragbox,'visible'),'on'))
    tr.selected(:) = false;
end
tr.selected(hp) = true;
tr.sel2root = path2root(tr,tr.selected);
set(tr.hdragbox,'Visible','off')
set(gcbf,'userdata',tr)
updateTree(gcbf,[],[])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectNode(h,varargin) %#ok
% Callback function to select a Node. 
% Entry points: from 1) the dots context menu or 2) toggle node
tr = get(gcbf,'userdata');
% set a virtual grid to get the point
xThres=diff(get(tr.ha,'Xlim'))/100;
yThres=diff(get(tr.ha,'Ylim'))/100;
cp = get(tr.ha,'CurrentPoint');
xPos = cp(1,1); yPos = cp(1,2);
hp = find(tr.x<(xPos+xThres) & tr.x>(xPos-xThres) & ...
          tr.y<(yPos+yThres) & tr.y>(yPos-yThres)) ;

if numel(hp)
     set(tr.hdragbox,'visible','off')
     temp = tr.selected(hp(1));
     switch get(gcbf,'SelectionType')
         case 'normal'; tr.selected(:) = false;
         case 'alt'; if ~temp 
                        tr.selected(:) = false; 
                     end
                     temp=false;
     end
     tr.selected(hp(1)) = ~temp;
     tr.sel2root = path2root(tr,tr.selected);
     set(gcbf,'userdata',tr)
     updateTree(gcbf,[],[])
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inspectNode(h,varargin) %#ok
% Callback function to inspect the reference Node. 
% Entry points: from 1) the dots context menu or 2) toggle node
tr = get(gcbf,'userdata');
% set a virtual grid to get the point
xThres=diff(get(tr.ha,'Xlim'))/100;
yThres=diff(get(tr.ha,'Ylim'))/100;
cp = get(tr.ha,'CurrentPoint');
xPos = cp(1,1); yPos = cp(1,2);
hp = find(tr.x<(xPos+xThres) & tr.x>(xPos-xThres) & ...
          tr.y<(yPos+yThres) & tr.y>(yPos-yThres)) ;
      
if numel(hp)
     temp = tr.selected(hp(1));
     tr.selected(:) = false;
     tr.selected(hp(1)) = ~temp;
     tr.sel2root = path2root(tr,tr.selected);
     set(gcbf,'userdata',tr)
     updateTree(gcbf,[],[])
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function collapseExpand(h,varargin) %#ok
% Callback function to Collapse/Expand a branch. 
% Entry points: from 1) the dots context menu or 2) toggle node

tr = get(gcbf,'userdata');
if strcmp(get(h,'Type'),'uimenu') % come from a context menu
    hp = find(tr.selected(tr.numLeaves+1:tr.numLabels));
else
    % set a virtual grid to get the point
    xThres=diff(get(tr.ha,'Xlim'))/100;
    yThres=diff(get(tr.ha,'Ylim'))/100;
    cp = get(tr.ha,'CurrentPoint');
    xPos = cp(1,1); yPos = cp(1,2);
    hp = find(tr.x<(xPos+xThres) & tr.x>(xPos-xThres) & ...
              tr.y<(yPos+yThres) & tr.y>(yPos-yThres)) ;
    hp=hp(hp>tr.numLeaves)-tr.numLeaves;
    if numel(hp)
        hp=hp(1); %just in case it picked two points
    end
end

if numel(hp)
    for ind = 1:numel(hp)
        tr.activeBranches(hp(ind))=~tr.activeBranches(hp(ind));
        activeBranches=find(tr.activeBranches)';
        % find active nodes by expanding active Branches
        tr.activeNodes(:)=false;
        tr.activeNodes(tr.numLabels,1)=true;
        for k = activeBranches(end:-1:1)
            tr.activeNodes(tr.tree(k,:))=tr.activeNodes(k+tr.numLeaves);
        end
        tr.selected(:) = false;
        tr.selected(hp(ind)+tr.numLeaves) = true;
    end
    tr.sel2root = path2root(tr,tr.selected);
    set(gcbf,'userdata',tr)
    updateTree(gcbf,[],hp(end)+tr.numLeaves)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotateBranch(h,varargin) %#ok
% Callback function to rotate a branch reordering the leaves. 
% Entry points: from 1) the dots context menu or 2) toggle node
tr = get(gcbf,'userdata');
if strcmp(get(h,'Type'),'uimenu') % come from a context menu
    hp = find(tr.selected(tr.numLeaves+1:tr.numLabels));
else
    % set a virtual grid to get the point
    xThres=diff(get(tr.ha,'Xlim'))/100;
    yThres=diff(get(tr.ha,'Ylim'))/100;
    cp = get(tr.ha,'CurrentPoint');
    xPos = cp(1,1); yPos = cp(1,2);
    hp = find(tr.x<(xPos+xThres) & tr.x>(xPos-xThres) & ...
        tr.y<(yPos+yThres) & tr.y>(yPos-yThres)) ;
    hp=hp(hp>tr.numLeaves)-tr.numLeaves;
    if numel(hp)
        hp=hp(1); %just in case it picked two points
    end
end
if numel(hp)
    for ind = 1:numel(hp)
        %find Leaves for every child
        childrenA = false(1,tr.numLabels);
        childrenA(tr.tree(hp(ind),1)) = true;
        for k = tr.tree(hp(ind),1)-tr.numLeaves:-1:1
            if childrenA(k+tr.numLeaves)
                childrenA(tr.tree(k,:))=true;
            end
        end
        childrenB = false(1,tr.numLabels);
        childrenB(tr.tree(hp(ind),2)) = true;
        for k = tr.tree(hp(ind),2)-tr.numLeaves:-1:1
            if childrenB(k+tr.numLeaves)
                childrenB(tr.tree(k,:))=true;
            end
        end
        permuta = 1:tr.numLabels;
        chA = find(childrenA(1:tr.numLeaves));
        chB = find(childrenB(1:tr.numLeaves));
        if chA(1)<chB(1)
            permuta([chA chB])=[chB chA];
        else
            permuta([chB chA])=[chA chB];
        end
        ipermuta = zeros(1,tr.numLabels);
        ipermuta(permuta)=1:tr.numLabels;
        tr.names = tr.names(permuta);
        tr.dist = tr.dist(permuta);
        tr.tree = ipermuta(tr.tree);
        tr.par = tr.par(permuta(1:end-1));
        tr.selected = tr.selected(permuta);
        tr.activeNodes = tr.activeNodes(permuta);
        tr.sel2root = tr.sel2root(permuta);
    end
    set(gcbf,'userdata',tr)
    updateTree(gcbf,[],[])
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function renameNode(h,varargin) %#ok
% Renames a Node. Puts an uicontrol to input the new name.
% Entry points: from 1) the dots context menu or 2) toggle node
tr = get(gcbf,'userdata');
Xlim=get(tr.ha,'Xlim');Ylim=get(tr.ha,'Ylim');
aPos=get(tr.ha,'Position');
% set a virtual grid to get the point
xThres=diff(Xlim)/100;
yThres=diff(Ylim)/100;
cp = get(tr.ha,'CurrentPoint');
xPos = cp(1,1); yPos = cp(1,2);
hp = find(tr.x<(xPos+xThres) & tr.x>(xPos-xThres) & ...
          tr.y<(yPos+yThres) & tr.y>(yPos-yThres)) ;
if numel(hp)
    tr.previousMode = tr.editMode; tr.editMode = 'Renaming';
    xBoxPos = (.02+(xPos-Xlim(1))/diff(Xlim))*aPos(3)+aPos(1);
    yBoxPos = (.02+(Ylim(2)-yPos)/diff(Ylim))*aPos(4)+aPos(2);
    position=get(gcbf,'position');
    position=[position(3)*xBoxPos position(4)*yBoxPos 150 20];
    set(tr.editBox,'position',position);
    set(tr.editBox,'Visible','on','string',tr.names{hp(1)},'Value',hp(1))
    disableAllUI(h)
    disableMyContextMenus(h)
    set(gcbf,'WindowButtonMotionFcn',[]); % disable windows mouse motion
    set(gcbf,'userdata',tr)
 end
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doneRenaming(h,varargin) %#ok
% Output helper function to abandon the "Renaming" mode
tr = get(gcbf,'userdata');
tr.editMode = tr.previousMode;
tr.names{get(tr.editBox,'Value')} = get(tr.editBox,'String');
set(tr.editBox,'Visible','off')
set(gcbf,'userdata',tr)
updateTree(gcbf,[],[]);
enableAllUI(h)
enableMyContextMenus(h)
enableMyWindowButtonActions(h)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pruneTree(h,varargin) %#ok
% Callback function to prune the tree, this function is complex because not
% only the basic structure is updated but all the other handles are
% also updated to contain the new tree.

% if changing to Prune, verify the warning has been accepted
propsForFigure = getappdata(gcbf,'propsForFigure');
if isequal(propsForFigure.PruneWarning,'NotDone')
     warndlg(['Pruning nodes cannot be undone. Before continuing,',...
            ' you may want to export the current tree to a new tool.'],...
            'Warning','modal')
     setacceptedwarningtoothertools
     return % do not do this pruning
end

tr = get(gcbf,'userdata');
if strcmp(get(h,'Type'),'uimenu') % comes from a context menu
    hp = find(tr.selected);
else
    tr = get(gcbf,'userdata');
    Xlim=get(tr.ha,'Xlim');Ylim=get(tr.ha,'Ylim');
    % aPos=get(tr.ha,'Position');
    % set a virtual grid to get the point
    xThres=diff(Xlim)/100;
    yThres=diff(Ylim)/100;
    cp = get(tr.ha,'CurrentPoint');
    xPos = cp(1,1); yPos = cp(1,2);
    hp = find(tr.x<(xPos+xThres) & tr.x>(xPos-xThres) & ...
        tr.y<(yPos+yThres) & tr.y>(yPos-yThres));
    hp=hp(1); %just in case it picked two points
end

hp(hp==tr.numLabels)=[]; %cannot delete the root

while numel(hp)
    %find all nodes to purge (i.e. all descendants)
    children = false(1,tr.numLabels);
    children(hp(1)) = true;
    for k = hp(1)-tr.numLeaves:-1:1
        if children(k+tr.numLeaves)
            children(tr.tree(k,:))=true;
        end
    end
    mypar = tr.par(hp(1));                                     % parent
    if mypar < tr.numLabels  % my parent is NOT the root
        % connect brother to granparent
        mygrpar = tr.par(mypar);                                 % grandparent
        myuncle = setxor(tr.tree(mygrpar-tr.numLeaves,:),mypar); % uncle
        mybro = setxor(tr.tree(mypar-tr.numLeaves,:),hp(1));   % brother
        tr.tree(mygrpar-tr.numLeaves,:) = [myuncle mybro];
        tr.dist(mybro) = tr.dist(mybro) + tr.dist(mypar);
        temp = get(tr.hlines(mygrpar-tr.numLeaves),'Xdata');
        temp([1 4])=tr.x(tr.tree(mygrpar-tr.numLeaves,:));
        set(tr.hlines(mygrpar-tr.numLeaves),'Xdata',temp);
        highlight = [mybro,mygrpar];
    else % if my parent is the root, now I am the new root
        temp=cell2mat(get(tr.hlines,'Xdata'))-tr.dist(end);
        for k = 1:tr.numBranches
            set(tr.hlines(k),'Xdata',temp(k,:));
        end
        highlight = setxor(tr.tree(mypar-tr.numLeaves,:),hp(1));
    end
    children(mypar) = true; %also delete my par
    % find indexes to change tree
    permuta = 1:tr.numLabels;
    permuta(children) = [];
    ipermuta = zeros(1,tr.numLabels);
    ipermuta(permuta) = 1:length(permuta);
    permutaBranches = permuta(permuta>tr.numLeaves)-tr.numLeaves;
    % update all tree structure fields
    tr.names = tr.names(permuta);
    tr.dist = tr.dist(permuta);
    tr.tree = tr.tree(permutaBranches,:);
    tr.tree = ipermuta(tr.tree);
    if isempty(tr.tree) 
        return; 
    end % one leaf, no branches !
    tr = doBasicCalculations(tr);
    hlines = tr.hlines;
    tr.hlines = tr.hlines(permuta);
    delete(setxor(hlines,tr.hlines));
    tr.activeNodes = tr.activeNodes(permuta);
    tr.activeBranches = tr.activeBranches(permutaBranches);
    tr.selected = tr.selected(permuta);
    tr.selected(:) = false;
    tr.selected(ipermuta(highlight)) = true;
    tr.sel2root = path2root(tr,tr.selected);
    % update the vector with nodes to prune (node index has changed)
    hp=ipermuta(hp);
    hp(1)=[];
    hp(hp==0)=[];
    set(gcbf,'userdata',tr)
    updateTree(gcbf,[],[])
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function colorDown(h,varargin) %#ok
% Color Down a branch.
% Entry points: from 1) the dots context menu or 2) toggle node
disp('Color Down a branch. Not implemented YET !!!')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function findNode(h,varargin) %#ok
treefig = gcbf;
tr = get(treefig,'userdata');
s = inputdlg('Regular Expression to match ?','Find Leaf/Branch',1);
if ~isempty(s)
    hc=regexpi(regexprep(tr.names,'_',' '),s);
    h = false(1,tr.numLabels);
    for ind = 1:tr.numLabels
        if ~isempty(hc{ind})
            h(ind)=true;
        end
    end
    hf = find(h);
    for ind = 1:length(hf)
        while ~tr.activeNodes(hf(ind))
            hf(ind)=tr.par(hf(ind));
        end
    end
    tr.selected(:) = false;
    tr.selected(hf) = true;
    tr.sel2root = path2root(tr,tr.selected); % update path to root
    set(treefig,'Userdata',tr);
    updateTree(treefig,[],[]) 

    % if selected are out of current view then fit the tree
    if (any(min(ylim(tr.ha))>tr.y(tr.selected)) || ...
        any(max(ylim(tr.ha))<tr.y(tr.selected)))
        autoFit(h)
    end
    
end
       
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function thresholdCut(h,varargin) %#ok
% sets the slider uicontrol and wait for the threshold cut to be entered.
expandAll(h);
autoFit(h);
changeEditMode(h)           % turns off any mode
tr = get(gcbf,'userdata');  
wS = get(gcbf,'Position');  % window dimensions
aP = get(tr.ha,'Position'); % axes position
set(tr.slider,  'Position',[aP(1)*wS(3) wS(4)-20 aP(3)*wS(3) 20])
set(tr.slidertx,'Position',[sum(aP([1 3]))*wS(3) wS(4)-20 60 20])
set(tr.sliderok,'Position',[sum(aP([1 3]))*wS(3)+60 wS(4)-20 30 20])
set([tr.slider tr.slidertx tr.sliderok],'visible','on')
if isempty(tr.lastThresholdValue)
    set(tr.slider,'max',max(tr.x),'value',max(tr.x)*.75)
else
    set(tr.slider,'max',max(tr.x),'value',tr.lastThresholdValue)
end
set(tr.slider,'TooltipString',get(h,'Label'))
sliderCallback(h,varargin)
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sliderCallback(h,varargin) %#ok
% this helper function serves both sliders (from root and from leaves), it
% hides selected nodes based on the current threshold cut
tr = get(gcbf,'userdata');
Value = get(tr.slider,'Value');
switch get(tr.slider,'TooltipString')
    case 'Distance to Leaves'
        tr.tocollapse = tr.dist2Leaf<(max(tr.x)-Value);
        set(tr.slidertx,'String',num2str(max(tr.x)-Value))
    case 'Distance to Root'
        tr.tocollapse = tr.x >= Value;
        set(tr.slidertx,'String',num2str(Value))
end
set(gcbf,'Userdata',tr);
toshow = [tr.tocollapse(tr.par(1:tr.numLabels-1));0]&tr.tocollapse;
mask = (1:tr.numLabels)'>tr.numLeaves;
% update light lines
set(tr.hlines(~toshow),'color','k')
set(tr.hlines(toshow),'color',[.87 .87 .87])
% update light dots
set(tr.hldots(1),'Ydata',tr.y(toshow&mask&tr.activeNodes),...
                 'Xdata',tr.x(toshow&mask&tr.activeNodes))
set(tr.hldots(2),'Ydata',tr.y(toshow&~mask&tr.activeNodes),...
                 'Xdata',tr.x(toshow&~mask&tr.activeNodes))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doThresholdCut(h,varargin) %#ok
% this helper function inactivates nodes based on the threshold cut selected
% with the slider. Entry point: the only way to get into this function is
% by the 'OK' uicontrol next to the slider.
tr = get(gcbf,'userdata');
tr.activeBranches = tr.activeBranches & ~tr.tocollapse(tr.numLeaves+1:tr.numLabels);
% find active nodes by expanding active Branches
activeBranches=find(tr.activeBranches)';
tr.activeNodes(:)=false;
tr.activeNodes(tr.numLabels,1)=true;
for ind = activeBranches(end:-1:1)
    tr.activeNodes(tr.tree(ind,:))=tr.activeNodes(ind+tr.numLeaves);
end
tr.lastThresholdValue = get(tr.slider,'Value');
set([tr.slider,tr.slidertx,tr.sliderok],'Visible','off')
tr.selected(:) = false;
tr.sel2root = path2root(tr,tr.selected);
set(gcbf,'userdata',tr)
changeEditMode(h);
updateTree(gcbf,[],[]);
autoFit(h)
enableAllUI(h)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cancelThresholdCut(h,varargin) %#ok
% this helper function cancels the threshold cut mode and returns to the 
% 'select' mode. Entry point: the ways to get into this function is
% by the 'CANCEL' uicontrol next to the slider (does not exist yet) or by
% mouse click over the axes diring the slider mode.
tr = get(gcbf,'userdata');
set([tr.slider,tr.slidertx,tr.sliderok],'Visible','off')
changeEditMode(h);
updateTree(gcbf,[],[]);
enableAllUI(h)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function expandAll(h,varargin) %#ok
% Callback function to expand all hidden nodes
tr = get(gcbf,'userdata');
x=[0 inf];
[dump,anchor]=min(abs((mean(ylim)-tr.y))+x(1+tr.activeNodes)'); %#ok
tr.activeBranches(:) = true;
tr.activeNodes(:) = true;
set(gcbf,'Userdata',tr);
updateTree(gcbf,[],anchor)
autoFit(h)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback function to save tree
function saveNewick(h,varargin) %#ok
tr = get(gcbf,'userdata');
newtr.tree = tr.tree;
newtr.dist = tr.dist;
newtr.names = tr.names;
phytreewrite(phytree(newtr),'GUI',true);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback function to restore the original tree
function restoreTree(h,varargin) %#ok
tr = getappdata(gcbf,'backupTree');
view(phytree(tr))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback function to load tree
function loadNewick(h,varargin) %#ok
if strcmp(get(h,'Type'),'uimenu') % if caller is the uimenu then needs 
    figtoclose = gcbf;
    tr = phytreeread;             % to pick a file 
    
else                              % if not, caller is the callback from get workspace var
    tr=[];
    pfig = get(h,'Parent');
    if strcmp(get(h,'string'),'Import') || ...
      (strcmp(get(h,'style'),'listbox') && strcmp(get(gcbf,'SelectionType'),'open'))     
        hp = get(pfig,'Userdata'); 
        figtoclose = hp(4);
        ops = get(hp(1),'string');
        if ~isempty(ops)
            tr = evalin('base',ops(get(hp(1),'value'),:));
        end
        close(pfig);
    elseif strcmp(get(h,'string'),'Cancel')
        close(pfig);
    end   
end
if ~isempty(tr)
    propsForFigure = getappdata(figtoclose,'propsForFigure');
    view(tr,[],propsForFigure);
    close(figtoclose)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback function to close tree
function closeNewick(h,varargin) %#ok
close(gcbf)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback function to copy the tree to a figure
function doPublishFigure(h,varargin,AllNodes,callerIsContextMenu) %#ok

pfig = get(h,'Parent');
if strcmp(get(h,'string'),'Cancel') 
    close(pfig); 
    return; 
end

% get options from the publishdlg window and close it
hp = get(pfig,'Userdata'); 
vfig = hp(7);
va = get(hp(1:6),'Value');
va = [va{:}];
switch find(va(1:3))
    case 1; args = {'type','square'};
    case 2; args = {'type','angular'};
    case 3; args = {'type','radial'};
end
args = {args{:},'bra',va(4),'lea',va(5),'ter',va(6)};
close(pfig);

% now select the branch to publish based on selected points
tr = get(vfig,'userdata');
if ~callerIsContextMenu & ~tr.selected %#ok % if called from the uimenu and 
    tr.selected(end) = true;           % nothing is selected pick the root
end
selected = find(tr.selected);
commonpath = true(tr.numLabels,1);
for ind = 1:numel(selected)
    commonpath = commonpath & path2root(tr,selected(ind));
end
branchtoexp = find(commonpath,1);
tr.selected(:) = false;
tr.selected(branchtoexp) = true;
set(vfig,'userdata',tr)
updateTree(vfig,[],[]);

hp = branchtoexp;
%find all nodes to export (i.e. all descendants)
children = false(1,tr.numLabels);
children(hp) = true;
if AllNodes
    for ind = hp-tr.numLeaves:-1:1
        if children(ind+tr.numLeaves)
            children(tr.tree(ind,:))=true;
        end
    end
    permuta = find(children);
else
    for ind = hp-tr.numLeaves:-1:1
        if children(ind+tr.numLeaves)
            children(tr.tree(ind,:))=tr.activeNodes(tr.tree(ind,:));
        end
    end
    braToLea = find(~tr.activeNodes(tr.tree(:,1)))+tr.numLeaves;
    expBran = find(children(tr.numLeaves+1:end)) + tr.numLeaves;
    permuta = [find(children(1:tr.numLeaves)) ...
        intersect(expBran,braToLea)];
    [dump,hs] = sort(tr.y(permuta)); %#ok
    permuta = [permuta(hs) setdiff(expBran,braToLea)];
end

if sum(children)>1 % enough leaves to export ?
    newtr = phytree;
    ipermuta(permuta) = 1:length(permuta);
    numLeaves = (ipermuta(end) + 1)/2;
    newtr.tree = ipermuta(tr.tree(permuta(numLeaves+1:end)-tr.numLeaves,:));
    newtr.dist = tr.dist(permuta);
    newtr.names = tr.names(permuta);
    plot(newtr,true(length(newtr.tree),1),args{:})
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback function for autofit 
function autoFit(h,varargin) %#ok
tr = get(gcbf,'userdata');
set(tr.ha,'Ylim',[min(tr.y(tr.activeNodes))-1,max(tr.y(tr.activeNodes))+1]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback function to Reset view
function myResetView(h,varargin) %#ok
tr = get(gcbf,'userdata');
set(tr.ha,'Ylim',tr.yLim);
set(tr.ha,'Xlim',tr.xLim);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common export function
function exportSubtree(h,varargin,AllNodes,ToWS,callerIsContextMenu) %#ok

tr = get(gcbf,'userdata');
% if called from the uimenu and nothing is selected pick the root
if ~callerIsContextMenu & ~tr.selected %#ok
    tr.selected(end) = true; 
end

selected = find(tr.selected);
commonpath = true(tr.numLabels,1);
for ind = 1:numel(selected)
    commonpath = commonpath & path2root(tr,selected(ind));
end
branchtoexp = find(commonpath,1);
tr.selected(:) = false;
tr.selected(branchtoexp) = true;

set(gcbf,'userdata',tr)
updateTree(gcbf,[],[]);
doExport(branchtoexp,AllNodes,ToWS)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common export function after having selected a point
function doExport(hp,AllNodes,ToWS) 

tr = get(gcbf,'userdata');

%find all nodes to export (i.e. all descendants)
children = false(1,tr.numLabels);
children(hp) = true;
if AllNodes
    for ind = hp-tr.numLeaves:-1:1
        if children(ind+tr.numLeaves)
            children(tr.tree(ind,:))=true;
        end
    end
    permuta = find(children);
else
    for ind = hp-tr.numLeaves:-1:1
        if children(ind+tr.numLeaves)
            children(tr.tree(ind,:))=tr.activeNodes(tr.tree(ind,:));
        end
    end
    braToLea = find(~tr.activeNodes(tr.tree(:,1)))+tr.numLeaves;
    expBran = find(children(tr.numLeaves+1:end)) + tr.numLeaves;
    permuta = [find(children(1:tr.numLeaves)) ...
        intersect(expBran,braToLea)];
    [dump,hs] = sort(tr.y(permuta)); %#ok
    permuta = [permuta(hs) setdiff(expBran,braToLea)];
end
if sum(children)>1 % enough leaves to export ?
    newtr=phytree;
    ipermuta(permuta) = 1:length(permuta);
    numLeaves = (ipermuta(end) + 1)/2;
    newtr.tree = ipermuta(tr.tree(permuta(numLeaves+1:end)-tr.numLeaves,:));
    newtr.dist = tr.dist(permuta);
    newtr.names = tr.names(permuta);
    if ToWS % export to workspace ?
        s = inputdlg('Workspace variable name ?','Export to Workspace',1);
        while ~(isempty(s) || isvarname(s{1}) || isempty(s{1}))
            s = inputdlg('Not a valid variable name, type a MATLAB variable name ?','Export to Workspace',1);
        end
        if ~(isempty(s) || isempty(s{1}))
            assignin('base',s{1},newtr)
        end
    else % no, then export to other viewer
        view(newtr);
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateTree(h,highLight,anchor) %#ok
% Redraws the tree depending on the active Branches
% rather than erase and redraw, we only change specific fields in hlines
% and hdots.

tr = get(h,'userdata');
activeBranches=find(tr.activeBranches)';
oldPos = tr.y(anchor); 
    
% propagate last leaf
lastleaf = 1:tr.numLabels;
for ind = tr.numBranches:-1:1
    if ~tr.activeNodes(tr.tree(ind,1))
        lastleaf(tr.tree(ind,:))=lastleaf(ind+tr.numLeaves);
    end
end

% find x coordinates of branches
tr.x = tr.dist; 
for ind = tr.numBranches:-1:1
    tr.x(tr.tree(ind,:)) = tr.x(tr.tree(ind,:)) + tr.x(ind+tr.numLeaves);
end

% find y coordinates of branches
dummy = lastleaf([true,diff(lastleaf(1:tr.numLeaves))~=0]);
tr.y=zeros(tr.numLabels,1);
tr.y(dummy)=1:length(dummy);
for ind = activeBranches
    tr.y(ind+tr.numLeaves) = mean(tr.y(tr.tree(ind,:)));
end

% update right labels
todis = tr.names(dummy);
set(tr.ha,'ytick',1:length(dummy),'yticklabel',todis)

% show only active branches
set(tr.hlines,'Visible','off')
set(tr.hlines(tr.activeNodes),'Visible','on')

% update coordinates in lines
for ind = 1:tr.numLabels-1
     set(tr.hlines(ind),'Ydata',tr.y([ind,ind,tr.par(ind)]))
     set(tr.hlines(ind),'Xdata',tr.x([ind,tr.par([ind ind])]))
end
set(tr.hlines(tr.numLabels),'Ydata',tr.y(tr.numLabels)*[1 1 1])
set(tr.hlines(tr.numLabels),'Xdata',tr.x(tr.numLabels)*[1 1 1])
          
% update dots
mask = false(tr.numLabels,1); mask(1:tr.numLeaves) = true;
set(tr.hdots(1),'Ydata',tr.y(tr.activeNodes&~mask),'Xdata',tr.x(tr.activeNodes&~mask))
set(tr.hdots(2),'Ydata',tr.y(tr.activeNodes&mask),'Xdata',tr.x(tr.activeNodes&mask))            

% update red dots
set(tr.hseldots(1),'Ydata',tr.y(tr.activeNodes&~mask&tr.selected),'Xdata',tr.x(tr.activeNodes&~mask&tr.selected))
set(tr.hseldots(2),'Ydata',tr.y(tr.activeNodes&mask&tr.selected),'Xdata',tr.x(tr.activeNodes&mask&tr.selected))            

% set the axis holders          
set(tr.axhold,'Ydata',[0.5,max(tr.y(tr.activeNodes))+0.5])
if numel(oldPos)
    set(tr.ha,'ylim',get(tr.ha,'ylim')+tr.y(anchor)-oldPos);
else
    set(tr.ha,'ylim',get(tr.ha,'ylim')) % just touch 'YLim' such that the listener is triggered
end

% turn on indicative modes
tr.indicativeMode = false;
set([tr.datatip tr.hpathline],'visible','off')
set(tr.hlines,'color','black')
set(tr.hldots(1),'Xdata',[],'Ydata',[])
set(tr.hldots(2),'Xdata',[],'Ydata',[])

% save figure data
set(h,'Userdata',tr)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function path2r = path2root(tr,from)
% helper function, finds path to root
path2r = false(tr.numLabels,1);
if (numel(from)~=1 && sum(from)~=1) 
    return; 
end
path2r(from) = true;
temp = find(path2r);
if numel(temp)
    while temp~=tr.numLabels;
        temp = tr.par(temp);
        path2r(temp) = true;
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pathA,pathB] = path2sel(tr,from)
% helper function, finds path to selected node
path2rt = path2root(tr,from);
commonPath = tr.sel2root & path2rt;
commonPath(find(commonPath,1)) = false;
pathB = tr.sel2root & ~commonPath;
pathA = path2rt & ~commonPath;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tr = doBasicCalculations(tr)
% helper function to compute and find some features of the tree
tr = struct(tr);
tr.numBranches = size(tr.tree,1);
tr.numLeaves = tr.numBranches + 1;
tr.numLabels = tr.numBranches + tr.numLeaves; 

% obtain parents for every node
tr.par(tr.tree(:)) = tr.numLeaves + [1:tr.numBranches 1:tr.numBranches];

% calculate the distance to the closest leaf for every node
% needed for fast threshold cut
tr.dist2Leaf = zeros(tr.numLabels,1);
for ind = 1:tr.numBranches
    tr.dist2Leaf(ind+tr.numLeaves) = ...
       min(tr.dist2Leaf(tr.tree(ind,:))+tr.dist(tr.tree(ind,:)));
end

% calculate drawing coordinates for the tree: x coordinated will never
% change, but y coordinates may change depending on the active branches and
% nodes. 
tr.x = tr.dist; tr.y=[1:tr.numLeaves zeros(1,tr.numBranches)]';
for ind = tr.numBranches:-1:1
    tr.x(tr.tree(ind,:)) = tr.x(tr.tree(ind,:)) + tr.x(ind+tr.numLeaves);
end
for ind =1:tr.numBranches
    tr.y(ind+tr.numLeaves) = mean(tr.y(tr.tree(ind,:)));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function correctFigureSize(fig,recommendedHeight)
% helper function to increase initial figure size depending on the screen &
% tree sizes
screenSize = diff(reshape(get(0,'ScreenSize'),2,2),[],2)-100;
            % 100 gives extra space for the figure header and win toolbar
position = get(fig,'Position');
if recommendedHeight > position(4)
    if recommendedHeight < sum(position([2 4]))
        position(2) = sum(position([2 4])) - recommendedHeight;
        position(4) = recommendedHeight;
    elseif recommendedHeight < screenSize(2)
        position(2) = 30; 
        position(4) = recommendedHeight;
    else 
        position(2) = 30; 
        position(4) = screenSize(2);
    end
    set(fig,'Position',position)
end    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hToggleToolbar = makePhyTreeViewerToolbar(fig)
% helper function to set the toolbar
%
% hToggleToolbar contains handles to easy change the state on/off when
%                changing modes

oldSH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

set(fig,'toolbar','figure')  % needs to update because uicontrols turn it off

% Fix toolbar options, we keep: ZoomIn,ZoomOut,Pan
hw = findall(fig,'type','uitoolbar');
hf = get(hw,'Children');
h1 = findall(hf,'Tag','Exploration.Pan');
h2 = findall(hf,'Tag','Exploration.ZoomOut');
h3 = findall(hf,'Tag','Exploration.ZoomIn');
delete(setxor(hf,[h1,h2,h3]))
set([h1 h2 h3],'Separator','off','clickedCallback',@changeEditMode);

% load icons
load(fullfile(matlabroot,'toolbox','bioinfo','bioinfo','@phytree','phytreeicons'))

h4 = uitoggletool('ToolTip','Inspect Tool Mode','separator','on',...
                  'Tag','Inspect',        'CData',icons(:,:,1:3)); %#ok
h5 = uitoggletool('ToolTip','Collapse/Expand Branch Mode',...
                  'Tag','Collapse/Expand','CData',icons(:,:,4:6)); %#ok
h6 = uitoggletool('ToolTip','Rotate Branch Mode',...
                  'Tag','Rotate Branch',  'CData',icons(:,:,7:9)); %#ok
h7 = uitoggletool('ToolTip','Rename Leaf/Branch Mode',...
                  'Tag','Rename',         'CData',icons(:,:,10:12)); %#ok
h8 = uitoggletool('ToolTip','Prune (delete) Leaf/Branch Mode',...
                  'Tag','Prune',          'CData',icons(:,:,13:15)); %#ok
set([h4 h5 h6 h7 h8],'clickedCallback',@changeEditMode,'state','off',...
                     'Serializable','off','HandleVisibility','off');   
hToggleToolbar = [h4 h5 h6 h7 h8 h1 h2 h3];
set(0,'ShowHiddenHandles',oldSH)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hToggleUIMenu = makePhyTreeViewerUIMenus(fig)
% helper function to set UI menus
% 
% hToggleUIMenu contains handles to easy check on/off modes in the UI menu

oldSH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

% delete figure menus not used
h1 = findall(fig,'Type','uimenu', 'Label','&Edit');
h2 = findall(fig,'Type','uimenu', 'Label','&View');
h3 = findall(fig,'Type','uimenu', 'Label','&Insert');
h4 = findall(fig,'Type','uimenu', 'Label','&Desktop');
delete([h1,h2,h3,h4])

% Repair "File" menu
hw = findall(fig,'Type','uimenu', 'Label','&File');
hf = get(hw,'children');
h1 = findall(hw,'Label','Expo&rt Setup...');
h3 = findall(hw,'Label','Print Pre&view...');
h4 = findall(hw,'Label','&Print...');
delete(setxor(hf,[h1,h3,h4]))
uimenu(hw,'Label','New Tool...', 'Position',1,'Callback','phytreetool')
uimenu(hw,'Label','Open...',     'Position',2,'Callback',@loadNewick)
uimenu(hw,'Label','Import from Workspace...','Position',3,'Callback',@importfromwsdlg)
uimenu(hw,'Label','Open Original in New Tool','Position',4,'Callback',@restoreTree)
uimenu(hw,'Label','Save As...',  'Position',5,'Callback',@saveNewick,'Separator','on')
item0 = uimenu(hw ,'Label','Print to Figure','Position',6);
    uimenu(item0,'Label','With Hidden Nodes...','Callback',{@publishdlg,1,0});
    uimenu(item0,'Label','Only Displayed...',   'Callback',{@publishdlg,0,0});
item1 = uimenu(hw,'Label','Export to New Tool','Position',7);
    uimenu(item1,'Label','With Hidden Nodes...','Callback',{@exportSubtree,1,0,0});
    uimenu(item1,'Label','Only Displayed...',   'Callback',{@exportSubtree,0,0,0});
item2 = uimenu(hw,'Label','Export to Workspace','Position',8);
    uimenu(item2,'Label','With Hidden Nodes...','Callback',{@exportSubtree,1,1,0});
    uimenu(item2,'Label','Only Displayed...',   'Callback',{@exportSubtree,0,1,0});
uimenu(hw,'Label','Exit','Separator','on','Position',12,'Callback',@closeNewick)
set(h1,'Separator','on')
    
% Repair "Tools" menu
hw = findall(fig,'Type','uimenu','Label','&Tools');
hf = get(hw,'children');
h1 = findall(hw,'Tag','figMenuZoomIn');    set(h1,'Callback',{@changeEditMode,gcbo});
h2 = findall(hw,'Tag','figMenuZoomOut');   set(h2,'Callback',{@changeEditMode,gcbo});
h3 = findall(hw,'Tag','figMenuPan');       set(h3,'Callback',{@changeEditMode,gcbo});
h4 = findall(hw,'Tag','figMenuResetView'); set(h4,'Callback',{@myResetView,gcbo});
h5 = findall(hw,'Tag','figMenuOptions');
set([h1,h4],'separator','off')
delete(setxor(hf,[h1,h2,h3,h4,h5]))
delete(findall(h5,'Tag','figMenuOptionsDatatip'))
delete(findall(h5,'Tag','figMenuOptionsDataBar'))
h6 = uimenu(hw,'Label','Inspect','Position',1,'Callback',@changeEditMode);
h7 = uimenu(hw,'Label','Collapse/Expand','Position',2,...
                                              'Callback',@changeEditMode);
h8 = uimenu(hw,'Label','Rotate Branch','Position',3,...
                                              'Callback',@changeEditMode);
h9 = uimenu(hw,'Label','Rename','Position',4, 'Callback',@changeEditMode);
h10 = uimenu(hw,'Label','Prune','Position',5, 'Callback',@changeEditMode);
item3 = uimenu(hw,'Label','Threshold Collapse','Position',9,'Separator','on');
   uimenu(item3,'Label','Distance to Leaves',  'Callback',@thresholdCut);
   uimenu(item3,'Label','Distance to Root',    'Callback',@thresholdCut);
   uimenu(hw,'Label','Expand All','Position',10,'Callback',@expandAll);
uimenu(hw,'Label','Find Leaf/Branch...','Position',11,'Callback',@findNode); 
uimenu(hw,'Label','Fit to Window','Position',12,'Separator','on',...
                                                 'Callback',@autoFit);
set(h1,'Separator','on')

% Repair "Help" menu
hw = findall(fig,'Type','uimenu','Label','&Help');
delete(get(hw,'children'));
uimenu(hw,'Label','Bioinformatics Toolbox Help','Position',1,'Callback',...
       'web([docroot ''/toolbox/bioinfo/bioinfo_product_page.html''])')
uimenu(hw,'Label','Phylogenetic Tree Tool Help','Position',2,'Callback',...
       ['helpview(fullfile(docroot,''toolbox'',''bioinfo'', ''bioinfo.map'')' ...
        ',''phytreetool_reference'')' ] )
uimenu(hw,'Label','Demos','Position',3,'Separator','on',...
       'Callback','demo(''toolbox'',''bioinfo'')')   
tlbx = ver('bioinfo');
mailstr = ['web(''mailto:bioinfofeedback@mathworks.com?subject=',...
           'Feedback%20for%20Phytreetool%20in%20Bioinformatics',...
           '%20Toolbox%20',tlbx(1).Version,''')'];
uimenu(hw,'Label','Send Feedback','Position',4,'Separator','on',...
       'Callback',mailstr);
set(0,'ShowHiddenHandles',oldSH)
hToggleUIMenu = [h6 h7 h8 h9 h10 h1 h2 h3];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hToggleContextMenu,hAxisContextMenu,hDotsContextMenu] = ...
         makePhyTreeViewerContextMenus(fig) %#ok
% helper function to set context menus
% 
% hToggleContextMenu contains handles to easy check on/off modes in the
%                    Context Menu
% hAxisContextMenu   contains handle to the Axis Context Menu to used to
%                    reactivate it
% hDotsContextMenu   contains handle to the Dots Context Menu to used to
%                    reactivate it

% set up context menu for the axes (when mouse is not over a node)
hcm1 = uicontextmenu('Callback',@hideActiveIndicators);

h1 = uimenu(hcm1,'Label','Inspect',        'Callback',@changeEditMode);
h2 = uimenu(hcm1,'Label','Collapse/Expand','Callback',@changeEditMode);
h3 = uimenu(hcm1,'Label','Rotate Branch',  'Callback',@changeEditMode);
h4 = uimenu(hcm1,'Label','Rename',         'Callback',@changeEditMode);
h5 = uimenu(hcm1,'Label','Prune',          'Callback',@changeEditMode);
item4 = uimenu(hcm1 , 'Label', 'Threshold Collapse','Separator','on');
   uimenu(item4,'Label','Distance to Leaves',    'Callback',@thresholdCut);
   uimenu(item4,'Label','Distance to Root',      'Callback',@thresholdCut);
   uimenu(hcm1,'Label','Expand All','Callback',@expandAll);
uimenu(hcm1 ,'Label','Find Leaf/Branch...',           'Callback',@findNode);
uimenu(hcm1 ,'Label','Fit to Window','Separator','on','Callback',@autoFit);
uimenu(hcm1 ,'Label','Reset to Original View',        'Callback',@myResetView);

% context menu for dots (when mouse over a dot and right mouse button pressed)
hcm2 = uicontextmenu('Callback',@enableOptionsContextMenu);
   uimenu(hcm2,'Label','Collapse/Expand', 'Callback',@collapseExpand);
   uimenu(hcm2,'Label','Rotate Branch',   'Callback',@rotateBranch);
   uimenu(hcm2,'Label','Rename',          'Callback',@renameNode);
   uimenu(hcm2,'Label','Prune',           'Callback',@pruneTree);
item0 = uimenu(hcm2,'Label','Print to Figure','Separator','on');   
   uimenu(item0,'Label','With Hidden Nodes...','Callback',{@publishdlg,1,1});
   uimenu(item0,'Label','Only Displayed...',   'Callback',{@publishdlg,0,1});
item1 = uimenu(hcm2,'Label','Export to New Tool');
   uimenu(item1,'Label','With Hidden Nodes...','Callback',{@exportSubtree,1,0,1});
   uimenu(item1,'Label','Only Displayed...',   'Callback',{@exportSubtree,0,0,1});
item2 = uimenu(hcm2,'Label','Export to Workspace');
   uimenu(item2,'Label','With Hidden Nodes...','Callback',{@exportSubtree,1,1,1});
   uimenu(item2,'Label','Only Displayed...',   'Callback',{@exportSubtree,0,1,1});
   %uimenu(hcm2,'Label','Color Down',    'Callback',@colorDown);
   
% save context menus in my data structure to later restore them if desired
hToggleContextMenu = [h1 h2 h3 h4 h5];
hAxisContextMenu = hcm1;
hDotsContextMenu = hcm2;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enableOptionsContextMenu(h,varargin)
hideActiveIndicators
selectNode(h,varargin)
tr = get(gcbf,'Userdata');
hc = get(tr.hDotsContextMenu,'Children');
set(hc,'Enable','on')
if sum(tr.selected)~=1 
    set(hc(5),'Enable','off'); 
end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setupYLabelsListeners  
% helper function to setsup listeners for the ylables, so we can detect if
% we would need to change the fontsize
hgp     = findpackage('hg');
axesC   = findclass(hgp,'axes');
figureC = findclass(hgp,'figure');
% listens when the Ylim of axes has changed
YLimListener = handle.listener(gca,axesC.findprop('YLim'),...
               'PropertyPostSet',{@ylabelsListener,gcf,gca});
% listens when Position of Figure has changed
PositionListener = handle.listener(gcf,figureC.findprop('Position'),...
                   'PropertyPostSet',{@ylabelsListener,gcf,gca});
% store the listeners
setappdata(gcf,'PhyTreeListeners',[YLimListener, PositionListener]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function publishdlg(h,varargin,AllNodes,callerIsContextMenu) %#ok
% dialog window to select options for publishing
vfig = gcbf;
c = get(0,'ScreenSize')*[1 0;0 1;.5 0;0 .5];
fig = figure('WindowStyle','modal','Color',[0.831373 0.815686 0.784314],...
             'Position',[c-[150 100] 300 200],'Resize','off','NumberTitle','off',...
             'Name','Print Phylogenetic Tree to Figure','IntegerHandle','off' );
h1=uibuttongroup;h2=uibuttongroup;
set(h1,'Position',[.08 .35 .35 .55],'Title','Rendering Type','backgroundcolor',[0.831373 0.815686 0.784314])
set(h2,'Position',[.52 .35 .39 .55],'Title','Display Labels','backgroundcolor',[0.831373 0.815686 0.784314])
ui1=uicontrol(h1,'style','radiobutton','Position',[5 70 90 20],'string','Square','value',1,'backgroundcolor',[0.831373 0.815686 0.784314]);
ui2=uicontrol(h1,'style','radiobutton','Position',[5 40 90 20],'string','Angular','backgroundcolor',[0.831373 0.815686 0.784314]);
ui3=uicontrol(h1,'style','radiobutton','Position',[5 10 90 20],'string','Radial','backgroundcolor',[0.831373 0.815686 0.784314]);
ui4=uicontrol(h2,'style','checkbox','Position',[5 70 109 20],'string','Branch Nodes','backgroundcolor',[0.831373 0.815686 0.784314]);
ui5=uicontrol(h2,'style','checkbox','Position',[5 40 109 20],'string','Leaf Nodes','backgroundcolor',[0.831373 0.815686 0.784314]);
ui6=uicontrol(h2,'style','checkbox','Position',[5 10 109 20],'string','Terminal Nodes','value',1,'backgroundcolor',[0.831373 0.815686 0.784314]);
uicontrol(fig,'style','pushbutton','Position',[70 20 60 30],'string','Print','Callback',{@doPublishFigure,AllNodes,callerIsContextMenu});
uicontrol(fig,'style','pushbutton','Position',[155 20 60 30],'string','Cancel','Callback',{@doPublishFigure,AllNodes,callerIsContextMenu});
set(fig,'Userdata',[ui1 ui2 ui3 ui4 ui5 ui6 vfig]);
set(h1,'SelectionChangeFcn',{@toggleCheckBoxs,ui3,ui6})

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggleCheckBoxs(h,event,h3,h6) %#ok
if get(h3,'value')
    set(h6,'enable','off')
else
    set(h6,'enable','on')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function importfromwsdlg(h,varargin) %#ok
% dialog window to select variable from workspace
mvars = evalin('base','whos');
mvars = mvars(strmatch('phytree',{mvars(:).class})); %#ok
vfig = gcbf;
c = get(0,'ScreenSize')*[1 0;0 1;.5 0;0 .5];
fig = figure('WindowStyle','modal','Color',[0.831373 0.815686 0.784314],...
             'Position',[c-[80 100] 160 220],'Resize','off','NumberTitle','off',...
             'Name','Get Phytree Object','IntegerHandle','off' );
ui1=uicontrol(fig,'style','list','Position',[22 70 120 120],'string',char(mvars(:).name),'backgroundcolor','w','Callback',@loadNewick);         
ui2=uicontrol(fig,'style','pushbutton','Position',[15 20 60 30],'string','Import','Callback',@loadNewick);
ui3=uicontrol(fig,'style','pushbutton','Position',[90 20 60 30],'string','Cancel','Callback',@loadNewick);         
uicontrol(fig,'style','text','Position',[20 190 140 20],'string','Select phytree object:','Horizontal','left','backgroundcolor',[0.831373 0.815686 0.784314])
set(fig,'Userdata',[ui1 ui2 ui3 vfig]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = getphytreetoolnumber()
% Computes the index number for this particular tool

% first, finds the used numbers so far
allFigs = findall(0,'tag','PhyTreeTool');
usedNumbers = zeros(1,numel(allFigs)+1);
baseName = 'Phylogenetic Tree Tool ';
baseLen = length(baseName);
for i = 1:numel(allFigs)
    str = get(allFigs(i),'Name');
    usedNumbers(i) = str2double(str(baseLen:end));
end

% This is how we find the next index.  The rule is that we find the lowest
% integer value (non-zero and positive) not yet prescribed to a phytree
% tool, This is the same way MATLAB figures behave.
n = num2str(min(setdiff(1:(max(usedNumbers)+1),usedNumbers)));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = getacceptedwarningfromothertools(cfig) %#ok
% Finds out if the pruning warning has been accepted
% first, finds the used numbers so far
allFigs = findall(0,'tag','PhyTreeTool');
n = 'NotDone';
for i = 1:numel(allFigs)
    propsForFigure = getappdata(allFigs(i),'propsForFigure');
    if isequal(propsForFigure.PruneWarning,'Done')
        n = 'Done';
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setacceptedwarningtoothertools()
% Set the pruning warning to all the other open tools as 'Done'
allFigs = findall(0,'tag','PhyTreeTool');
for i = 1:numel(allFigs)
    propsForFigure = getappdata(allFigs(i),'propsForFigure');
    propsForFigure.PruneWarning = 'Done';
    setappdata(allFigs(i),'propsForFigure',propsForFigure);
end
    
