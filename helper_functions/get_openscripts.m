function ScriptNames = get_openscripts

%% get open scripts names
openFiles = matlab.desktop.editor.getAll;
ScriptNames = {openFiles.Filename};
temp = cellfun(@(x) strsplit(x,'\'),ScriptNames,'uniformoutput',false);
ScriptNames = cellfun(@(x) x(end),temp);
rmindx = cell2mat(cellfun(@(x) strcmp(x,'get_openscripts.m'), ScriptNames,'uniformoutput',false));
ScriptNames(rmindx) = [];
%% save names of open scripts
cd('G:\My Drive\MATLAB\Code\Single firefly with motion cuing')
save('ScriptNames','ScriptNames')