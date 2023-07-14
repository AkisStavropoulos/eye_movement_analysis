function openscripts(varargin)
%% Open scripts when you want to run or edit something

if isempty(varargin)
    
    % Open pre-determined scripts    
    open FinalSensoryAnalysis
    open export_lfp_struct_EDOARDO
%     open PreliminaryAnalysis
%     open LookAllSubjects
%     open SensoryPaperAnalysis
%     open AnalyseEyemovement
    % open EyeMovementsAnalysis
%     open AccuracyAnalysis
%     open JerkAnalysis
%     open ConditionalVarianceAnalysis
%     open StrategyAnalysis
%     open CoasterBrakerCriteria
    
else
    
    % Open scripts that you saved from Editor
   temp = cat(2,varargin{:});
   cellfun(@(x) open(x),temp);

end