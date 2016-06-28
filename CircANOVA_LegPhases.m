% Circular ANOVA for testing actual vs. predicted leg phases

run TransformLegPhases_Spreadsheet.m
clear c_* lsN* bpN*

legPhaseSpreadsheet_vel = legPhaseSpreadsheet(:,3);
[legPhaseSpreadsheet_SLOW, legPhaseSpreadsheet_MED, legPhaseSpreadsheet_FAST] = SetSpeedBins(legPhaseSpreadsheet,legPhaseSpreadsheet_vel);

% [pval table] = circ_hktest(alpha, idp, idq, inter, fn)

% [pval, stats] = circ_hktest(alpha, idp, idq, inter, fn)
%   Parametric two-way ANOVA for circular data with interations.
%
%   Input:
%     alpha   angles in radians
%     idp     indicates the level of factor 1 (1:p)
%     idq     indicates the level of factor 2 (1:q)
%     inter   0 or 1 - whether to include effect of interaction or not
%     fn      cell array containing strings with the names of the factors
%               
%
%   Output:
%     pval    vector of pvalues testing column, row and interaction effects
%     table   cell array containg the anova table
%
%   The test assumes underlying von-Mises distributrions.
%   All groups are assumed to have a common concentration parameter k,
%   between 0 and 2.

actualLegPhases_s = legPhaseSpreadsheet_SLOW(:,2); % actual leg phases % alpha
LvsREffect_s = legPhaseSpreadsheet_SLOW(:,4); % LvsR % idp
segNumEffect_s = legPhaseSpreadsheet_SLOW(:,5); % segNum %idq

actualLegPhases_m = legPhaseSpreadsheet_MED(:,2); % actual leg phases % alpha
LvsREffect_m = legPhaseSpreadsheet_MED(:,4); % LvsR % idp
segNumEffect_m = legPhaseSpreadsheet_MED(:,5); % segNum %idq

actualLegPhases_f = legPhaseSpreadsheet_FAST(:,2); % actual leg phases % alpha
LvsREffect_f = legPhaseSpreadsheet_FAST(:,4); % LvsR % idp
segNumEffect_f = legPhaseSpreadsheet_FAST(:,5); % segNum %idq

fn = {'LvsR','segNum'};

[pval_s table_s] = circ_hktest(actualLegPhases_s, LvsREffect_s, segNumEffect_s, 0, fn); % ANOVA SLOW

[pval_m table_m] = circ_hktest(actualLegPhases_m, LvsREffect_m, segNumEffect_m, 0, fn); %ANOVA MED

[pval_f table_f] = circ_hktest(actualLegPhases_f, LvsREffect_f, segNumEffect_f, 0, fn); %ANOVA FAST





