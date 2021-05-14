function MAC = getMAC(Ph_1,Ph_2)
%
% MAC = getMAC(Ph1,Ph2)
%
% INPUTS:
%   Ph1: mode shapes 1, each column is one mode shape.
%   Ph2: mode shapes 2, each column is one mode shape.
%
% OUTPUTS:
%   MAC:
%
% %%%%%%%%%%%%%%%%%%
% %%% Mateo G.H. %%%
% %%% 2021/05/12 %%%
% %%%%%%%%%%%%%%%%%%
tic_MAC = tic;
%%% -----------------------------------------------------------------------
NumCoordinates_1 = size(Ph_1,1);
NumCoordinates_2 = size(Ph_2,1);
%%% -----------------------------------------------------------------------
if NumCoordinates_1 == NumCoordinates_2
    %%% -------------------------------------------------------------------
    NumModeShapes_1  = size(Ph_1,2);
    NumModeShapes_2  = size(Ph_2,2);
    MAC = zeros(NumModeShapes_1,NumModeShapes_2);
    %%% -------------------------------------------------------------------
    for COUNTER_Ph_1 = 1:NumModeShapes_1
        PhByCounter_1 = Ph_1(:,COUNTER_Ph_1);
        %%% ---------------------------------------------------------------
        for COUNTER_Ph_2 = 1:NumModeShapes_2
            PhByCounter_2 = Ph_2(:,COUNTER_Ph_2);
            %%%
            MAC(COUNTER_Ph_1,COUNTER_Ph_2) = ...
                (abs(PhByCounter_1'*PhByCounter_2))^2/((PhByCounter_1'*PhByCounter_1)*(PhByCounter_2'*PhByCounter_2));
        end
        %%% ---------------------------------------------------------------
    end
    %%% -------------------------------------------------------------------
else
    error('Error: El número de filas de Ph_1 y Ph_2 no coincide. La función getMAC no se ha ejecutado.')
end
%%% -----------------------------------------------------------------------
disp(['MAC: ',num2str(toc(tic_MAC),'%.2f')])
end