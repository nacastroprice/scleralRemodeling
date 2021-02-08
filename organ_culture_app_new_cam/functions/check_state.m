function [state,currentDate] = check_state(startDate,endDate)
% Checks the current time and changes state
    % update currentDate
    currentDate = datetime('now');
%     dateTaken = datetime('now');
    % Check if state has changed and change if has
    if (currentDate < startDate)
        state = 1; % before it is supposed to start
    elseif ((currentDate >= startDate) && (currentDate <= endDate))
        state = 2; % while it is supposed to take pictures
    elseif (currentDate > endDate)
        state = 3; % after it is supposed to have ended
    end
end

