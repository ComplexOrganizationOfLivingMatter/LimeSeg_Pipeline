function [output_txt] = pickCell(obj,event_obj)
%PICKCELL Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');

labelledImage = getappdata(0, 'labelledImageTemp_Resized');
selectedCell = labelledImage(pos(1), pos(2), getappdata(0, 'selectedZ'));
output_txt = {['ID Cell: ',num2str(selectedCell)]};

myhandles = guidata(get(event_obj, 'Target'));
%myhandles.tbCellId.String = num2str(selectedCell);
set(myhandles.tbCellId, 'String', num2str(selectedCell));

% hObject = get(event_obj, 'Target');
% 
% notify(@(hObject,eventdata)window('tbCellId_Callback',hObject,[],guidata(hObject)))
end

