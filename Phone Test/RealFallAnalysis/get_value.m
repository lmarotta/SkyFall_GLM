function value=get_value(text) % returns corresponding label value for text

expression= {'like', 'sli', 'tri','rig', 'lef', 'walk', 'stair', 'up','down', 'sit', 'stand', 'start', 'end'};
lm=cellfun(@(x) regexp(x,expression,'match','ignorecase'),text,'UniformOutput',false);
lm=cellfun(@(x) cell2mat(cellfun(@(y) ~isempty(y),x,'UniformOutput',false)),lm,'UniformOutput',false);
lm=cell2mat(lm);

ind_values=[4 1 2 3 4 9 11 13 15 17 19 0 1]';

value=lm*ind_values;
