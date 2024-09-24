function cb_catch_button_press(varargin)

for k = 1:length(varargin)
varargin{k}
end

currentChar = get(varargin{1},'CurrentCharacter')

size(currentChar)

disp('here!')

end