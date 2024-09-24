octave-signal-ORIGINAL
 - contains the original signal package for GNU OCTAVE
octave-signal
 - contains modified code to work with MATLAB


=====================================
SYNTAX DIFFERENCES THAT CAN BE AUTOMATICALLY CHANGED USING "octave2matlab.bash"
find and replace
OCTAVE  	MATLAB
#		%
!		~
endif		end
endswitch	end
endfor		end
endfunction	end


**some functions begin with two underscores, which is INVALID in MATLAB.